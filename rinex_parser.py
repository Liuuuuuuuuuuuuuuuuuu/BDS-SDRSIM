import numpy as np
import re
from datetime import datetime

class BDSRinexParser:
    def __init__(self):
        self.data = []

    def parse_file(self, file_path):
        """
        Parses a RINEX Navigation file (v3.02+ preferred for BDS)
        and extracts ephemeris parameters for BeiDou satellites (Cxx).
        """
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: File {file_path} not found.")
            return []

        header_end = False
        current_sat = None
        
        # Buffer to hold lines for a single record
        record_buffer = []
        
        print(f"Parsing {file_path}...")
        
        for line in lines:
            if "END OF HEADER" in line:
                header_end = True
                continue
            
            if not header_end:
                continue

            # RINEX 3.0+ starts records with the Satellite ID (e.g., C01)
            # Regex to find lines starting with 'C' followed by digits
            if line.startswith('C'):
                # If we have a previous record buffered, process it
                if record_buffer:
                    self._process_record(record_buffer)
                    record_buffer = []
                
                record_buffer.append(line)
            
            # Continuation lines for the current record (usually start with space)
            elif line.startswith(' ') and record_buffer:
                record_buffer.append(line)

        # Process the final record
        if record_buffer:
            self._process_record(record_buffer)

        print(f"Successfully parsed {len(self.data)} ephemeris blocks.")
        return np.array(self.data)

    def _process_record(self, lines):
        """
        Extracts Keplerian parameters from a single RINEX record block.
        Standard RINEX 3.04 Navigation File Format for Beidou/Galileo/GPS.
        """
        # Line 0: SatID, Time, SV Clock Bias, Drift, Drift Rate
        # C01 2023 01 01 00 00 00  .1234e-04  .1234e-11  .0000e+00
        l1 = lines[0].replace('D', 'e').replace('d', 'e') # Handle Fortran scientific notation
        
        try:
            # Parse Epoch Time
            year = int(l1[4:8])
            month = int(l1[9:11])
            day = int(l1[12:14])
            hour = int(l1[15:17])
            minute = int(l1[18:20])
            sec = int(l1[21:23])
            
            # Clock params
            clock_bias = float(l1[23:42])
            clock_drift = float(l1[42:61])
            clock_drift_rate = float(l1[61:80])
            
            # We need to flatten the rest of the parameters into a standard feature vector
            # The specific layout depends on RINEX version, but standard usually follows:
            
            # Line 1: IODE, Crs, Delta N, M0
            l2 = lines[1].replace('D', 'e')
            iode = float(l2[4:23])
            crs = float(l2[23:42])
            delta_n = float(l2[42:61])
            m0 = float(l2[61:80]) # Mean Anomaly
            
            # Line 2: Cuc, e, Cus, sqrt(A)
            l3 = lines[2].replace('D', 'e')
            cuc = float(l3[4:23])
            e = float(l3[23:42])  # Eccentricity
            cus = float(l3[42:61])
            sqrt_a = float(l3[61:80]) # Sqrt Semi-major Axis
            
            # Line 3: Toe, Cic, Omega0, Cis
            l4 = lines[3].replace('D', 'e')
            toe = float(l4[4:23]) # Time of Ephemeris
            cic = float(l4[23:42])
            omega0 = float(l4[42:61]) # Longitude of Ascending Node
            cis = float(l4[61:80])
            
            # Line 4: i0, Crc, omega, Omega_dot
            l5 = lines[4].replace('D', 'e')
            i0 = float(l5[4:23])  # Inclination
            crc = float(l5[23:42])
            omega = float(l5[42:61]) # Argument of Perigee
            omega_dot = float(l5[61:80]) # Rate of Right Ascension
            
            # Line 5: IDOT, L2 Codes, Week, L2 P flag
            l6 = lines[5].replace('D', 'e')
            idot = float(l6[4:23]) # Rate of Inclination
            
            # Feature Vector for AI (15 params usually)
            # [ClockBias, ClockDrift, ClockDriftRate, Crs, DeltaN, M0, Cuc, e, Cus, SqrtA, Toe, Cic, Omega0, Cis, i0, Crc, Omega, OmegaDot, IDOT]
            features = [
                clock_bias, clock_drift, clock_drift_rate,
                crs, delta_n, m0,
                cuc, e, cus, sqrt_a,
                toe, cic, omega0, cis,
                i0, crc, omega, omega_dot,
                idot
            ]
            
            self.data.append(features)
            
        except (ValueError, IndexError) as e:
            print(f"Skipping malformed line: {e}")
            pass

if __name__ == "__main__":
    # Test with a dummy file if run directly
    parser = BDSRinexParser()
    print("RINEX Parser initialized. Import this class to use it.")
