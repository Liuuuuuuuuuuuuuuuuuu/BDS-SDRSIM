# BDS-SDRSIM

This project generates synthetic BeiDou B1I baseband samples for SDR
experiments.  Navigation data is read from a RINEX navigation file and
converted into the binary subframe format required by receivers.

## System Overview

BDS‑SDRSIM parses BeiDou ephemeris from the RINEX navigation file,
computes satellite positions and Doppler shifts for the requested start
time and user location, and builds the B1I navigation subframes.  Each
enabled satellite channel spreads these bits with the appropriate PRN
code.  The 50 bps D1 navigation message is further modulated by the
standard 20‑bit Neumann–Hoffman sequence so that the resulting signal
matches the BeiDou B1I specification.  Finally the channels are summed to
produce complex baseband samples ready for SDR playback.
The output is ready to be transmitted by an SDR.

## Detailed System Description

BDS-SDRSIM reads BeiDou ephemeris from a RINEX navigation file and
initialises a set of satellite channels. For each simulation step the
tool computes the user position (static or from a path file), derives
satellite geometry and Doppler, and updates the channel state. The
navigation message for subframes 1–5 is generated with correct
time-of-week. Each channel spreads the bits with its PRN code and the
samples are summed into a 16‑bit I/Q stream at 5 MHz.

## Build

```
make
```

The optional command `make prn_test` builds a small utility that prints
the first chips of PRN codes using the bundled RINEX file.

The build produces `bds-sim` which outputs a signed `beidou_b1i.bin`
file containing interleaved **16‑bit little‑endian** I/Q samples at
5 MHz.  Each sample is a pair of 16‑bit integers `(I,Q)` so be
careful not to interpret the file as 8‑bit data—otherwise the Q channel
may appear to contain only zeros or `-1` values.

## Usage

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --duration 60 \
          --gain 1.0 \
          --llh lat,lon,height
```

`--gain` scales the output amplitude.  With the default gain of `1.0`
the composite signal uses most of the 16‑bit range.  Larger values boost
the level but may cause clipping.

`--start` specifies the UTC start time; `--duration` is in seconds and
`--llh` defines the user location in degrees and meters.
Run `./bds-sim --help` to see all command options. A brief configuration
summary is printed before signal generation begins.

### Dynamic user trajectory

Instead of a fixed location you may supply a 1 Hz path file. Three formats are supported:

```
--xyz       ECEF coordinates (x y z in metres)
--llh-file  Geodetic coordinates (lat lon h)
--nmea      NMEA GGA sentences
```

Example path files are provided in the `examples/` directory. Usage:
The files contain one position per line at a 1 Hz rate. Coordinates are
either ECEF XYZ in metres, latitude/longitude/height in degrees and
metres, or NMEA GGA sentences.

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --xyz examples/path_xyz.txt \
          --duration 60
```

## Code layout

- `navbits.c`  – builds subframes 1–5 from ephemeris data.
- `channel.c`  – generates PRN code and navigation modulation for each
  satellite channel.
- `bdssim.c`   – combines all channels into the final sample stream.
- `rinex.c`    – small parser for RINEX navigation files.

Subframes 2 and 3 are pre‑built from ephemeris while subframes 1, 4 and
5 are generated on demand so that the Time‑of‑Week is consistent with
the simulation start.

## SDR playback

The file `beidou_b1i.bin` contains interleaved **16‑bit little‑endian**
I/Q samples at 5 MHz.  Ensure any analysis or playback software
reads the file as 16‑bit integers; using 8‑bit interpretation will yield
Q samples that look like zeros.

The signal is at baseband (zero‑IF), so tune the
SDR’s RF centre frequency to the desired transmit frequency – for B1I
this is typically **1561.098 MHz** – and play the samples at the same
rate.  The following command illustrates playback with a HackRF:

```bash
hackrf_transfer -t beidou_b1i.bin -f 1561098000 -s 5000000 -x 0
```

Use the `-R` option for continuous looping if needed.  Other SDRs can
transmit the file in a similar manner as long as they support the same
sample rate.
