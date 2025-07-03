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
samples are summed into a 16‑bit I/Q stream at a configurable sample
rate (default 5 MHz).

## Build

```
make
```

The optional command `make prn_test` builds a small utility that prints
the first chips of PRN codes using the bundled RINEX file.

The build produces `bds-sim` which outputs a signed `beidou_b1i.bin`
file containing interleaved **16‑bit little‑endian** I/Q samples. By
default the file is written at 5 MHz.  Each sample is a pair of
16‑bit integers `(I,Q)` so be
careful not to interpret the file as 8‑bit data—otherwise the Q channel
may appear to contain only zeros or `-1` values.
Using the `--byte` flag sets the sample rate to **25 MHz** and additionally
writes `beidou_b1i_u8.bin` with signed **8‑bit** samples for quick
inspection.
The legacy option `-byte` is still recognised as an alias for `--byte`.

## Usage

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --duration 60 \
          --srate 5000000 \
          --gain 1.0 \
          --llh lat,lon,height \
          --byte
```

`--gain` scales the output amplitude.  With the default gain of `1.0`
the composite signal uses most of the 16‑bit range.  Larger values boost
the level but may cause clipping.

The amplitude itself is derived from a simple link budget.  Each orbit
type is assigned a nominal transmit power (about 52 dBm for GEO,
53 dBm for IGSO and 55 dBm for MEO).  Path loss is computed from the
current slant range including a fixed 2 dB atmospheric term.  The gain
factor multiplies this result before limiting to the 16‑bit output
range.

`--noise` adds complex AWGN with the given standard deviation (in 16‑bit
sample units).  The default is `0` (no noise).
`--seed` specifies the random seed for noise generation. The default is `1`.

`--srate` sets the I/Q sample rate in Hertz. The default is `5000000`.
`--byte`  forces a 25 MHz sample rate and also saves an 8‑bit version of
the baseband samples.

`--start` specifies the UTC start time; `--duration` is in seconds and
`--llh` defines the user location in degrees and meters. These static
coordinates are mutually exclusive with the dynamic path options
`--xyz`, `--llh-file` and `--nmea`.
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
I/Q samples written at the configured sample rate (default 5 MHz).
Ensure any analysis or playback software reads the file as 16‑bit
integers; using 8‑bit interpretation will yield
Q samples that look like zeros.
When `--byte` is used, the sample rate is fixed to **25 MHz** and the
companion file `beidou_b1i_u8.bin` stores the same samples in signed
8‑bit format.

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
