# BDS-SDRSIM

This project generates synthetic BeiDou B1I baseband samples for SDR
experiments.  Navigation data is read from a RINEX navigation file and
converted into the binary subframe format required by receivers.

## Build

```
make
```

The build produces `bds-sim` which outputs a signed `beidou_b1i.bin`
file containing interleaved I/Q samples at 8.184 MHz.

## Usage

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --duration 60 \
          --llh lat,lon,height
```

`--start` specifies the UTC start time; `--duration` is in seconds and
`--llh` defines the user location in degrees and meters.

### Dynamic user trajectory

Instead of a fixed location you may supply a 1 Hz path file. Three formats are supported:

```
--xyz       ECEF coordinates (x y z in metres)
--llh-file  Geodetic coordinates (lat lon h)
--nmea      NMEA GGA sentences
```

Example path files are provided in the `examples/` directory. Usage:

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
