# BDS-SDRSIM

This project generates synthetic BeiDou B1I baseband samples for SDR
experiments.  Navigation data is read from a RINEX navigation file and
converted into the binary subframe format required by receivers.

## Project Summary
BDS-SDRSIM converts RINEX ephemeris into BeiDou B1I navigation frames and outputs complex baseband samples for software-defined radios. It supports both D1 and D2 messages, static or dynamic user motion and configurable power with optional AWGN. See `docs/frame_layout.md` for an overview of the D1/D2 structure.


## System Overview

BDS‑SDRSIM parses BeiDou ephemeris from the RINEX navigation file,
computes satellite positions and Doppler shifts for the requested start
time and user location, and builds the B1I navigation subframes.  Each
enabled satellite channel spreads these bits with the appropriate PRN
code.  The 50 bps D1 navigation message is further modulated by the
standard 20‑bit Neumann–Hoffman sequence so that the resulting signal
matches the BeiDou B1I specification.  Optionally a 500 bps D2 message
without secondary coding can be interleaved every other millisecond.
Finally the channels are summed to
produce complex baseband samples ready for SDR playback.
The output is ready to be transmitted by an SDR.

## 功能總覽

程式以 `main.c` 處理命令列與基本參數，並由 `bdssim.c` 統整整個
模擬流程。`rinex.c` 讀取星曆，`navbits.c` 依時間產生子幀位元，
`channel.c` 依各衛星計算擴頻與都卜勒，最後由 `bdssim.c` 將所有
通道的 I/Q 樣本累加成輸出檔。座標轉換與使用者路徑處理分別
由 `coord.c` 與 `path.c` 負責，其間的資料透過 `ephemeris_t` 與
`channel_t` 結構傳遞。

## Detailed System Description

BDS-SDRSIM reads BeiDou ephemeris from a RINEX navigation file and
initialises a set of satellite channels. For each simulation step the
tool computes the user position (static or from a path file), derives
satellite geometry and Doppler, and updates the channel state. The
navigation message for subframes 1–5 is generated with correct
time-of-week. Each channel spreads the bits with its PRN code and the
samples are summed into a 16‑bit I/Q stream at a fixed 6.144 MHz sample
rate.

## Build

```
make
```

The optional command `make tests/prn_test` builds a small utility that prints
each PRN number and the first 16 bits of its code in binary using the
bundled RINEX file.

Run `make check` to build and execute the self test
(`tests/test_prn_d1d2`).

The build produces `bds-sim` which outputs a signed `beidou_b1i.bin`
file containing interleaved **16‑bit little‑endian** I/Q samples. By
default the file is written at **6.144 MHz**.  Each sample is a pair of
16‑bit integers `(I,Q)` so be
careful not to interpret the file as 8‑bit data—otherwise the Q channel
may appear to contain only zeros or `-1` values.
Using the `--byte` flag writes `beidou_b1i_u8.bin` containing the **8‑bit**
I channel only, sampled at **25 MHz** for use with GNSS‑SDR.  The samples are
generated directly from the internal accumulator rather than converting the
16‑bit stream.  All output values are normalised to remain within ±1 before
quantisation so the integers never exceed the format limits.
The legacy option `-byte` is still recognised as an alias for `--byte`.

## 訊號型別

- GEO PRN → D2 (500 bps, 無 NH 二次碼)
- IGSO/MEO PRN → D1 (50 bps, 有 NH 二次碼)
- `--geo-first` 在挑選模擬衛星時會優先加入可見的 GEO 衛星
- `--no-geo` 排除所有 GEO 衛星
- `--prn`  可只輸出指定 PRN

MEO 與 IGSO 依據星曆中的 `sqrtA`(軌道半長軸平方根) 分辨。大於 6000 的視
為與 GEO 相同的軌道高度，即 IGSO，否則為 MEO。

## Usage

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --duration 60 \
          --gain 1.0 \
          --llh lat,lon,height \
          --byte
```

Example prioritising visible GEO satellites:

```
./bds-sim --rinex BRDM00DLR_S_20251870000_01D_MN.rnx \
          --geo-first --duration 120
```

Example excluding GEO satellites:

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --no-geo --duration 60
```

Example generating only PRN 6:

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --prn 6 --duration 60
```

`--gain` scales the output amplitude.  With the default gain of `1.0`
the composite signal uses most of the 16‑bit range.  Larger values boost
the level but may cause clipping.

The amplitude itself is derived from a simple link budget.  Each orbit
type is assigned a nominal transmit power (about 52 dBm for GEO,
53 dBm for IGSO and 55 dBm for MEO).  Orbit type is detected by
checking `sqrtA` so that IGSO and MEO receive the proper power.
Path loss is computed from the
current slant range including a fixed 2 dB atmospheric term.  The gain
factor multiplies this result before limiting to the 16‑bit output
range.  The computed power is further adjusted toward the target
C/N₀ value (45 dB‑Hz by default) so that adding more channels does not
reduce signal strength excessively.

`--noise` adds complex AWGN with the given standard deviation (in 16‑bit
sample units).  The default is `0` (no noise).
`--seed` specifies the random seed for both noise generation and the
initial carrier phase of each channel. The default is `1`.

`--byte`  saves an 8‑bit file containing only the I samples.

`--start` specifies the UTC start time; `--duration` is in seconds and
`--llh` defines the user location in degrees and meters. The start time may
extend up to 24 hours past the last ephemeris in the RINEX file; the simulator
continues using the nearest ephemeris block for interpolation. These static
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

### Testing

```
make check
```

This builds and runs `tests/test_prn_d1d2` to verify D1/D2 generation.

## SDR playback

The file `beidou_b1i.bin` contains interleaved **16‑bit little‑endian**
I/Q samples written at a fixed **6.144 MHz** sample rate.
Ensure any analysis or playback software reads the file as 16‑bit
integers; using 8‑bit interpretation will yield
Q samples that look like zeros.
When `--byte` is used, the output file `beidou_b1i_u8.bin` stores only the
I samples in signed 8‑bit format at **25 MHz**.  Only the real component is retained;
Q is not saved.

The signal is at baseband (zero‑IF), so tune the
SDR’s RF centre frequency to the desired transmit frequency – for B1I
this is typically **1561.098 MHz** – and play the samples at the same
rate.  The following command illustrates playback with a HackRF:

```bash
hackrf_transfer -t beidou_b1i.bin -f 1561098000 -s 5120000 -x 0
```

Use the `-R` option for continuous looping if needed.  Other SDRs can
transmit the file in a similar manner as long as they support the same
sample rate.

## Using GNSS-SDR

When analyzing `beidou_b1i_u8.bin` with GNSS-SDR, ensure that the PRNs
configured in your `gnss-sdr.conf` match the satellites generated by
`bds-sim`. The standard `gnss-sdr_BDS_B1I_byte.conf` file is preset for
PRNs 1–17 and may not align with the PRN list printed by the simulator.
Adjust the `ChannelX.satellite` lines accordingly or set
`Channels.in_acquisition` to the total channel count so that GNSS-SDR
searches for the correct PRNs automatically.

## v0.3.0 (2025-07-07)
* Fix UTC→BDT +4 s
* GEO satellites enabled
* Added subframe 4/5 template
* Default Fs → 6.144 MHz
* Power scaling by target CN₀
* Basic D2 (500 bps) support with D1 interleaving

## v0.3.1 (unreleased)
* Correct subframe 4/5 constant words
* Fix Time-of-Week handling for D2 frames
* Minor API cleanup
* Orbit propagation now applies the standard GPS/BDS RAAN
  expression `Ω(t) = Ω₀ + (Ω̇ − Ωₑ)·tk − Ωₑ·Toe` for all BeiDou
  satellites. Results were cross-checked against the provided SP3
  precise orbit file.

## License

BDS-SDRSIM is released under the MIT License. See the [LICENSE](LICENSE)
file for details.

