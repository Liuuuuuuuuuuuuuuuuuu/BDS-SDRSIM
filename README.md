# BDS-SDRSIM

This project generates synthetic BeiDou B1I baseband samples for SDR
experiments.  Navigation data is read from a RINEX navigation file and
converted into the binary subframe format required by receivers.
（中文）本程式可為 SDR 實驗產生北斗 B1I 基帶訊號，並從 RINEX 檔轉換導航數據成接收機需要的子幀格式。

## System Overview

BDS‑SDRSIM parses BeiDou ephemeris from the RINEX navigation file,
computes satellite positions and Doppler shifts for the requested start
time and user location, and builds the B1I navigation subframes.  Each
enabled satellite channel spreads these bits with the appropriate PRN
code.  The 50 bps D1 navigation message is further modulated by the
standard 20‑bit Neumann–Hoffman sequence so that the resulting signal
matches the BeiDou B1I specification.
Finally the channels are summed to
produce complex baseband samples ready for SDR playback.
The output is ready to be transmitted by an SDR.
（中文）BDS‑SDRSIM 會從 RINEX 檔讀取星曆並計算衛星位置與多普勒，產生 B1I 子幀並展碼，最後合成複數基帶取樣供 SDR 播放。

## Detailed System Description

BDS-SDRSIM reads BeiDou ephemeris from a RINEX navigation file and
initialises a set of satellite channels. For each simulation step the
tool computes the user position (static or from a path file), derives
satellite geometry and Doppler, and updates the channel state. The
navigation message for subframes 1–5 is generated with correct
time-of-week. Each channel spreads the bits with its PRN code and the
samples are summed into a 16‑bit I/Q stream at a fixed 5.120 MHz sample
rate.
（中文）程式讀取 RINEX 星曆後初始化各衛星通道，依序計算使用者與衛星幾何與多普勒，並產生對應的導航訊息，在 5.120u2009MHz 取樣率下輸出 16 位元 I/Q 流。

## Build

```
make
```

The optional command `make prn_test` builds a small utility that prints
the first chips of PRN codes using the bundled RINEX file.
（中文）可用 `make prn_test` 編譯小工具，列印隨附 RINEX 檔中 PRN 碼的起始片段。

The build produces `bds-sim` which outputs a signed `beidou_b1i.bin`
file containing interleaved **16‑bit little‑endian** I/Q samples. By
default the file is written at 5.120 MHz.  Each sample is a pair of
16‑bit integers `(I,Q)` so be
careful not to interpret the file as 8‑bit data—otherwise the Q channel
may appear to contain only zeros or `-1` values.
Using the `--byte` flag writes `beidou_b1i_u8.bin` containing the **8‑bit**
I channel only for quick inspection.  The Q samples are discarded
entirely.
（中文）`bds-sim` 會產生 `beidou_b1i.bin`，內容為交錯的 16 位元小端 I/Q 取樣，預設取樣率 5.120u2009MHz；若誤以 8 位元讀取將看到 Q 通道僅為 0 或 -1。加上 `--byte` 可輸出 8 位元 I 取樣，舊旗標 `-byte` 仍可使用。
The legacy option `-byte` is still recognised as an alias for `--byte`.

### Signal Types
GEO PRN (1–5, 59–63) → D2 (500 bps, no NH)
GEO 衛星（PRN 1–5, 59–63）使用 D2（500 bps，無 Neumann–Hoffman 二次碼）
MEO/IGSO PRN (6–58) → D1 (50 bps, with NH)
MEO/IGSO 衛星（PRN 6–58）使用 D1（50 bps，包含 Neumann–Hoffman 二次碼）
--force-d2 forces D2 on every PRN.
--force-d2 旗標可將所有 PRN 強制改為 D2。

## Usage

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --duration 60 \
          --gain 1.0 \
          --llh lat,lon,height \
          --byte
```

（中文）以上範例顯示基本用法，可輸入 RINEX 檔名、起始時間與其他參數。
`--gain` scales the output amplitude.  With the default gain of `1.0`
the composite signal uses most of the 16‑bit range.  Larger values boost
the level but may cause clipping.
（中文）`--gain` 用來調整輸出振幅，預設 1.0 已足夠，過大可能造成飽和。

The amplitude itself is derived from a simple link budget.  Each orbit
type is assigned a nominal transmit power (about 52 dBm for GEO,
53 dBm for IGSO and 55 dBm for MEO).  Path loss is computed from the
current slant range including a fixed 2 dB atmospheric term.  The gain
factor multiplies this result before limiting to the 16‑bit output
range.  The computed power is further adjusted toward the target
C/N₀ value (45 dB‑Hz by default) so that adding more channels does not
（中文）振幅依簡易鏈路預算計算，包含各軌道標稱功率、路徑損耗及 2dB 大氣衰減，再依目標 C/N₀ 調整，避免通道數多時功率下降。
reduce signal strength excessively.

`--noise` adds complex AWGN with the given standard deviation (in 16‑bit
sample units).  The default is `0` (no noise).
（中文）`--noise` 設定 AWGN 標準差，預設 0；`--seed` 控制噪聲與初相位的隨機性，預設 1。
`--seed` specifies the random seed for both noise generation and the
（中文）`--byte` 只輸出 I 分量的 8 位元檔案。
initial carrier phase of each channel. The default is `1`.

`--byte`  saves an 8‑bit file containing only the I samples.

`--start` specifies the UTC start time; `--duration` is in seconds and
`--llh` defines the user location in degrees and meters. The start time may
extend up to 24 hours past the last ephemeris in the RINEX file; the simulator
continues using the nearest ephemeris block for interpolation. These static
（中文）`--start` 設定 UTC 起始時間；`--duration` 為秒數；`--llh` 指定使用者位置。亦可用動態路徑檔 `--xyz`、`--llh-file` 或 `--nmea`。
coordinates are mutually exclusive with the dynamic path options
`--xyz`, `--llh-file` and `--nmea`.
Run `./bds-sim --help` to see all command options. A brief configuration
summary is printed before signal generation begins.

### Dynamic user trajectory

Instead of a fixed location you may supply a 1 Hz path file. Three formats are supported:
（中文）若非固定位置，可提供 1Hz 路徑檔，支援三種格式。

```
--xyz       ECEF coordinates (x y z in metres)
--llh-file  Geodetic coordinates (lat lon h)
--nmea      NMEA GGA sentences
```
（中文）以上三種格式分別為 ECEF XYZ、經緯度高程，以及 NMEA GGA。

Example path files are provided in the `examples/` directory. Usage:
The files contain one position per line at a 1 Hz rate. Coordinates are
（中文）範例路徑檔位於 `examples/` 目錄，每行代表 1Hz 位置，可為 XYZ、LLH 或 NMEA GGA。
either ECEF XYZ in metres, latitude/longitude/height in degrees and
metres, or NMEA GGA sentences.

```
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --xyz examples/path_xyz.txt \
          --duration 60
```

（中文）此範例展示使用路徑檔模擬動態使用者。
## Code layout

- `navbits.c`  – builds subframes 1–5 from ephemeris data.
- `channel.c`  – generates PRN code and navigation modulation for each
  satellite channel.
- `bdssim.c`   – combines all channels into the final sample stream.
- `rinex.c`    – small parser for RINEX navigation files.
（中文）程式碼主要分布於以上檔案。

Subframes 2 and 3 are pre‑built from ephemeris while subframes 1, 4 and
5 are generated on demand so that the Time‑of‑Week is consistent with
（中文）其中子幀 2 與 3 先由星曆建立，其餘則根據模擬時間動態產生。
the simulation start.

## SDR playback

The file `beidou_b1i.bin` contains interleaved **16‑bit little‑endian**
I/Q samples written at a fixed 5.120 MHz sample rate.
Ensure any analysis or playback software reads the file as 16‑bit
integers; using 8‑bit interpretation will yield
Q samples that look like zeros.
When `--byte` is used, the output file `beidou_b1i_u8.bin` stores only the
I samples in signed 8‑bit format.  Only the real component is retained;
Q is not saved.
（中文）輸出的 `beidou_b1i.bin` 為交錯 16 位元 I/Q 取樣，播放軟體需以 16 位元載入；使用 `--byte` 時僅保存 8 位元 I 取樣。

The signal is at baseband (zero‑IF), so tune the
SDR’s RF centre frequency to the desired transmit frequency – for B1I
this is typically **1561.098 MHz** – and play the samples at the same
rate.  The following command illustrates playback with a HackRF:

```bash
hackrf_transfer -t beidou_b1i.bin -f 1561098000 -s 5120000 -x 0
（中文）信號為零中頻，請將 SDR 中心頻率設為 1561.098MHz 並以同樣取樣率播放；以下範例使用 HackRF。
```

Use the `-R` option for continuous looping if needed.  Other SDRs can
（中文）如需持續循環可加上 -R；其他 SDR 亦可在支援相同速率下傳輸。
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
（中文）解析 `beidou_b1i_u8.bin` 時，請確保 GNSS-SDR 的 PRN 設定與產生的衛星一致，或開啟自動搜尋。

## v0.3.0 (2025-07-07)
* Fix UTC→BDT +14 s
* GEO satellites enabled
* Added subframe 4/5 template
* Default Fs → 5.120 MHz
* Power scaling by target CN₀
* Basic D2 (500 bps) support

## v0.3.1 (unreleased)
* Correct subframe 4/5 constant words
* Minor API cleanup

