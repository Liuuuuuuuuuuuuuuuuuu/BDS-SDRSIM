<p align="center">
  <a href="https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset=".github/assets/bds-sdrsim-card-dark.png">
      <img src=".github/assets/bds-sdrsim-card-light.png" alt="BDS-SDRSIM" width="880">
    </picture>
  </a>
</p>

# BDS-SDRSIM

BDS-SDRSIM is an open-source C simulator that synthesizes BeiDou B1I baseband signals for software-defined radio experiments.（BDS-SDRSIM 是一個以 C 語言撰寫的開源模擬器，可產生北斗 B1I 基帶訊號供 SDR 測試）

## Features
- **Signal output**: Generates interleaved I/Q samples in 16-bit (5.2 MHz) or 8-bit (25 MHz) format.（輸出 16 位元或 8 位元的 I/Q 樣本）
- **Navigation data**: Builds D1 subframes from RINEX ephemeris and inserts 20‑bit Neumann–Hoffman codes.（從 RINEX 星曆產生 D1 子幀並加入 20 位元 Neumann–Hoffman 序列）
- **Multi-satellite channels**: Automatically selects visible PRNs or restricts to user-specified sets.（自動選取可視衛星或依參數限制 PRN）
- **User motion**: Supports fixed LLH coordinates or 1 Hz trajectory files (`--xyz`, `--llh-file`, `--nmea`).（可模擬靜止或動態使用者）
- **Power control**: Simple link budget model with adjustable gain and target C/N₀.（內建簡易鏈結損失模型，可調整增益與 C/N₀）
- **No external dependencies**: Portable C code compiled with `make`.（純 C 實作，透過 `make` 即可編譯）

## Architecture and Signal Flow（整體架構與流程）
1. **Ephemeris parsing** – `rinex.c` reads broadcast parameters from RINEX files.（由 `rinex.c` 解析 RINEX 星曆）
2. **Satellite state** – `orbits.c` solves Kepler’s equation  
   `E - e\sin E = M` to obtain the ECEF position and velocity. The range  
   is `r = A (1 - e\cos E)` and the RAAN is  
   `Ω(t) = Ω₀ + (Ω̇ - Ωₑ)·tk - Ωₑ·Toe`.（`orbits.c` 使用 Kepler 方程計算衛星軌道及 RAAN）
3. **Receiver frame** – `coord.c` converts between LLA and ECEF using WGS‑84:  
   `N = a / \sqrt{1 - e² \sin²φ}`,  
   `x = (N + h) \cosφ \cosλ`, etc.（`coord.c` 進行 WGS‑84 座標轉換）
4. **Geometry** – `bdssim.c` iterates to compute pseudorange `ρ` and range rate `ṙ`,
   applying Earth rotation `Ωₑτ` for Sagnac correction:
   `ρ = ‖r_sat - r_usr‖ - c·Δt_sv`.（`bdssim.c` 求得幾何距離與地球自轉修正）
5. **Ionospheric delay** – `iono.c` applies the Klobuchar model to correct `ρ`.（`iono.c` 以 Klobuchar 模型修正電離層延遲）
6. **Channel model** – `channel.c` derives carrier and code dynamics:
   `f_d = -f_B1I·ṙ/c`,
   `f_code = f_chip·(1 - ṙ/c)`. Amplitude uses a link budget
   `A ∝ 10^{(C/N₀ - 45)/20} · (1/ρ) · G_ant · gain`.（`channel.c` 依多普勒與功率模型產生通道）
7. **Navigation bits** – `navbits.c` assembles five 6‑s subframes and applies BCH encoding.（`navbits.c` 建構 6 秒的導航子幀並進行 BCH 編碼）
8. **Sample synthesis** – `bdssim.c` mixes PRN, NH code, navigation bits and carrier to produce interleaved I/Q output.（`bdssim.c` 混合 PRN、NH 以及導航資料產生 I/Q 樣本）

## Quick Start（快速開始）
1. **Build**
   ```bash
   make
   ```
   Produces executable `bds-sim`.（編譯生成 `bds-sim`）
2. **Fetch ephemeris**
   ```bash
   make download-brdm2025176
   ```
   Downloads MGEX BRDM file.（下載 MGEX 星曆）
3. **Run simulation**
   ```bash
   ./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
             --start 2025/06/25,00:00:00 \
             --duration 60 \
             --llh 25.0,121.0,30 \
             --gain 1.0
   ```
   Generates `beidou_b1i.bin` (16‑bit). Add `--byte` for `beidou_b1i_u8.bin`.（預設輸出 16 位元，加入 `--byte` 轉為 8 位元）
4. **Playback**
   ```bash
   hackrf_transfer -t beidou_b1i.bin -f 1561098000 -s 5200000 -x 0
   ```
   Sends the samples with HackRF.（以 HackRF 播放樣本）

## Common Options（常用參數）

| Option | Description |
|--------|-------------|
| `--rinex file` | RINEX navigation file (required).（指定 RINEX 星曆） |
| `--start YYYY/MM/DD,hh:mm:ss` | Simulation start time UTC.（模擬起始 UTC 時刻） |
| `--duration sec` | Simulation length 1–3600 s.（模擬秒數） |
| `--llh lat,lon,h` | Fixed user position in degrees/meters.（固定使用者位置） |
| `--xyz/--llh-file/--nmea file` | 1 Hz trajectory file for moving user.（匯入路徑檔模擬移動） |
| `--gain amp` | Output amplitude gain (>0).（輸出振幅倍率） |
| `-cn0 value` | Target C/N₀ to balance channel power.（目標 C/N₀） |
| `--byte` | Output 8‑bit I/Q at 25 MHz.（以 8 位元輸出 25 MHz 樣本） |
| `--meo-only` | Use only MEO satellites.（僅模擬 MEO 衛星） |
| `--prn N` | Simulate only PRN N.（指定單一 PRN） |
| `--prn37` | Restrict to PRN 1–37.（限制 PRN 1–37） |
| `--no-iono` | Disable ionospheric delay model.（停用電離層延遲模型） |
| `--seed n` | Random seed.（設定亂數種子） |

All options: `./bds-sim --help`.（完整參數請見 `./bds-sim --help`）

## Dynamic Path Example（動態路徑範例）
```bash
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --xyz examples/path_xyz.txt \
          --duration 60
```
The user position updates every second from `examples/path_xyz.txt`.（使用者位置將依檔案每秒更新）

## Testing（測試）
Build and run unit tests:
```bash
make check
```
Executes `tests/test_prn` and `tests/test_nh_prn`.（執行 PRN 與 NH 序列測試）

Additional test resources:
- `tests/gnss_sdr/` – GNSS-SDR configuration example and multiple analysis reports.（GNSS-SDR 使用設定範例與多個分析報告）
- `tests/gnuradio_usrp/` – GNU Radio USRP settings and a smartphone test using the Android “GPS Test” app.（GNU Radio USRP 設定與 Android「GPS Test」應用程式的測試結果）

## Project Structure（專案結構）
- `bdssim.c` – Main simulation loop and I/Q writer.（核心模擬流程與輸出）
- `channel.c` – Per‑satellite signal generation and power control.（通道計算與功率模型）
- `orbits.c` – Satellite orbit and clock model.（衛星軌道與時鐘模型）
- `navbits.c` – B1I subframe assembly.（導航資料組合）
- `rinex.c` – RINEX navigation parser.（星曆解析）
- `coord.c`, `path.c` – Coordinate transforms and user trajectory.（座標轉換與路徑處理）
- `iono.c` – Klobuchar ionospheric delay correction.（Klobuchar 電離層延遲修正）
- `examples/` – Sample trajectory files.（範例路徑）
- `tests/` – Unit tests.（測試程式）

## License
Released under the MIT License. See [LICENSE](LICENSE).（以 MIT 授權發布）
