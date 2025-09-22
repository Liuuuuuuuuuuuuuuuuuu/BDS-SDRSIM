<p align="center">
  <a href="https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="assets/bds-sdrsim-card-dark.png?v=1">
      <source media="(prefers-color-scheme: light)" srcset="assets/bds-sdrsim-card-light.png?v=1">
      <img alt="BDS-SDRSIM" src="assets/bds-sdrsim-card-light.png?v=1" width="880">
    </picture>
  </a>
</p>

<p align="center">
  <a href="README.en.md">English</a> ·
  <a href="README.zh-TW.md">繁體中文</a>
</p>

---

# BDS-SDRSIM（繁體中文）

## 簡介
BDS-SDRSIM 是以 C 語言撰寫的北斗 B1I 基帶模擬器，可輸出交錯的 I/Q 樣本供 SDR 播放或離線分析。系統支援 MGEX/IGS RINEX 星曆，模擬衛星動態並產生含 Neumann–Hoffman 與 BCH 編碼的導航資料。

## 特色
- **完整導航鏈**：支援 RINEX 3.0x NAV 解析、D1 子幀組裝、NH20 與 BCH 30-bit 字。
- **多通道訊號模型**：可同時模擬 8 顆衛星，含多普勒、碼速、Klobuchar 電離層延遲與可調 C/N₀。
- **彈性使用者軌跡**：提供固定 LLH、1 Hz ECEF/LLH 路徑檔與 NMEA GGA 重播。
- **輸出可調**：預設 5.2 Msps 16 位元，啟用 `--byte` 後輸出 25 Msps 8 位元基帶。
- **純 C99 編譯**：僅需 libm 與 OpenMP，易於跨平台建置。

## 系統架構
1. `rinex.c` 解析星曆與電離層參數。
2. `orbits.c` 計算衛星位置／速度與時鐘修正。
3. `coord.c` 完成 LLH、ECEF、ENU 座標轉換。
4. `bdssim.c` 管理衛星選取、通道更新與 I/Q 合成。
5. `channel.c` 產生每毫秒 PRN、NH、導航位元與載波。
6. `navbits.c` 將星曆轉為北斗 D1 子幀並進行 BCH/NH 編碼。
7. `iono.c` 計算 Klobuchar 電離層延遲。
8. `path.c` 載入及內插使用者路徑。

## 需求環境
- 具 `make`、`gcc`/`clang`、`curl`、`gunzip` 的 POSIX 類系統。
- 支援 OpenMP (`-fopenmp`) 的編譯器。macOS 請先安裝 `libomp`。
- 若要透過 CDDIS 自動下載 MGEX 檔案，需事先於 `cookies.txt` 儲存登入資訊。

## 環境建置範例

### Ubuntu / Debian / Linux Mint（含 WSL）
```bash
sudo apt update
sudo apt install -y build-essential git curl gzip libomp-dev
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make
make check                     # 可選：執行自我檢查
make download-brdm2025176      # 下載示範 MGEX 星曆
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \\
          --start 2025/06/25,00:00:00 \\
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

### Fedora / RHEL / Rocky Linux
```bash
sudo dnf groupinstall -y "Development Tools"
sudo dnf install -y git curl gzip libomp
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make CC=gcc
make check
make download-brdm2025176
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \\
          --start 2025/06/25,00:00:00 \\
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

### Arch Linux / Manjaro
```bash
sudo pacman -S --needed base-devel git curl gzip libomp
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make CC=gcc
make check
make download-brdm2025176
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \\
          --start 2025/06/25,00:00:00 \\
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

### macOS（Homebrew 工具鏈）
```bash
xcode-select --install                           # 第一次使用時安裝開發工具
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"  # 尚未安裝 Homebrew 時使用
brew install gcc libomp git curl
export CC=$(brew --prefix)/bin/$(ls $(brew --prefix)/bin | grep -E '^gcc-[0-9]+$' | sort -V | tail -1)
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make CC="$CC"
make check
make download-brdm2025176
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \\
          --start 2025/06/25,00:00:00 \\
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

> 若偏好使用 Clang，可先安裝 `llvm` 與 `libomp`，再執行 `make CC=$(brew --prefix llvm)/bin/clang CFLAGS="-fopenmp" LDFLAGS="-lomp"`。

### Windows 11 / 10（WSL2）
```powershell
wsl --install                                       # 預設啟用 Ubuntu，若提示請重新開機
wsl
sudo apt update
sudo apt install -y build-essential git curl gzip libomp-dev
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make
make check
make download-brdm2025176
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \\
          --start 2025/06/25,00:00:00 \\
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

## 建置
```bash
make            # 產生 ./bds-sim 執行檔
make clean      # 清除目標檔、可執行檔與暫存檔
```

## 下載星曆
```bash
make download-brdm2025176   # 下載 BRDM00DLR_S_20251760000_01D_MN.rnx.gz
```
`download-%` 規則可將 `brdmYYYYDDD` 或 `brdcDDD0.YYn` 擴展為完整檔名。若因權限失敗，請手動下載並解壓後放在專案根目錄。

## 快速開始
```bash
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --llh 25.0,121.0,30 \
          --duration 60 \
          --gain 1.0
```
- 預設輸出 `beidou_b1i.bin`（16 位元 I/Q，5.2 Msps）。
- 加入 `--byte` 可輸出 `beidou_b1i_u8.bin`（8 位元 I/Q，25 Msps，0 Hz IF）。

### SDR 播放範例
```bash
hackrf_transfer -t beidou_b1i.bin -f 1561098000 -s 5200000 -x 0
```
請依需求調整 TX 增益，以達到預期的 C/N₀。

## 指令參考
| 參數 | 說明 |
| --- | --- |
| `--rinex file` | 指定 RINEX 星曆檔（必填）。 |
| `--start YYYY/MM/DD,hh:mm:ss` | 模擬起始 UTC 時刻。 |
| `--duration sec` | 模擬秒數（1–3600 秒）。 |
| `--llh lat,lon,h` | 固定使用者位置（度、度、公尺）。 |
| `--xyz file` | 1 Hz ECEF 軌跡檔。 |
| `--llh-file file` | 1 Hz LLH 軌跡檔（度、度、公尺）。 |
| `--nmea file` | NMEA GGA 紀錄檔（1 Hz）。 |
| `--gain amp` | 線性輸出增益 (>0)。 |
| `--fs MHz` | 取樣率（MHz）。預設 5.2（16 位元）或 25（搭配 `--byte`）。 |
| `--byte` | 8 位元交錯輸出，25 Msps，0 Hz IF。 |
| `--seed n` | 隨機載波初相位的種子值。 |
| `--meo-only` | 僅模擬 MEO 衛星。 |
| `--prn N` | 指定單一 PRN。 |
| `--prn37` | 將星座限制為 PRN 1–37。 |
| `--no-iono` | 停用電離層延遲模型。 |
| `-cn0 value` | 舊版 C/N₀ 覆寫（套用於全部通道）。 |

> `--llh` 與路徑檔參數二擇一，請勿同時指定。`--fs` 以 MHz 為單位，程式會自動換算為 Hz。

## 使用者軌跡
- `examples/path_xyz.txt` – 1 Hz ECEF 座標（公尺）。
- `examples/path_llh.txt` – 1 Hz 緯度、經度、海拔（度、公尺）。
- `examples/path_nmea.nmea` – 範例 NMEA GGA 記錄。
程式會在線性內插兩筆樣點之間的座標，產生子秒級更新。

## 輸出檔案
- `beidou_b1i.bin` – 16 位元小端序 I/Q。
- `beidou_b1i_u8.bin` – 8 位元有號 I/Q，適用於 HackRF/USRP 播放。
終端機會顯示模擬時間進度。

## 測試
```bash
make check               # 執行 PRN / NH / Iono 單元測試
```
額外資源：
- `tests/gnss_sdr/` – GNSS-SDR 設定檔與分析報告。
- `tests/gnuradio_usrp/` – GNU Radio 流程圖與場測紀錄。

## 專案結構
- `main.c` – 命令列參數解析與檢查。
- `bdssim.c` – 模擬主程式與輸出。
- `channel.c` – 單衛星通道模型與樣本生成。
- `orbits.c` – 星曆軌道與時鐘推算。
- `navbits.c` – 北斗 D1 子幀組合與 BCH/NH 編碼。
- `rinex.c` – RINEX NAV 解析器。
- `coord.c`、`path.c` – 座標轉換與軌跡內插。
- `iono.c` – Klobuchar 電離層延遲模型。
- `examples/`、`tests/`、`assets/` – 範例資料、測試套件、圖像資源。

## 疑難排解
- **OpenMP 無法使用**：請安裝 `libomp` 或改用 GCC。若必須停用 OpenMP，可移除 `-fopenmp` 並註解 `bdssim.c` 中的 `#pragma omp`。
- **「start outside nav window」**：調整 `--start` 至距離 TOE ±1 天內。
- **輸出功率不足**：調整 `--gain`、`-cn0` 或 SDR 發射增益。

## 授權條款
詳見 [LICENSE](LICENSE.md)。
