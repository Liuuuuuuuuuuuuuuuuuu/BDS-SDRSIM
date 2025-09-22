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

# BDS-SDRSIM (English)

## Overview
BDS-SDRSIM is an open-source BeiDou B1I baseband simulator written in portable C. It synthesizes interleaved I/Q samples that can be replayed with software-defined radios (SDR) or recorded for offline analysis. The simulator ingests MGEX/IGS RINEX navigation files, models satellite dynamics, and emits navigation data with Neumann–Hoffman and BCH encoding.

## Key Features
- **Full navigation pipeline**: Parses RINEX 3.0x NAV files, assembles D1 subframes, applies NH20 and BCH 30-bit words.
- **Multi-channel signal model**: Up to 8 simultaneous satellites with Doppler, code rate dynamics, ionospheric delay (Klobuchar), and adjustable C/N₀ budgeting.
- **Flexible receiver motion**: Static LLH, 1 Hz ECEF or LLH trajectories, or recorded NMEA GGA tracks.
- **Configurable output**: 16-bit I/Q at 5.2 Msps or 8-bit I/Q at 25 Msps baseband when `--byte` is enabled.
- **Self-contained build**: C99 + OpenMP, no external libraries required beyond libm and optional libomp.

## Architecture
1. `rinex.c` parses broadcast ephemeris and ionospheric parameters.
2. `orbits.c` propagates satellite position/velocity and clock correction.
3. `coord.c` converts between LLH, ECEF, and ENU frames.
4. `bdssim.c` schedules satellite selection, updates channel dynamics, and mixes I/Q samples.
5. `channel.c` generates 1 ms chunks of PRN code, NH sequence, navigation bits, and carrier.
6. `navbits.c` converts ephemerides to BeiDou D1 frames with BCH/NH coding.
7. `iono.c` evaluates Klobuchar ionospheric delay.
8. `path.c` loads/interpolates user trajectories.

## Requirements
- POSIX-like environment with `make`, `gcc`/`clang`, `curl`, and `gunzip`.
- OpenMP support (`-fopenmp`). macOS users should install `libomp` via Homebrew.
- CDDIS credentials if downloading MGEX files directly (`cookies.txt`).

## Platform Setup Examples

### Ubuntu / Debian / Linux Mint (includes WSL)
```bash
sudo apt update
sudo apt install -y build-essential git curl gzip libomp-dev
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make
make check                     # optional self-test
make download-brdm2025176      # sample MGEX nav file
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
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
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
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
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

### macOS (Homebrew toolchain)
```bash
xcode-select --install                           # one-time, provides developer tools
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"  # install Homebrew if needed
brew install gcc libomp git curl
export CC=$(brew --prefix)/bin/$(ls $(brew --prefix)/bin | grep -E '^gcc-[0-9]+$' | sort -V | tail -1)
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make CC="$CC"
make check
make download-brdm2025176
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

> If you prefer Clang, install `llvm` and `libomp`, then invoke `make CC=$(brew --prefix llvm)/bin/clang CFLAGS="-fopenmp" LDFLAGS="-lomp"`.

### Windows 11 / 10 (WSL2)
```powershell
wsl --install                                       # enables Ubuntu by default, reboot if prompted
wsl
sudo apt update
sudo apt install -y build-essential git curl gzip libomp-dev
git clone https://github.com/Liuuuuuuuuuuuuuuuuuu/BDS-SDRSIM.git
cd BDS-SDRSIM
make
make check
make download-brdm2025176
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --llh 25.0,121.0,30 --duration 60 --gain 1.0
```

## Build
```bash
make            # builds ./bds-sim
make clean      # removes objects, binaries, and temporary files
```

## Fetch Ephemeris
```bash
make download-brdm2025176   # downloads BRDM00DLR_S_20251760000_01D_MN.rnx.gz
```
The `download-%` rule expands shorthand such as `brdmYYYYDDD` or `brdcDDD0.YYn`. If automated download fails (e.g., authentication), manually place the decompressed RINEX file in the project root.

## Quick Start
```bash
./bds-sim --rinex BRDM00DLR_S_20251760000_01D_MN.rnx \
          --start 2025/06/25,00:00:00 \
          --llh 25.0,121.0,30 \
          --duration 60 \
          --gain 1.0
```
- Outputs `beidou_b1i.bin` (16-bit I/Q, 5.2 Msps).
- Add `--byte` to emit `beidou_b1i_u8.bin` (8-bit I/Q, 25 Msps at 0 Hz IF).

### SDR Playback Example
```bash
hackrf_transfer -t beidou_b1i.bin -f 1561098000 -s 5200000 -x 0
```
Adjust TX gain to achieve the desired C/N₀ at the victim receiver.

## Command Reference
| Option | Description |
| --- | --- |
| `--rinex file` | RINEX navigation file (mandatory). |
| `--start YYYY/MM/DD,hh:mm:ss` | UTC start time. |
| `--duration sec` | Simulation length (1–3600 s). |
| `--llh lat,lon,h` | Fixed user position (deg, deg, m). |
| `--xyz file` | 1 Hz ECEF trajectory file. |
| `--llh-file file` | 1 Hz LLH trajectory file (deg, deg, m). |
| `--nmea file` | NMEA GGA log (1 Hz). |
| `--gain amp` | Linear amplitude scaling (>0). |
| `--fs MHz` | Output sample rate. Defaults to 5.2 (16-bit) or 25 (`--byte`). |
| `--byte` | 8-bit interleaved output at 25 Msps, 0 Hz IF. |
| `--seed n` | RNG seed for random carrier phase. |
| `--meo-only` | Restrict to MEO satellites. |
| `--prn N` | Simulate only PRN *N*. |
| `--prn37` | Limit constellation to PRN 1–37. |
| `--no-iono` | Disable ionospheric delay model. |
| `-cn0 value` | Legacy C/N₀ override (applies to all channels). |

> Provide either `--llh` or one trajectory option, not both. `--fs` is specified in MHz and converted internally to Hz.

## Receiver Trajectories
- `examples/path_xyz.txt` – Cartesian ECEF (meters) sampled at 1 Hz.
- `examples/path_llh.txt` – Latitude/longitude (degrees) and altitude (meters) at 1 Hz.
- `examples/path_nmea.nmea` – Sample NMEA GGA trace.
The simulator interpolates between samples for sub-second updates.

## Output Files
- `beidou_b1i.bin` – 16-bit little-endian I/Q pairs.
- `beidou_b1i_u8.bin` – 8-bit signed I/Q pairs suitable for HackRF/USRP playback.
Console progress reports elapsed simulation time.

## Testing
```bash
make check               # runs PRN/NH/Iono unit tests
```
Additional resources:
- `tests/gnss_sdr/` – GNSS-SDR configurations and analysis reports.
- `tests/gnuradio_usrp/` – GNU Radio flowgraphs and field-test notes.

## Directory Layout
- `main.c` – CLI parsing and configuration validation.
- `bdssim.c` – Simulation driver and sample writer.
- `channel.c` – Per-satellite channel model and sample synthesis.
- `orbits.c` – Orbit/clock propagation from ephemerides.
- `navbits.c` – BeiDou D1 subframe generation.
- `rinex.c` – RINEX NAV parser.
- `coord.c`, `path.c` – Coordinate transforms and trajectory interpolation.
- `iono.c` – Klobuchar ionospheric delay.
- `examples/`, `tests/`, `assets/` – Sample data, validation suites, artwork.

## Troubleshooting
- **OpenMP unsupported**: Install `libomp` or rebuild with GCC. If OpenMP must be disabled, remove `-fopenmp` and the `#pragma omp` block in `bdssim.c`.
- **"start outside nav window"**: Choose a start time within ±1 day of the ephemeris TOE.
- **Low playback power**: Increase `--gain`, adjust `-cn0`, or raise SDR transmit gain.

## License
See [LICENSE](LICENSE.md).
