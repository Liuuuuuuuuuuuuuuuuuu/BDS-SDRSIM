## v0.3.0 (2025-07-07)
* Fix UTC→BDT +14 s
* GEO satellites enabled
* Added subframe 4/5 template
* Default Fs → 8.192 MHz
* Power scaling by target CN₀
* Basic D2 (500 bps) support with D1 interleaving

## v0.3.1 (unreleased)
* Subframe 4/5 words set to official 0xAAAAAAAA pattern
* Remove unused GEO helper
* Amplitude API no longer depends on channel count
* Allow --start up to 24 h past last ephemeris
* Improved utc_to_bdt() handling of negative offsets
* Orbit classification now distinguishes IGSO and MEO by
  checking the ephemeris `sqrtA` value
* Orbit propagation now uses the unified GPS/BDS RAAN
  expression `Ω(t)=Ω₀+(Ω̇−Ωₑ)·tk−Ωₑ·Toe` for all BeiDou
  satellites
* Added `--meo-only` flag to simulate only MEO satellites
* Skip satellites marked unhealthy in ephemeris and
  always emit `SatH1=0` in navigation data

