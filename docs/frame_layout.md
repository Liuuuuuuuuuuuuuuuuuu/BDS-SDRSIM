# BeiDou B1I D1 Frame Layout

This document summarises the BeiDou B1I navigation message according to the official standard. The simulator implements a subset of these fields when building the D1 message.

## Overview

The navigation message uses BCH(15,11,1) coding to produce 30‑bit words arranged into subframes and frames. The D1 message uses a Neumann–Hoffman secondary code.

## D1 Frame Layout

```
Superframe (12 min)
 └─ 24 Frames (30 s each)
     └─ 5 Subframes (300 b)
         └─ 10 Words (30 b)
```

Subframe 1–3 carry each satellite's ephemeris and clock parameters and repeat every 30 s. Subframes 4 and 5 contain almanac pages, health data and UTC parameters. The page number (PNUM) for subframe 4/5 is located at bits 44–50. The simulator cycles this field from 1–24 while keeping the remaining payload words filled with the official `0x2AAAAA` placeholder pattern.

## Index Conversions

URA, UDREI, RURAI and GIVEI fields map index values to physical error estimates. For example URAI index 0 corresponds to a user range accuracy of 2 m (1 σ), while 15 indicates no prediction.
