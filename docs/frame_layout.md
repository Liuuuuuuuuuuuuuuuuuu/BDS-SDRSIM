# BeiDou B1I D1/D2 Frame Layout

This document summarises the BeiDou B1I navigation message according to the official standard. The simulator implements a subset of these fields when building the D1 and D2 messages.

## Overview

Both message types share the same error‑correction and word structure. BCH(15,11,1) coding produces 30‑bit words, arranged into subframes and frames. The D1 message uses a Neumann–Hoffman secondary code, while the D2 message is transmitted without it.

| Parameter          | D1                     | D2                 |
|--------------------|-----------------------|--------------------|
| Data rate          | 50 bps (NH 1 kbps)    | 500 bps            |
| Subframe length    | 300 bits (6.0 s)      | 300 bits (0.6 s)   |
| Frame length       | 1500 bits (30 s)      | 1500 bits (3.0 s)  |
| Superframe length  | 36000 bits (12 min)   | 180000 bits (6 min)|

## D1 Frame Layout

```
Superframe (12 min)
 └─ 24 Frames (30 s each)
     └─ 5 Subframes (300 b)
         └─ 10 Words (30 b)
```

Subframe 1–3 carry each satellite's ephemeris and clock parameters and repeat every 30 s. Subframes 4 and 5 contain almanac pages, health data and UTC parameters. The page number (PNUM) for subframe 4/5 is located at bits 44–50.

## D2 Frame Layout

```
Superframe (6 min)
 └─ 120 Frames (3.0 s each)
     └─ 5 Subframes (0.6 s)
         └─ 10 Words (30 b)
```

GEO satellites transmit D2 messages. Subframe assignments are:

1. Pages 1–10 – basic navigation (PNUM1 field, 4 bits)
2. Pages 1–6  – integrity / differential (PNUM2, 4 bits)
3. Pages 1–6  – continuation of differential data
4. Pages 1–6  – reserved
5. Pages 1–120 – almanac, time offsets and ionospheric grid

### Example Page Map (Subframe 5)

Pages 1–13 contain ionospheric grid data, pages 37–60 carry almanac information, and pages 101–102 provide time offsets relative to GPS, Galileo and GLONASS.

## Index Conversions

URA, UDREI, RURAI and GIVEI fields map index values to physical error estimates. For example URAI index 0 corresponds to a user range accuracy of 2 m (1 σ), while 15 indicates no prediction.

