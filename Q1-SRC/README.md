# Multirate Audio Sample‑Rate Converter (SRC)

## Authors:  Jonathan Akuaku Edem, Felicia Archeresuah, Grace Wiredu

## Overview
- This project includes MATLAB examples showing sample‑rate conversion 96 kHz → 44.1 kHz using upsample → FIR → downsample.
- The script files are: `single-stage.m` (single‑stage SRC) and `two-stage.m` (two‑stage SRC).

Quick run
1. Open MATLAB and set the current folder to the repository root.
2. Run `single-stage.m` or `two-stage.m`.

Outputs
- Interactive stem plots for input, upsampled, filtered and downsampled signals.
- Console diagnostics: designed filter lengths and resulting sampling rates.

Requirements
- MATLAB R2020a or newer.
- Signal Processing Toolbox (for `designfilt`, optional `fvtool`).

Troubleshooting (short)
- Missing toolbox functions → install Signal Processing Toolbox.
- Plots/files not produced → check write permissions and run interactively (GUI required for `fvtool`).

