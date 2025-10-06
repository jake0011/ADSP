# Advanced Digital Signal Processing (ADSP) 

This repository contains the solutions and code implementations for three practical Digital Signal Processing (DSP) design problems, as part of an academic project in Advanced DSP.

## Overview

This project addresses three DSP tasks:

1. **Two-Stage Polyphase Sample-Rate Converter (SRC) 96 kHz to 44.1 kHz**  
   Design and implement an efficient two-stage polyphase sample-rate converter meeting strict filter specifications and computational constraints for deployment on an ARM-A55 processor.

2. **Fx-LMS Adaptive Active Noise Cancellation (ANC) Headset**  
   Develop an adaptive noise-cancelling system based on the Fx-LMS algorithm to suppress non-stationary road noise, with objectives for adaptation speed and minimal speech distortion, targeting fixed-point DSP deployment.

3. **64-QAM OFDM Baseband Receiver Chain**  
   Simulate a baseband OFDM physical layer receiver including timing synchronization, carrier frequency offset estimation, pilot-based channel estimation (LS and MMSE), equalization, and performance evaluation under multipath and mobility effects.

Each section includes a design rationale, algorithmic details, simulation results, and implementation considerations. Code scripts and plots supporting the work are included in the repository.

## Contents

- `Q1_SRC/` - MATLAB scripts for single and two-stage polyphase sample-rate conversion, technical memo, filter design scripts, and filter response plots.
- `Q2_ANC/` - Python Jupyter Notebook implementing the Fx-LMS ANC algorithm, simulation results including attenuation plots and spectrograms, plus a technical memo.
- `Q3_OFDM/` - MATLAB scripts for OFDM baseband receiver simulation, BER and channel estimation performance plots, and a technical memo.

## Features

- Iterative filter design using Kaiser and Parks-McClellan methods for SRC with numerical verification.
- Polyphase decomposition for computational efficiency analysis and practical MAC counts.
- Fx-LMS algorithm with variable step-size normalization and noise environment simulation.
- Channel estimation comparison with LS and MMSE methods, including mobility and Doppler effects.
- Fixed-point implementation guidance including coefficient quantization and accumulator sizing.
- Clear mapping from question specifications to engineering design choices.

## Usage

- MATLAB is required to run the SRC and OFDM simulation scripts.
- Python with NumPy and Jupyter Notebook support is needed for the ANC simulations.

## Deliverables

- Technical reports summarizing design approaches, results, and considerations.
- Source code for simulation and evaluation.
  


