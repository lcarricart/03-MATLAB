# LAB 2: Signal Decomposition and Fourier Series

## Overview
This lab contains solutions to signal decomposition and Fourier series problems, organized into separate folders.

## Folder Structure

### Problem1/
Signal decomposition into even and odd parts
- `decompose_signal.m` - Main decomposition function
- `test_decomposition.m` - Test script with plots
- `problem_1a_solution.txt` - Manual solution
- `README.md` - Detailed documentation

### Problem2/
Fourier series of periodic square wave
- `fourier_square_wave.m` - Computes Fourier coefficients
- `plot_fourier_series.m` - Interactive script with user input
- `problem_2a_derivation.txt` - Mathematical derivation
- `README.md` - Detailed documentation

## Quick Start

### Problem 1 (Signal Decomposition)
```matlab
cd Problem1
test_decomposition
```

### Problem 2 (Fourier Series)
```matlab
cd Problem2
plot_fourier_series  % Opens dialog to enter N, f0, and time vector
```

## Key Concepts

### Problem 1: Even/Odd Decomposition
- x_e(t) = 0.5[x(t) + x(-t)] - Even part
- x_o(t) = 0.5[x(t) - x(-t)] - Odd part
- x(t) = x_e(t) + x_o(t) - Reconstruction

### Problem 2: Fourier Series
- Square wave: x(t) toggles between -1 and +1
- a₀ = 0, aₖ = 0 for all k
- bₖ = 4/(πk) for odd k only
- Gibbs phenomenon at discontinuities
