# Problem 2: Fourier Series of Square Wave

## Overview
This folder contains MATLAB scripts to compute and plot the Fourier series of a periodic square wave that toggles between -1 and +1.

## Theory
For a square wave x(t) toggling between -1 and +1 with period T₀ = 1/f₀:

**Fourier Series:**
```
x(t) = a₀ + Σ[aₖ·cos(2πkf₀t) + bₖ·sin(2πkf₀t)]
```

**Coefficients for this square wave:**
- a₀ = 0 (DC component is zero - signal averages to zero)
- aₖ = 0 for all k (cosine coefficients are zero)
- bₖ = 4/(πk) for odd k, bₖ = 0 for even k

**Simplified form:**
```
x(t) = Σ[4/(πk)·sin(2πkf₀t)]  for odd k only
```

## Files

### 1. `fourier_square_wave.m` (Main Function)
Computes Fourier series coefficients and reconstructs the signal.

**Usage:**
```matlab
[a_k, b_k, x_fourier] = fourier_square_wave(N, f0, t)
```

**Parameters:**
- `N`: Number of Fourier coefficients (try 10, 50, 100)
- `f0`: Fundamental frequency in Hz (try 0.5, 1, 2, 5)
- `t`: Time vector for signal reconstruction

**Returns:**
- `a_k`: Cosine coefficients [a₀, a₁, ..., aₙ]
- `b_k`: Sine coefficients [b₁, ..., bₙ]
- `x_fourier`: Reconstructed signal

### 2. `plot_fourier_series.m` (Solution to Problem 2b & 2c)
Main script that:
- Prompts user to enter parameters via input dialog
- Computes Fourier series coefficients
- Plots 2 side-by-side graphs:
  - Original square wave
  - Fourier series approximation
- Prints coefficients to console

**To run:** Execute `plot_fourier_series` in MATLAB

**Input parameters (via pop-up dialog):**
- `N`: Number of Fourier coefficients (try 1, 3, 10, 50, 100)
- `f0`: Fundamental frequency in Hz (try 0.5, 1, 2, 5)
- `Start time`: Beginning of time vector in seconds (default: -2)
- `End time`: End of time vector in seconds (default: 2)
- `Number of time points`: Resolution of time vector (default: 1000)

**Try different values to observe:**
- Low N (1-5): Poor approximation
- High N (50+): Excellent approximation with Gibbs phenomenon
- Different f₀: Changes oscillation frequency
- Larger time range: See more periods
- More time points: Smoother curves

## Key Observations

### Effect of N (Number of Harmonics)
- **N = 1**: Very rough approximation (just fundamental frequency)
- **N = 3**: Basic square wave shape appears
- **N = 10**: Good approximation with visible Gibbs phenomenon
- **N = 50+**: Excellent approximation, Gibbs overshoot ~9% remains

### Effect of f₀ (Fundamental Frequency)
- **Higher f₀**: Faster oscillation, shorter period (T₀ = 1/f₀)
- **Lower f₀**: Slower oscillation, longer period
- f₀ only affects time scale, not approximation quality

### Gibbs Phenomenon
At discontinuities, the Fourier series exhibits ~9% overshoot that persists regardless of N. This is a fundamental property of Fourier series approximations of discontinuous functions.

## Example Results

For N=10, f₀=1 Hz:
- b₁ = 4/π ≈ 1.273
- b₃ = 4/(3π) ≈ 0.424
- b₅ = 4/(5π) ≈ 0.255
- b₇ = 4/(7π) ≈ 0.182
- Even harmonics (b₂, b₄, b₆...) = 0

## How to Use

1. Open MATLAB
2. Navigate to the Problem2 folder
3. Run `plot_fourier_series`
4. Enter desired values for N and f₀ in the pop-up dialog
5. Observe the original and approximation plots
6. Experiment with different N and f₀ values to see effects
