# Problem 1: Signal Decomposition into Even and Odd Parts

## Overview
This folder contains MATLAB scripts to decompose signals into even and odd components.

## Theory
Any signal x(t) can be decomposed into:
- **Even part**: x_e(t) = 0.5 × [x(t) + x(-t)]
- **Odd part**: x_o(t) = 0.5 × [x(t) - x(-t)]
- **Verification**: x(t) = x_e(t) + x_o(t)

## Files

### 1. `decompose_signal.m` (Solution to Problem 1b)
Function that decomposes a signal vector into even and odd components.

**Usage:**
```matlab
[x_e, x_o] = decompose_signal(x)
```

### 2. `test_decomposition.m` (Solution to Problem 1c & 1d)
Test script that solves Problem 1 using the signal from part (a).

**Signal definition:**
- x(t) = 0 for t < 0
- x(t) = 0.5 for t = 0
- x(t) = 1 for t > 0

**Features:**
- Decomposes the signal
- Creates 4 plots with proper axes
- Shows discontinuities with open/filled circles
- Displays numerical results
- Verifies x = x_e + x_o

**To run:** Execute `test_decomposition` in MATLAB

### 3. `problem_1a_solution.txt`
Step-by-step manual solution for Problem 1(a) showing graphical determination of even and odd parts.

## Expected Results

For the given signal:
- **Even part**: x_e(t) = 0.5 (constant for all t)
- **Odd part**: x_o(t) = -0.5 for t<0, x_o(0)=0, x_o(t)=0.5 for t>0
- **Sum**: x_e(t) + x_o(t) = original signal ✓
