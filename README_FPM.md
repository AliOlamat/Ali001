# FPM Visualization

Simple Fourier Ptychography Microscopy (FPM) reconstruction visualization.

## Overview

**Fourier Ptychography Microscopy (FPM)** is a computational imaging technique that captures high-resolution images by combining multiple low-resolution measurements taken at different illumination angles.

### Key Concepts

1. **Low-Resolution Measurements**: Images captured with different angle illuminations
2. **Fourier Space Stitching**: Each measurement covers a different region in Fourier space
3. **High-Resolution Reconstruction**: Combining all measurements reconstructs a high-resolution image

## What This Script Does

- **Simulates FPM**: Creates a synthetic sample object and simulates measurements from different illumination angles
- **Reconstructs**: Uses frequency stitching to combine measurements into a high-resolution image
- **Visualizes**: Displays the original object, all measurements, and the reconstruction

## Usage

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Run Visualization
```bash
python fpm_visualization.py
```

This will:
- Generate a sample object with circles and lines
- Simulate 9 measurements from different illumination angles (with noise)
- Reconstruct a high-resolution image
- Display results in a figure and save as `fpm_reconstruction.png`

## Parameters

In `main()`, you can adjust:

- **object_size**: Size of the test object (default: 64x64)
- **num_illumination_angles**: Number of different illumination angles (default: 9)
- **aperture_radius**: Detection aperture size in Fourier space (default: 0.3)
- **noise_level**: Gaussian noise added to measurements (default: 0.02)

## How FPM Works

```
Original Object
    ↓
Illuminate from different angles θ
    ↓
Capture low-res images (limited by aperture)
    ↓
FFT each measurement
    ↓
Place FFT in high-res Fourier space based on illumination angle
    ↓
Combine all frequency components
    ↓
IFFT to get high-resolution reconstruction
```

## Output

The visualization shows:
- **Top-left**: Original test object
- **Top-right**: Reconstructed high-resolution image
- **Bottom**: Low-resolution measurements from different illumination angles (colored heat map)

## References

Fourier Ptychography: A microscopy technique for high-resolution imaging using computational methods.

G. Zheng et al. "Wide-field, high-resolution, two-dimensional microscopy" *Nature Photonics* (2013)
