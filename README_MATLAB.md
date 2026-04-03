# FPM Visualization - MATLAB Version

MATLAB implementation of Fourier Ptychography Microscopy (FPM) reconstruction visualization.

## Overview

This is a direct port of the Python FPM visualization script to MATLAB, maintaining the same functionality and algorithm.

## Files

- **fpm_visualization.m** - Main MATLAB script containing all functions

## Requirements

- MATLAB R2016b or later
- Image Processing Toolbox (optional, for enhanced visualization)
- Signal Processing Toolbox (for FFT operations)

## Usage

1. Open MATLAB
2. Navigate to the directory containing `fpm_visualization.m`
3. Run the script:

```matlab
fpm_visualization
```

Or execute from command line:
```bash
matlab -r "fpm_visualization"
```

## What the Script Does

- **create_sample_object()**: Generates a test object with circles and lines
- **simulate_fpm_measurement()**: Simulates low-resolution measurements from different illumination angles
- **simple_fpm_reconstruction()**: Reconstructs high-resolution image using Fourier space stitching
- **visualize_fpm()**: Creates comprehensive visualization

## Parameters

You can modify these parameters at the top of the script:

```matlab
object_size = 64;              % Size of the test object
num_illumination_angles = 9;   % Number of different illumination angles
aperture_radius = 0.3;         % Detection aperture size in Fourier space
```

## Output

The script generates:
- A figure window showing all results
- **fpm_reconstruction.png** - Saved visualization

The visualization includes:
- **Top-left**: Original test object
- **Top-right**: Reconstructed high-resolution image
- **Bottom 8 panels**: Low-resolution measurements from different illumination angles

## Algorithm

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

## Key Differences from Python Version

- Uses `padarray()` instead of numpy's `pad()`
- Uses `circshift()` instead of numpy's `roll()`
- Uses `cat(3, ...)` to build 3D measurement array
- MATLAB's FFT functions work directly on complex arrays
- Plotting uses MATLAB's `subplot()` and `sgtitle()`

## Performance Notes

- For larger object sizes or more angles, computation may take longer
- Memory usage scales with object_size² and num_illumination_angles
- Consider reducing `object_size` or `num_illumination_angles` for faster computation

## Troubleshooting

- **Out of memory**: Reduce `object_size` from 64 to 32
- **Slow performance**: Use fewer angles or smaller object size
- **Visualization issues**: Ensure Image Processing Toolbox is installed for `imagesc()`

## References

Fourier Ptychography: A microscopy technique for high-resolution imaging using computational methods.

G. Zheng et al. "Wide-field, high-resolution, two-dimensional microscopy" *Nature Photonics* (2013)
