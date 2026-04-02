"""
Simple FPM (Fourier Ptychography Microscopy) Reconstruction Visualization
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import ndimage
from scipy.fft import fft2, ifft2, fftshift


def create_sample_object(size=128):
    """Create a simple sample object with features."""
    x = np.linspace(-1, 1, size)
    y = np.linspace(-1, 1, size)
    X, Y = np.meshgrid(x, y)

    # Create a pattern with circles and lines
    circle1 = np.exp(-((X-0.3)**2 + (Y-0.3)**2) / 0.05)
    circle2 = np.exp(-((X+0.3)**2 + (Y+0.3)**2) / 0.05)
    line = np.exp(-(X**2) / 0.01) * (np.abs(Y) < 0.3)

    obj = circle1 + circle2 + line
    return np.abs(obj)


def simulate_fpm_measurement(obj, num_angles=9, aperture_radius=0.3, noise_level=0.05):
    """
    Simulate FPM measurements from different illumination angles.

    Parameters:
    - obj: Object to measure
    - num_angles: Number of illumination angles
    - aperture_radius: Radius of detection aperture (normalized)
    - noise_level: Gaussian noise standard deviation

    Returns:
    - measurements: List of low-resolution images
    - illumination_angles: Angles used
    """
    size = obj.shape[0]
    fft_size = size * 2  # Padding for convolution

    # Create aperture in Fourier space
    freq = np.fft.fftfreq(fft_size)[:fft_size//2]
    freq_2d = np.sqrt(np.add.outer(freq**2, freq**2))
    aperture = (freq_2d <= aperture_radius).astype(float)

    measurements = []
    angles = np.linspace(0, 2*np.pi, num_angles, endpoint=False)

    # Pad object
    padded_obj = np.pad(obj, ((size//2, size//2), (size//2, size//2)), mode='constant')

    # FFT of object
    obj_fft = fft2(padded_obj)

    for angle in angles:
        # Create illumination wave
        ill_kx = 0.4 * np.cos(angle)
        ill_ky = 0.4 * np.sin(angle)

        # Spatial domain coordinates
        y_coord = np.arange(fft_size)
        x_coord = np.arange(fft_size)
        X, Y = np.meshgrid(x_coord - fft_size//2, y_coord - fft_size//2)

        # Illumination pattern
        illumination = np.exp(1j * 2 * np.pi * (ill_kx * X + ill_ky * Y) / fft_size)

        # Multiply object by illumination
        illuminated = padded_obj * illumination

        # FFT
        illuminated_fft = fft2(illuminated)

        # Apply aperture
        aperture_2d = np.tile(aperture, (2, 2))
        aperture_centered = fftshift(aperture_2d)
        filtered = illuminated_fft * aperture_centered

        # IFFT back to spatial domain
        measurement = np.abs(ifft2(filtered))**2

        # Crop to original size and extract low-res
        measurement_cropped = measurement[size//4:size//4+size, size//4:size//4+size]

        # Downsample
        downsampled = measurement_cropped[::2, ::2]

        # Add noise
        noisy = downsampled + np.random.normal(0, noise_level * downsampled.max(), downsampled.shape)
        noisy = np.maximum(noisy, 0)  # Remove negative values

        measurements.append(noisy)

    return np.array(measurements), angles


def simple_fpm_reconstruction(measurements, num_angles, aperture_radius=0.3):
    """
    Simple FPM reconstruction using frequency stitching.

    Parameters:
    - measurements: Array of low-resolution images
    - num_angles: Number of illumination angles
    - aperture_radius: Aperture radius used

    Returns:
    - reconstructed: High-resolution reconstruction
    """
    size = measurements[0].shape[0]
    hires_size = size * 2

    # Initialize high-resolution Fourier space
    hires_fft = np.zeros((hires_size, hires_size), dtype=complex)
    weight = np.zeros((hires_size, hires_size))

    angles = np.linspace(0, 2*np.pi, num_angles, endpoint=False)

    for m, angle in enumerate(angles):
        # FFT of measurement
        meas_fft = fft2(np.pad(measurements[m], ((size//2, size//2), (size//2, size//2)), mode='constant'))

        # Illumination wavevector
        ill_kx = 0.4 * np.cos(angle)
        ill_ky = 0.4 * np.sin(angle)

        # Shift in Fourier space
        shift_x = int(ill_kx * hires_size)
        shift_y = int(ill_ky * hires_size)

        # Place measurement in high-res Fourier space
        if abs(shift_x) < hires_size//2 and abs(shift_y) < hires_size//2:
            hires_fft = np.roll(np.roll(meas_fft, shift_y, axis=0), shift_x, axis=1)
            hires_fft += hires_fft
            weight += 1

    # Normalize
    hires_fft = hires_fft / (weight + 1e-10)

    # IFFT to get reconstruction
    reconstruction = np.abs(ifft2(hires_fft))**2

    return np.sqrt(reconstruction)


def visualize_fpm(obj, measurements, reconstruction, angles):
    """Visualize FPM process and results."""
    num_measurements = len(measurements)

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(3, 4, figure=fig, hspace=0.3, wspace=0.3)

    # Original object
    ax = fig.add_subplot(gs[0, :2])
    im = ax.imshow(obj, cmap='gray')
    ax.set_title('Original Object', fontsize=12, fontweight='bold')
    ax.axis('off')
    plt.colorbar(im, ax=ax)

    # Reconstruction
    ax = fig.add_subplot(gs[0, 2:])
    im = ax.imshow(reconstruction, cmap='gray')
    ax.set_title('FPM Reconstruction', fontsize=12, fontweight='bold')
    ax.axis('off')
    plt.colorbar(im, ax=ax)

    # Show sample measurements
    num_display = min(8, num_measurements)
    for i in range(num_display):
        row = 1 + i // 4
        col = i % 4
        ax = fig.add_subplot(gs[row, col])
        im = ax.imshow(measurements[i], cmap='hot')
        angle_deg = np.degrees(angles[i])
        ax.set_title(f'θ = {angle_deg:.1f}°', fontsize=10)
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046)

    plt.suptitle('FPM (Fourier Ptychography Microscopy) Reconstruction',
                 fontsize=14, fontweight='bold', y=0.995)

    return fig


def main():
    """Main visualization function."""
    print("Generating FPM visualization...")

    # Parameters
    object_size = 64
    num_illumination_angles = 9
    aperture_radius = 0.3

    # Generate sample object
    print("Creating sample object...")
    obj = create_sample_object(object_size)

    # Simulate FPM measurements
    print(f"Simulating FPM measurements ({num_illumination_angles} angles)...")
    measurements, angles = simulate_fpm_measurement(
        obj,
        num_angles=num_illumination_angles,
        aperture_radius=aperture_radius,
        noise_level=0.02
    )

    # Reconstruct
    print("Reconstructing high-resolution image...")
    reconstruction = simple_fpm_reconstruction(
        measurements,
        num_illumination_angles,
        aperture_radius=aperture_radius
    )

    # Visualize
    print("Creating visualization...")
    fig = visualize_fpm(obj, measurements, reconstruction, angles)

    # Save and show
    plt.savefig('fpm_reconstruction.png', dpi=150, bbox_inches='tight')
    print("Saved: fpm_reconstruction.png")
    plt.show()

    # Print statistics
    print("\n--- Reconstruction Statistics ---")
    print(f"Object size: {obj.shape}")
    print(f"Measurement size: {measurements[0].shape}")
    print(f"Reconstruction size: {reconstruction.shape}")
    print(f"Number of illumination angles: {num_illumination_angles}")
    print(f"Aperture radius: {aperture_radius}")


if __name__ == "__main__":
    main()
