%% Simple FPM (Fourier Ptychography Microscopy) Reconstruction Visualization

clear all; close all; clc;

fprintf('Generating FPM visualization...\n');

% Parameters
object_size = 64;
num_illumination_angles = 9;
aperture_radius = 0.3;

% Generate sample object
fprintf('Creating sample object...\n');
obj = create_sample_object(object_size);

% Simulate FPM measurements
fprintf('Simulating FPM measurements (%d angles)...\n', num_illumination_angles);
[measurements, angles] = simulate_fpm_measurement(...
    obj, ...
    num_illumination_angles, ...
    aperture_radius, ...
    0.02);

% Reconstruct
fprintf('Reconstructing high-resolution image...\n');
reconstruction = simple_fpm_reconstruction(...
    measurements, ...
    num_illumination_angles, ...
    aperture_radius);

% Visualize
fprintf('Creating visualization...\n');
fig = visualize_fpm(obj, measurements, reconstruction, angles);

% Save figure
saveas(fig, 'fpm_reconstruction.png');
fprintf('Saved: fpm_reconstruction.png\n');

% Print statistics
fprintf('\n--- Reconstruction Statistics ---\n');
fprintf('Object size: %d x %d\n', size(obj, 1), size(obj, 2));
fprintf('Measurement size: %d x %d\n', size(measurements, 1), size(measurements, 2));
fprintf('Reconstruction size: %d x %d\n', size(reconstruction, 1), size(reconstruction, 2));
fprintf('Number of illumination angles: %d\n', num_illumination_angles);
fprintf('Aperture radius: %.1f\n', aperture_radius);

%% Function: Create sample object
function obj = create_sample_object(size_obj)
    % Create a simple sample object with features.

    x = linspace(-1, 1, size_obj);
    y = linspace(-1, 1, size_obj);
    [X, Y] = meshgrid(x, y);

    % Create a pattern with circles and lines
    circle1 = exp(-((X-0.3).^2 + (Y-0.3).^2) / 0.05);
    circle2 = exp(-((X+0.3).^2 + (Y+0.3).^2) / 0.05);
    line = exp(-(X.^2) / 0.01) .* (abs(Y) < 0.3);

    obj = abs(circle1 + circle2 + line);
end

%% Function: Simulate FPM measurements
function [measurements, angles] = simulate_fpm_measurement(...
    obj, num_angles, aperture_radius, noise_level)
    % Simulate FPM measurements from different illumination angles.
    %
    % Parameters:
    % - obj: Object to measure
    % - num_angles: Number of illumination angles
    % - aperture_radius: Radius of detection aperture (normalized)
    % - noise_level: Gaussian noise standard deviation
    %
    % Returns:
    % - measurements: Array of low-resolution images
    % - angles: Angles used

    size_obj = size(obj, 1);
    fft_size = size_obj * 2;  % Padding for convolution

    % Create aperture in Fourier space
    freq = fftshift(fftfreq(fft_size));
    freq = freq(1:fft_size);
    [freq_x, freq_y] = meshgrid(freq, freq);
    freq_2d = sqrt(freq_x.^2 + freq_y.^2);
    aperture = (freq_2d <= aperture_radius);

    measurements = [];
    angles = linspace(0, 2*pi, num_angles+1);
    angles = angles(1:num_angles);  % Remove last point

    % Pad object
    padded_obj = padarray(obj, [size_obj/2, size_obj/2], 'both');

    for m = 1:num_angles
        angle = angles(m);

        % Create illumination wave
        ill_kx = 0.4 * cos(angle);
        ill_ky = 0.4 * sin(angle);

        % Spatial domain coordinates
        y_coord = 1:fft_size;
        x_coord = 1:fft_size;
        [X, Y] = meshgrid(x_coord - fft_size/2, y_coord - fft_size/2);

        % Illumination pattern
        illumination = exp(1j * 2 * pi * (ill_kx * X + ill_ky * Y) / fft_size);

        % Multiply object by illumination
        illuminated = padded_obj .* illumination;

        % FFT
        illuminated_fft = fft2(illuminated);

        % Apply aperture
        aperture_2d = [aperture, aperture; aperture, aperture];
        aperture_centered = fftshift(aperture_2d);
        filtered = illuminated_fft .* aperture_centered;

        % IFFT back to spatial domain
        measurement = abs(ifft2(filtered)).^2;

        % Crop to original size and extract low-res
        crop_start = floor(size_obj/4) + 1;
        crop_end = crop_start + size_obj - 1;
        measurement_cropped = measurement(crop_start:crop_end, crop_start:crop_end);

        % Downsample
        downsampled = measurement_cropped(1:2:end, 1:2:end);

        % Add noise
        noise = noise_level * max(downsampled(:)) * randn(size(downsampled));
        noisy = downsampled + noise;
        noisy = max(noisy, 0);  % Remove negative values

        measurements = cat(3, measurements, noisy);
    end
end

%% Function: Simple FPM reconstruction
function reconstruction = simple_fpm_reconstruction(measurements, num_angles, aperture_radius)
    % Simple FPM reconstruction using frequency stitching.
    %
    % Parameters:
    % - measurements: Array of low-resolution images
    % - num_angles: Number of illumination angles
    % - aperture_radius: Aperture radius used
    %
    % Returns:
    % - reconstructed: High-resolution reconstruction

    size_meas = size(measurements, 1);
    hires_size = size_meas * 2;

    % Initialize high-resolution Fourier space
    hires_fft = zeros(hires_size, hires_size);
    weight = zeros(hires_size, hires_size);

    angles = linspace(0, 2*pi, num_angles+1);
    angles = angles(1:num_angles);  % Remove last point

    for m = 1:num_angles
        angle = angles(m);

        % FFT of measurement
        padded_meas = padarray(measurements(:,:,m), [size_meas/2, size_meas/2], 'both');
        meas_fft = fft2(padded_meas);

        % Illumination wavevector
        ill_kx = 0.4 * cos(angle);
        ill_ky = 0.4 * sin(angle);

        % Shift in Fourier space
        shift_x = round(ill_kx * hires_size);
        shift_y = round(ill_ky * hires_size);

        % Place measurement in high-res Fourier space
        if abs(shift_x) < hires_size/2 && abs(shift_y) < hires_size/2
            meas_fft_shifted = circshift(meas_fft, [shift_y, shift_x]);
            hires_fft = hires_fft + meas_fft_shifted;
            weight = weight + 1;
        end
    end

    % Normalize
    hires_fft = hires_fft ./ (weight + 1e-10);

    % IFFT to get reconstruction
    reconstruction_temp = abs(ifft2(hires_fft)).^2;
    reconstruction = sqrt(reconstruction_temp);
end

%% Function: Visualize FPM
function fig = visualize_fpm(obj, measurements, reconstruction, angles)
    % Visualize FPM process and results.

    num_measurements = size(measurements, 3);

    fig = figure('Position', [100, 100, 1400, 900]);

    % Original object
    subplot(3, 4, [1, 2]);
    imagesc(obj);
    colormap('gray');
    colorbar;
    title('Original Object', 'FontSize', 12, 'FontWeight', 'bold');
    axis off;

    % Reconstruction
    subplot(3, 4, [3, 4]);
    imagesc(reconstruction);
    colormap('gray');
    colorbar;
    title('FPM Reconstruction', 'FontSize', 12, 'FontWeight', 'bold');
    axis off;

    % Show sample measurements
    num_display = min(8, num_measurements);
    for i = 1:num_display
        row = floor((i-1)/4) + 2;
        col = mod(i-1, 4) + 1;
        idx = (row - 2) * 4 + col;

        subplot(3, 4, idx);
        imagesc(measurements(:,:,i));
        colormap('hot');
        colorbar;
        angle_deg = rad2deg(angles(i));
        title(sprintf('\\theta = %.1f°', angle_deg), 'FontSize', 10);
        axis off;
    end

    sgtitle('FPM (Fourier Ptychography Microscopy) Reconstruction', ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% Helper function: FFT frequency
function freq = fftfreq(N)
    % Returns frequency array for FFT (similar to numpy.fft.fftfreq)
    freq = [(0:N/2-1), (-N/2:-1)] / N;
end
