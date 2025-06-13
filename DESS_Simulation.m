clear; clc; close all;

%% Phantom Dimensions and Initialization
nx = 256; ny = 256;
phantom = zeros(nx, ny); 
[X, Y] = meshgrid(1:nx, 1:ny); 
S0 = 100;
% Initialize T1 and T2 Maps
T1_map = zeros(nx, ny);
T2_map = zeros(nx, ny);

%% Define Tissue Masks
% Bone (Femur and Tibia)
femur_mask = ((X - 128).^2 / 800 + (Y - 200).^2 / 4000) < 1;
tibia_mask = ((X - 128).^2 / 800 + (Y - 60).^2 / 4000) < 1;

% Cartilage (Femur and Tibia)
cartilage_mask_femur = ((X - 128).^2 / 1000 + (Y - 200).^2 / 4200) < 1 & ~femur_mask;
cartilage_mask_tibia = ((X - 128).^2 / 1000 + (Y - 60).^2 / 4200) < 1 & ~tibia_mask;

% Meniscus
meniscus_mask = (Y > 120 & Y < 180 & abs(X - 128) < 20);

%% Assign T1 and T2 Values
T1_map(femur_mask | tibia_mask) = 246;  % Bone
T2_map(femur_mask | tibia_mask) = 0.04;

T1_map(cartilage_mask_femur | cartilage_mask_tibia) = 900;  % Cartilage
T2_map(cartilage_mask_femur | cartilage_mask_tibia) = 39;

T1_map(meniscus_mask) = 960;  % Meniscus
T2_map(meniscus_mask) = 26.7;

%% Simulation Parameters
FA = deg2rad(30);  % Flip angle in radians
TR_DESS = 20;      % Repetition time for DESS in ms
TE_DESS = 5;       % Echo time for DESS in ms
M0 = 100;  % Proton Density (Assumed Constant)
sigma = 1; % Noise level (absolute noise standard deviation)

%% Initialize DESS Signal Matrices for Whole Phantom
signal_DESS_FID = zeros(nx, ny);
signal_DESS_SE = zeros(nx, ny);
noisy_signal_DESS_FID = zeros(nx, ny);
noisy_signal_DESS_SE = zeros(nx, ny);

%% DESS Signal Calculations with Updated Equations (Whole Phantom)
for x = 1:nx
    for y = 1:ny
        if T1_map(x, y) > 0 && T2_map(x, y) > 0
            % Compute relaxation factors
            E1 = exp(-TR_DESS / T1_map(x, y));
            E2 = exp(-TE_DESS / T2_map(x, y));
            
            % Compute r factor and S2/S1 ratio
            r = sqrt((1 - E2^2) / ((1 - E1 * cos(FA))^2 - E2^2 * (E1 - cos(FA))^2));
            S2_S1_ratio = exp(-2 * (TR_DESS - TE_DESS) / T2_map(x,y)) * sin(FA/2)^2 * ((1 + exp(-TR_DESS / T1_map(x,y))) / (1 - cos(FA) * exp(-TR_DESS / T1_map(x,y))));
            
            % Calculate FID and SE signals using updated equations
            signal_DESS_FID(x, y) = M0 * tan(FA/2) * (1 - (E1 - cos(FA)) * r) * exp(-TE_DESS / T2_map(x, y));
            signal_DESS_SE(x, y) = signal_DESS_FID(x, y) * S2_S1_ratio;
            
            % Add Gaussian noise
            noisy_signal_DESS_FID(x, y) = signal_DESS_FID(x, y) + sigma * randn;
            noisy_signal_DESS_SE(x, y) = signal_DESS_SE(x, y) + sigma * randn;
        end
    end
end

%% Visualize Noisy DESS FID Signal
figure;
imagesc(noisy_signal_DESS_FID);
colormap('hot');
colorbar;
title('Noisy DESS FID Signal with Flip Angle');

%% Visualize Noisy DESS SE Signal
figure;
imagesc(noisy_signal_DESS_SE);
colormap('hot');
colorbar;
title('Noisy DESS SE Signal with Flip Angle');

%% Simulate K-Space Acquisition and Image Reconstruction
% Fourier Transform (K-Space)
k_space_S1 = fftshift(fft2(noisy_signal_DESS_FID));  % Apply 2D FFT
k_space_S2 = fftshift(fft2(noisy_signal_DESS_SE));

% Visualize K-Space
figure;
imagesc(log(abs(k_space_S1)+1));
colormap('gray');
colorbar;
title('K-Space of S1 (FID)');

figure;
imagesc(log(abs(k_space_S2)+1));
colormap('gray');
colorbar;
title('K-Space of S2 (SE)');

% Reconstruct Images using Inverse FFT
reconstructed_S1 = abs(ifft2(ifftshift(k_space_S1)));
reconstructed_S2 = abs(ifft2(ifftshift(k_space_S2)));

% Visualize Reconstructed Images
figure;
imagesc(reconstructed_S1);
colormap('hot');
colorbar;
title('Reconstructed S1 (FID) Image');

figure;
imagesc(reconstructed_S2);
colormap('hot');
colorbar;
title('Reconstructed S2 (SE) Image');