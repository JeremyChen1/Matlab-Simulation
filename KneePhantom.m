clear; clc; close all
% Knee Phantom with T1 and T2 Maps
nx = 256; ny = 256;  % Number of pixels in x&y direction
phantom = zeros(nx, ny);  % Start a blank phantom
[X, Y] = meshgrid(1:nx, 1:ny);  % Create a coordinate grid

% Initialize T1 and T2 maps
T1_map = zeros(nx, ny);
T2_map = zeros(nx, ny);

% Bone (Femur and Tibia)
femur_mask = ((X - 128).^2 / 800 + (Y - 200).^2 / 4000) < 1;
tibia_mask = ((X - 128).^2 / 800 + (Y - 60).^2 / 4000) < 1;

% Assign T1 and T2 values for bone
T1_map(femur_mask) = 246;
T2_map(femur_mask) = 0.04; 
T1_map(tibia_mask) = 246;
T2_map(tibia_mask) = 0.04;

% Cartilage (Femur and Tibia)
cartilage_mask_femur = ((X - 128).^2 / 1000 + (Y - 200).^2 / 4200) < 1 & ~femur_mask;  % Femoral cartilage
cartilage_mask_tibia = ((X - 128).^2 / 1000 + (Y - 60).^2 / 4200) < 1 & ~tibia_mask;  % Tibial cartilage

% Assign T1 and T2 values for cartilage
T1_map(cartilage_mask_femur) = 900;
T2_map(cartilage_mask_femur) = 39;
T1_map(cartilage_mask_tibia) = 900;
T2_map(cartilage_mask_tibia) = 39;

% Meniscus
meniscus_mask = (Y > 120 & Y < 180 & abs(X - 128) < 20);

% Assign T1 and T2 values for meniscus
T1_map(meniscus_mask) = 960;
T2_map(meniscus_mask) = 26.7; 

%% Simulate Signal Decay (T2 Relaxation Equation)
TE = linspace(10, 200, 11);  % Echo times (10 ms to 200 ms, 10 steps)
S0 = 100;  % Initial signal intensity
sigma = 5; % Standard deviation for Gaussian noise

% Initialize Signal Matrix
signal_T2 = zeros(nx, ny, length(TE));  % 3D matrix: (x, y, TE)
noisy_signal_T2 = zeros(nx, ny, length(TE));

% Calculation
for i = 1:length(TE)
    signal_T2(:, :, i) = S0 * exp(-TE(i) ./ T2_map);
    signal_T2(isnan(signal_T2)) = 0;  % Handle division by zero for bone T2=0
    % Add Gaussian noise
    noisy_signal_T2(:, :, i) = signal_T2(:, :, i) + sigma * randn(size(signal_T2(:, :, i)));
end
% % Visualisation
% for i = 1:length(TE)
%     figure;
%     imagesc(noisy_signal_T2(:, :, i));
%     colormap('hot');
%     colorbar;
%     title(['Noisy Signal Decay at TE = ', num2str(TE(i)), ' ms']);
%     clim([0 S0]);  % Adjust color axis to the signal range
% end

%% Simulate Signal Recovery (T1 Relaxation Equation)
TR = linspace(10, 1000, 11);  % Repetition times (10 ms to 1000 ms, 10 steps)
M0 = 100;  % Initial magnetization

% Initialize Signal Matrix
signal_T1 = zeros(nx, ny, length(TR));  % 3D matrix: (x, y, TR)
noisy_signal_T1 = zeros(nx, ny, length(TR));

% Calculation
for i = 1:length(TR)
    signal_T1(:, :, i) = M0 * (1 - exp(-TR(i) ./ T1_map));
    signal_T1(isnan(signal_T1)) = 0;  % Handle division by zero for unassigned regions
    % Add Gaussian noise
    noisy_signal_T1(:, :, i) = signal_T1(:, :, i) + sigma * randn(size(signal_T1(:, :, i)));
end

% Visualisation
% for i = 1:length(TR)
%     figure;
%     imagesc(noisy_signal_T1(:, :, i));
%     colormap('hot');
%     colorbar;
%     title(['Noisy Signal Recovery at TR = ', num2str(TR(i)), ' ms']);
%     clim([0 M0]);  % Adjust color axis to the magnetization range
% end

%% Simulate DESS Sequence
FA = deg2rad(30);  % Flip angle in radians
TR_DESS = 20;      % Repetition time for DESS in ms
TE_DESS = 5;       % Echo time for DESS in ms

% Initialize DESS Signal Matrices
signal_DESS_FID = zeros(nx, ny);
signal_DESS_SE = zeros(nx, ny);
noisy_signal_DESS_FID = zeros(nx, ny);
noisy_signal_DESS_SE = zeros(nx, ny);

% DESS Signal Calculations
for x = 1:nx
    for y = 1:ny
        if T1_map(x, y) > 0 && T2_map(x, y) > 0
            % Calculate FID and SE signals
            signal_DESS_FID(x, y) = S0 * exp(-TE_DESS / T2_map(x, y));
            signal_DESS_SE(x, y) = S0 * exp(-2 * (TR_DESS - TE_DESS) / T2_map(x, y));
            % Add Gaussian noise
            noisy_signal_DESS_FID(x, y) = signal_DESS_FID(x, y) + sigma * randn;
            noisy_signal_DESS_SE(x, y) = signal_DESS_SE(x, y) + sigma * randn;
        end
    end
end

% Visualize Noisy DESS FID Signal
figure;
imagesc(noisy_signal_DESS_FID);
colormap('hot');
colorbar;
title('Noisy DESS FID Signal');

% Visualize Noisy DESS SE Signal
figure;
imagesc(noisy_signal_DESS_SE);
colormap('hot');
colorbar;
title('Noisy DESS SE Signal');


%% Plot DESS Signal for Cartilage and Meniscus with Theoretical Decay
% Calculate Mean Experimental Data for Cartilage and Meniscus
cartilage_DESS_FID = mean(noisy_signal_DESS_FID(cartilage_mask_femur));
cartilage_DESS_SE = mean(noisy_signal_DESS_SE(cartilage_mask_femur));

meniscus_DESS_FID = mean(noisy_signal_DESS_FID(meniscus_mask));
meniscus_DESS_SE = mean(noisy_signal_DESS_SE(meniscus_mask));

% Theoretical Curves for Cartilage
cartilage_T2_mean = mean(T2_map(cartilage_mask_femur));
TE_points = linspace(TE_DESS, TR_DESS, 100);
cartilage_FID_theoretical = S0 * exp(-TE_points / cartilage_T2_mean);
cartilage_SE_theoretical = S0 * exp(-(2 * TR_DESS - TE_points) / cartilage_T2_mean);

% Theoretical Curves for Meniscus
meniscus_T2_mean = mean(T2_map(meniscus_mask));
meniscus_FID_theoretical = S0 * exp(-TE_points / meniscus_T2_mean);
meniscus_SE_theoretical = S0 * exp(-(2 * TR_DESS - TE_points) / meniscus_T2_mean);

% Plot Cartilage Signals (FID and SE)
figure;
hold on;
plot(TE_points, cartilage_FID_theoretical, '--b', 'DisplayName', 'Cartilage FID (Theoretical)');
plot(TE_points, cartilage_SE_theoretical, '--r', 'DisplayName', 'Cartilage SE (Theoretical)');
plot([TE_DESS, TR_DESS], [cartilage_DESS_FID, cartilage_DESS_SE], '-ob', 'DisplayName', 'Cartilage FID (Experimental)');
plot([TE_DESS, TR_DESS], [cartilage_DESS_FID, cartilage_DESS_SE], '-or', 'DisplayName', 'Cartilage SE (Experimental)');
xlabel('Time (ms)');
ylabel('Signal Intensity');
title('Cartilage FID and SE Signals');
legend;
grid on;
hold off;

% Plot Meniscus Signals (FID and SE)
figure;
hold on;
plot(TE_points, meniscus_FID_theoretical, '--g', 'DisplayName', 'Meniscus FID (Theoretical)');
plot(TE_points, meniscus_SE_theoretical, '--k', 'DisplayName', 'Meniscus SE (Theoretical)');
plot([TE_DESS, TR_DESS], [meniscus_DESS_FID, meniscus_DESS_SE], '-xg', 'DisplayName', 'Meniscus FID (Experimental)');
plot([TE_DESS, TR_DESS], [meniscus_DESS_FID, meniscus_DESS_SE], '-xk', 'DisplayName', 'Meniscus SE (Experimental)');
xlabel('Time (ms)');
ylabel('Signal Intensity');
title('Meniscus FID and SE Signals');
legend;
grid on;
hold off;

%% Simulate K-Space Acquisition for DESS
% Fourier transform to simulate k-space acquisition
k_space_DESS_FID = fftshift(fft2(noisy_signal_DESS_FID));
k_space_DESS_SE = fftshift(fft2(noisy_signal_DESS_SE));

% Visualize K-Space
figure;
imagesc(log(abs(k_space_DESS_FID))); % Log scale for better visibility
colormap('gray');
colorbar;
title('K-Space of DESS FID Signal');

figure;
imagesc(log(abs(k_space_DESS_SE))); % Log scale for better visibility
colormap('gray');
colorbar;
title('K-Space of DESS SE Signal');

% Reconstruct the image from k-space (inverse Fourier transform)
reconstructed_DESS_FID = abs(ifft2(ifftshift(k_space_DESS_FID)));
reconstructed_DESS_SE = abs(ifft2(ifftshift(k_space_DESS_SE)));

% Visualize Reconstructed Images
figure;
imagesc(reconstructed_DESS_FID);
colormap('hot');
colorbar;
title('Reconstructed DESS FID Image');

figure;
imagesc(reconstructed_DESS_SE);
colormap('hot');
colorbar;
title('Reconstructed DESS SE Image');

