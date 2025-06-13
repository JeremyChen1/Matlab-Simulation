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

%% Monte Carlo Simulation for Parameter Optimization (One Voxel per Tissue)
% Define T1 and T2 values for a single voxel in each tissue
T1_meniscus = 960;   % T1 for meniscus (ms)
T2_meniscus = 26.7;  % T2 for meniscus (ms)

T1_cartilage = 900;  % T1 for cartilage (ms)
T2_cartilage = 39;   % T2 for cartilage (ms)

% Proton density (assumed constant)
M0 = 100;

% Noise level (sigma already defined)
% Define the range of MRI parameters to optimize
TE_range = linspace(3, 10, 30);  % Example range for TE (ms)
TR_range = linspace(15, 30, 30);  % Example range for TR (ms)
FA_range = linspace(10, 50, 30);  % Example range for flip angle (degrees)

% Number of Monte Carlo iterations
num_iterations = 5000;

% Initialize variables to store results
SD_T2_meniscus = zeros(length(TE_range), length(TR_range), length(FA_range));
SD_T2_cartilage = zeros(length(TE_range), length(TR_range), length(FA_range));

% Loop over parameter ranges
for te_idx = 1:length(TE_range)
    for tr_idx = 1:length(TR_range)
        for fa_idx = 1:length(FA_range)
            TE_DESS = TE_range(te_idx);
            TR_DESS = TR_range(tr_idx);
            FA = deg2rad(FA_range(fa_idx));
            
            % Calculate the constant factor C for each tissue type
            C_meniscus = sin(FA/2)^2 * ((1+exp(-TR_DESS/T1_meniscus))/(1-cos(FA)*exp(-TR_DESS/T1_meniscus)));
            C_cartilage = sin(FA/2)^2 * ((1+exp(-TR_DESS/T1_cartilage))/(1-cos(FA)*exp(-TR_DESS/T1_cartilage)));
            
            % Initialize arrays to store T2 estimates for each Monte Carlo iteration
            T2_estimates_meniscus = zeros(num_iterations, 1);
            T2_estimates_cartilage = zeros(num_iterations, 1);
            
            % Monte Carlo simulation for each parameter set
            for iter = 1:num_iterations
                % Compute relaxation factors for meniscus voxel
                E1_meniscus = exp(-TR_DESS / T1_meniscus);
                E2_meniscus = exp(-TE_DESS / T2_meniscus);
                r_meniscus = sqrt((1 - E2_meniscus^2) / ((1 - E1_meniscus * cos(FA))^2 - E2_meniscus^2 * (E1_meniscus - cos(FA))^2));
                % True S2/S1 ratio (without noise)
                S2_S1_ratio_meniscus = exp(-2*(TR_DESS-TE_DESS)/T2_meniscus);
                
                % True signals (before noise)
                signal_FID_meniscus = M0 * tan(FA/2) * (1 - (E1_meniscus - cos(FA)) * r_meniscus) * exp(-TE_DESS/T2_meniscus);
                signal_SE_meniscus = signal_FID_meniscus * S2_S1_ratio_meniscus * C_meniscus;
                
                % Add Gaussian noise
                noisy_FID_meniscus = signal_FID_meniscus + sigma * randn;
                noisy_SE_meniscus = signal_SE_meniscus + sigma * randn;
                
                % Measured ratio from noisy signals
                R_meas_meniscus = noisy_SE_meniscus / noisy_FID_meniscus;
                % Invert to estimate T2 (only if ratio positive and valid)
                if R_meas_meniscus > 0 && (R_meas_meniscus/C_meniscus) > 0
                    T2_estimates_meniscus(iter) = -2*(TR_DESS-TE_DESS)/log(R_meas_meniscus/C_meniscus);
                else
                    T2_estimates_meniscus(iter) = NaN;
                end
                
                % Compute relaxation factors for cartilage voxel
                E1_cartilage = exp(-TR_DESS / T1_cartilage);
                E2_cartilage = exp(-TE_DESS / T2_cartilage);
                r_cartilage = sqrt((1 - E2_cartilage^2) / ((1 - E1_cartilage * cos(FA))^2 - E2_cartilage^2 * (E1_cartilage - cos(FA))^2));
                % True S2/S1 ratio (without noise)
                S2_S1_ratio_cartilage = exp(-2*(TR_DESS-TE_DESS)/T2_cartilage);
                
                % True signals (before noise)
                signal_FID_cartilage = M0 * tan(FA/2) * (1 - (E1_cartilage - cos(FA)) * r_cartilage) * exp(-TE_DESS/T2_cartilage);
                signal_SE_cartilage = signal_FID_cartilage * S2_S1_ratio_cartilage * C_cartilage;
                
                % Add Gaussian noise
                noisy_FID_cartilage = signal_FID_cartilage + sigma * randn;
                noisy_SE_cartilage = signal_SE_cartilage + sigma * randn;
                
                % Measured ratio from noisy signals
                R_meas_cartilage = noisy_SE_cartilage / noisy_FID_cartilage;
                % Invert to estimate T2
                if R_meas_cartilage > 0 && (R_meas_cartilage/C_cartilage) > 0
                    T2_estimates_cartilage(iter) = -2*(TR_DESS-TE_DESS)/log(R_meas_cartilage/C_cartilage);
                else
                    T2_estimates_cartilage(iter) = NaN;
                end
            end
            
            % Remove NaN entries from the estimates
            T2_estimates_meniscus = T2_estimates_meniscus(~isnan(T2_estimates_meniscus));
            T2_estimates_cartilage = T2_estimates_cartilage(~isnan(T2_estimates_cartilage));
            
            % Calculate the standard deviation of T2 estimates
            if ~isempty(T2_estimates_meniscus)
                SD_T2_meniscus(te_idx, tr_idx, fa_idx) = std(T2_estimates_meniscus);
            else
                SD_T2_meniscus(te_idx, tr_idx, fa_idx) = NaN;
            end
            
            if ~isempty(T2_estimates_cartilage)
                SD_T2_cartilage(te_idx, tr_idx, fa_idx) = std(T2_estimates_cartilage);
            else
                SD_T2_cartilage(te_idx, tr_idx, fa_idx) = NaN;
            end
        end
    end
end

% Find optimal parameters for meniscus
[min_SD_meniscus, min_idx_meniscus] = min(SD_T2_meniscus(:));
[te_opt_idx_meniscus, tr_opt_idx_meniscus, fa_opt_idx_meniscus] = ind2sub(size(SD_T2_meniscus), min_idx_meniscus);
TE_opt_meniscus = TE_range(te_opt_idx_meniscus);
TR_opt_meniscus = TR_range(tr_opt_idx_meniscus);
FA_opt_meniscus = FA_range(fa_opt_idx_meniscus);

% Find optimal parameters for cartilage
[min_SD_cartilage, min_idx_cartilage] = min(SD_T2_cartilage(:));
[te_opt_idx_cartilage, tr_opt_idx_cartilage, fa_opt_idx_cartilage] = ind2sub(size(SD_T2_cartilage), min_idx_cartilage);
TE_opt_cartilage = TE_range(te_opt_idx_cartilage);
TR_opt_cartilage = TR_range(tr_opt_idx_cartilage);
FA_opt_cartilage = FA_range(fa_opt_idx_cartilage);

fprintf('Optimal Parameters for Meniscus: TE = %.2f ms, TR = %.2f ms, Flip Angle = %.2f degrees\n', TE_opt_meniscus, TR_opt_meniscus, FA_opt_meniscus);
fprintf('Optimal Parameters for Cartilage: TE = %.2f ms, TR = %.2f ms, Flip Angle = %.2f degrees\n', TE_opt_cartilage, TR_opt_cartilage, FA_opt_cartilage);

% Remove NaN values if any
T2_estimates_meniscus = T2_estimates_meniscus(~isnan(T2_estimates_meniscus));
T2_estimates_cartilage = T2_estimates_cartilage(~isnan(T2_estimates_cartilage));

%% Check Distribution of T2 Estimates Before SD Calculation
figure;
histogram(T2_estimates_cartilage, 50);
title('Distribution of T2 Estimates (Cartilage)');
xlabel('T2 (ms)');
ylabel('Frequency');

figure;
histogram(T2_estimates_meniscus, 50);
title('Distribution of T2 Estimates (Meniscus)');
xlabel('T2 (ms)');
ylabel('Frequency');

% Calculate mean and standard deviation of T2 estimates
mean_T2_meniscus = mean(T2_estimates_meniscus);
std_T2_meniscus = std(T2_estimates_meniscus);

mean_T2_cartilage = mean(T2_estimates_cartilage);
std_T2_cartilage = std(T2_estimates_cartilage);

fprintf('Meniscus: Mean T2 = %.2f ms, SD = %.2f ms\n', mean_T2_meniscus, std_T2_meniscus);
fprintf('Cartilage: Mean T2 = %.2f ms, SD = %.2f ms\n', mean_T2_cartilage, std_T2_cartilage);

%% Plot
% Replace NaN values with zero
SD_T2_cartilage(isnan(SD_T2_cartilage)) = 0;
SD_T2_meniscus(isnan(SD_T2_meniscus)) = 0;

% Flip Angle Index
fa_index = round(length(FA_range)/2);

% Extract SD data for chosen FA
SD_T2_cartilage_fixed_FA = SD_T2_cartilage(:, :, fa_index);
SD_T2_meniscus_fixed_FA  = SD_T2_meniscus(:, :, fa_index);

% Ensure values are reasonable for visualization
SD_T2_cartilage_fixed_FA(SD_T2_cartilage_fixed_FA > 1e4) = 1e4; % Clip extreme values
SD_T2_meniscus_fixed_FA(SD_T2_meniscus_fixed_FA > 1e4) = 1e4;
% Compute the constraint region
[TE_grid, TR_grid] = meshgrid(TE_range, TR_range);
TE_max_grid = TR_grid / 2 - 4.31; % Maximum allowable TE from scanner hardware

% Create a mask where TE exceeds TE_max
invalid_region = TE_grid > TE_max_grid;

% Plot heatmap for Cartilage
figure;
hold on;
imagesc(TE_range, TR_range, SD_T2_cartilage_fixed_FA); % Plot SD heatmap
set(gca, 'YDir', 'normal'); % Flip Y-axis for correct TR ordering
colormap('parula');
colorbar;
clim([0,50]);
xlabel('TE (ms)');
ylabel('TR (ms)');
title(['SD of T2 for Cartilage at FA = ', num2str(FA_range(fa_index)), '°']);

% Overlay invalid region (masking out the non-allowed TE-TR values)
contour(TE_grid, TR_grid, invalid_region, [0.5, 0.5], 'r', 'LineWidth', 2);

% Add legend for the contour line only
legend('Forbidden Region');

hold off;

% Plot heatmap for Meniscus
figure;
hold on;
imagesc(TE_range, TR_range, SD_T2_meniscus_fixed_FA);
set(gca, 'YDir', 'normal'); 
colormap('parula');
colorbar;
clim([0,50]);
xlabel('TE (ms)');
ylabel('TR (ms)');
title(['SD of T2 for Meniscus at FA = ', num2str(FA_range(fa_index)), '°']);

% Overlay invalid region (forbidden TE-TR)
contour(TE_grid, TR_grid, invalid_region, [0.5, 0.5], 'r', 'LineWidth', 2);

% Add legend for the contour line only
legend('Forbidden Region');

hold off;

%% 3D Surface Plot
% Cartilage
[TE_grid, TR_grid] = meshgrid(TE_range, TR_range);
figure;
surf(TE_grid, TR_grid, SD_T2_cartilage_fixed_FA); % Log scale
xlabel('TE (ms)');
ylabel('TR (ms)');
zlabel('log10(SD of T2) (Cartilage)');
title(['SD of T2 for Cartilage at FA = ', num2str(FA_range(fa_index)), '°']);
colormap('parula');  
colorbar;
shading interp;
clim([0,50]);
% Meniscus
figure;
surf(TE_grid, TR_grid, SD_T2_meniscus_fixed_FA); % Log scale0
xlabel('TE (ms)');
ylabel('TR (ms)');
zlabel('log10(SD of T2) (Meniscus)');
title(['SD of T2 for Meniscus at FA = ', num2str(FA_range(fa_index)), '°']);
colormap('parula');
colorbar;
shading interp;
clim([0,50]);
% 4.31