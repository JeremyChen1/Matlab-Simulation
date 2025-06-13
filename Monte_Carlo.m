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

%% Monte Carlo Simulation for Parameter Optimization (One Voxel per Tissue)
% Define T1 and T2 values for a single voxel in each tissue
T1_meniscus = 960;   % T1 for meniscus (ms)
T2_meniscus = 26.7;  % T2 for meniscus (ms)

T1_cartilage = 900;  % T1 for cartilage (ms)
T2_cartilage = 39;   % T2 for cartilage (ms)

% Define the range of MRI parameters to optimize
TE_range = linspace(3, 9, 30);
TR_range = linspace(20, 40, 30);
FA_range = linspace(10, 50, 10);

% Number of Monte Carlo iterations
num_iterations = 10000;

% Initialize variables to store results
SD_T2_meniscus = NaN(length(TE_range), length(TR_range), length(FA_range));
SD_T2_cartilage = NaN(length(TE_range), length(TR_range), length(FA_range));

% Create a mask where TE exceeds TE_max (invalid region)
[TE_grid, TR_grid] = meshgrid(TE_range, TR_range);
invalid_region_mask = TE_grid > ((TR_grid / 2) - 4.31);

% Loop over parameter ranges
for te_idx = 1:length(TE_range)
    for tr_idx = 1:length(TR_range)
        for fa_idx = 1:length(FA_range)
            % Skip invalid regions
            if invalid_region_mask(tr_idx, te_idx)
                continue;
            end
            
            TE_DESS = TE_range(te_idx);
            TR_DESS = TR_range(tr_idx);
            FA = deg2rad(FA_range(fa_idx)); % Convert FA to radians

            % Monte Carlo simulation
            T2_estimates_meniscus = zeros(num_iterations, 1);
            T2_estimates_cartilage = zeros(num_iterations, 1);
            
            for iter = 1:num_iterations
                % Meniscus Calculations
                E1_meniscus = exp(-TR_DESS / T1_meniscus);
                E2_meniscus = exp(-TE_DESS / T2_meniscus);
                r_meniscus = sqrt((1 - E2_meniscus^2) / ((1 - E1_meniscus * cos(FA))^2 - E2_meniscus^2 * (E1_meniscus - cos(FA))^2));
                S2_S1_ratio_meniscus = exp(-2*(TR_DESS-TE_DESS)/T2_meniscus);
                signal_FID_meniscus = M0 * tan(FA/2) * (1 - (E1_meniscus - cos(FA)) * r_meniscus) * exp(-TE_DESS/T2_meniscus);
                signal_SE_meniscus = signal_FID_meniscus * S2_S1_ratio_meniscus;
                
                noisy_FID_meniscus = signal_FID_meniscus + sigma * randn;
                noisy_SE_meniscus = signal_SE_meniscus + sigma * randn;
                
                R_meas_meniscus = noisy_SE_meniscus / noisy_FID_meniscus;
                if R_meas_meniscus > 0
                    T2_estimates_meniscus(iter) = -2*(TR_DESS-TE_DESS)/log(R_meas_meniscus);
                else
                    T2_estimates_meniscus(iter) = NaN;
                end
                
                % Cartilage Calculations
                E1_cartilage = exp(-TR_DESS / T1_cartilage);
                E2_cartilage = exp(-TE_DESS / T2_cartilage);
                r_cartilage = sqrt((1 - E2_cartilage^2) / ((1 - E1_cartilage * cos(FA))^2 - E2_cartilage^2 * (E1_cartilage - cos(FA))^2));
                S2_S1_ratio_cartilage = exp(-2*(TR_DESS-TE_DESS)/T2_cartilage);
                signal_FID_cartilage = M0 * tan(FA/2) * (1 - (E1_cartilage - cos(FA)) * r_cartilage) * exp(-TE_DESS/T2_cartilage);
                signal_SE_cartilage = signal_FID_cartilage * S2_S1_ratio_cartilage;
                
                noisy_FID_cartilage = signal_FID_cartilage + sigma * randn;
                noisy_SE_cartilage = signal_SE_cartilage + sigma * randn;
                
                R_meas_cartilage = noisy_SE_cartilage / noisy_FID_cartilage;
                if R_meas_cartilage > 0
                    T2_estimates_cartilage(iter) = -2*(TR_DESS-TE_DESS)/log(R_meas_cartilage);
                else
                    T2_estimates_cartilage(iter) = NaN;
                end
            end
            
            % Compute SD(T2)
            SD_T2_meniscus(te_idx, tr_idx, fa_idx) = nanstd(T2_estimates_meniscus);
            SD_T2_cartilage(te_idx, tr_idx, fa_idx) = nanstd(T2_estimates_cartilage);
        end
    end
end

%% Optimizing TE, TR, and FA based on SD(T2)

% Mask out invalid regions in SD data
SD_T2_meniscus_valid = SD_T2_meniscus;
SD_T2_meniscus_valid(invalid_region_mask) = NaN;
SD_T2_cartilage_valid = SD_T2_cartilage;
SD_T2_cartilage_valid(invalid_region_mask) = NaN;

% Remove NaN values and compute the 10th percentile
valid_data_meniscus = SD_T2_meniscus_valid(~isnan(SD_T2_meniscus_valid));
valid_data_cartilage = SD_T2_cartilage_valid(~isnan(SD_T2_cartilage_valid));

threshold_meniscus = prctile(valid_data_meniscus, 5);
threshold_cartilage = prctile(valid_data_cartilage, 5);

% Find the best valid region with max TE
valid_meniscus_indices = find(SD_T2_meniscus_valid <= threshold_meniscus);
valid_cartilage_indices = find(SD_T2_cartilage_valid <= threshold_cartilage);

[te_valid_meniscus_idx, tr_valid_meniscus_idx, fa_valid_meniscus_idx] = ind2sub(size(SD_T2_meniscus_valid), valid_meniscus_indices);
[te_valid_cartilage_idx, tr_valid_cartilage_idx, fa_valid_cartilage_idx] = ind2sub(size(SD_T2_cartilage_valid), valid_cartilage_indices);

% Find the maximum TE within the valid indices
[max_TE_meniscus, idx_max_TE_meniscus] = max(TE_range(te_valid_meniscus_idx));
[max_TE_cartilage, idx_max_TE_cartilage] = max(TE_range(te_valid_cartilage_idx));

% Get the indices of all TE values that match max_TE_meniscus
max_TE_indices_meniscus = find(TE_range(te_valid_meniscus_idx) == max_TE_meniscus);
max_TE_indices_cartilage = find(TE_range(te_valid_cartilage_idx) == max_TE_cartilage);

% Within these, find the minimum TR
[min_TR_meniscus, idx_min_TR_meniscus] = min(TR_range(tr_valid_meniscus_idx(max_TE_indices_meniscus)));
[min_TR_cartilage, idx_min_TR_cartilage] = min(TR_range(tr_valid_cartilage_idx(max_TE_indices_cartilage)));

% Get the corresponding optimal TE, TR indices
te_opt_idx_meniscus = te_valid_meniscus_idx(max_TE_indices_meniscus(idx_min_TR_meniscus));
tr_opt_idx_meniscus = tr_valid_meniscus_idx(max_TE_indices_meniscus(idx_min_TR_meniscus));

te_opt_idx_cartilage = te_valid_cartilage_idx(max_TE_indices_cartilage(idx_min_TR_cartilage));
tr_opt_idx_cartilage = tr_valid_cartilage_idx(max_TE_indices_cartilage(idx_min_TR_cartilage));

% Extract SD values at optimal TE and TR for all FA values
SD_meniscus_fa_values = squeeze(SD_T2_meniscus_valid(te_opt_idx_meniscus, tr_opt_idx_meniscus, :));
SD_cartilage_fa_values = squeeze(SD_T2_cartilage_valid(te_opt_idx_cartilage, tr_opt_idx_cartilage, :));

% Find the minimum SD FA index
[min_SD_meniscus, fa_opt_idx_meniscus] = nanmin(SD_meniscus_fa_values);
[min_SD_cartilage, fa_opt_idx_cartilage] = nanmin(SD_cartilage_fa_values);

% Extract the actual values
TE_opt_meniscus = TE_range(te_opt_idx_meniscus);
TR_opt_meniscus = TR_range(tr_opt_idx_meniscus);
FA_opt_meniscus = FA_range(fa_opt_idx_meniscus);

TE_opt_cartilage = TE_range(te_opt_idx_cartilage);
TR_opt_cartilage = TR_range(tr_opt_idx_cartilage);
FA_opt_cartilage = FA_range(fa_opt_idx_cartilage);

% Print results
fprintf('Optimized Parameters for Meniscus: TE = %.2f ms, TR = %.2f ms, FA = %.2f°\n', ...
    TE_opt_meniscus, TR_opt_meniscus, FA_opt_meniscus);

fprintf('Optimized Parameters for Cartilage: TE = %.2f ms, TR = %.2f ms, FA = %.2f°\n', ...
    TE_opt_cartilage, TR_opt_cartilage, FA_opt_cartilage);


%% Compute Combined SD using Weighting Factors
k_m = 0.5;  % Weight for meniscus
k_c = 0.5;  % Weight for cartilage

% Compute combined SD
SD_T2_combined = sqrt((SD_T2_meniscus.^2 * k_m) + (SD_T2_cartilage.^2 * k_c));

% Mask out invalid regions
SD_T2_combined_valid = SD_T2_combined;
SD_T2_combined_valid(invalid_region_mask) = NaN;

% Extract valid values
valid_data_combined = SD_T2_combined_valid(~isnan(SD_T2_combined_valid));

% Compute the 10th percentile threshold
threshold_combined = prctile(valid_data_combined, 10);

% Find valid indices where SD_combined is within the lowest 10%
valid_combined_indices = find(SD_T2_combined_valid <= threshold_combined);

% Convert to TE, TR, FA indices
[te_valid_combined_idx, tr_valid_combined_idx, fa_valid_combined_idx] = ...
    ind2sub(size(SD_T2_combined_valid), valid_combined_indices);

% Find the maximum TE within this region
[max_TE_combined, idx_max_TE_combined] = max(TE_range(te_valid_combined_idx));

% Get all indices with max TE
max_TE_indices_combined = find(TE_range(te_valid_combined_idx) == max_TE_combined);

% Within these, find the minimum TR
[min_TR_combined, idx_min_TR_combined] = min(TR_range(tr_valid_combined_idx(max_TE_indices_combined)));

% Get the corresponding optimal TE, TR indices
te_opt_idx_combined = te_valid_combined_idx(max_TE_indices_combined(idx_min_TR_combined));
tr_opt_idx_combined = tr_valid_combined_idx(max_TE_indices_combined(idx_min_TR_combined));

% Extract SD(T2) values for all FA at optimal TE and TR
SD_combined_fa_values = squeeze(SD_T2_combined_valid(te_opt_idx_combined, tr_opt_idx_combined, :));

% Find the FA index that minimizes SD(T2)
[min_SD_combined, fa_opt_idx_combined] = nanmin(SD_combined_fa_values);

% Extract actual values
TE_opt_combined = TE_range(te_opt_idx_combined);
TR_opt_combined = TR_range(tr_opt_idx_combined);
FA_opt_combined = FA_range(fa_opt_idx_combined);

% Print results
fprintf('Optimal Parameters for Combined SD: TE = %.2f ms, TR = %.2f ms, FA = %.2f°\n', ...
    TE_opt_combined, TR_opt_combined, FA_opt_combined);
fprintf('Minimum combined SD is %d \n', round(min_SD_combined,5));
fprintf('Minimum Cartilage SD is %d \n', round(min_SD_cartilage,5));
fprintf('Minimum Meniscus SD is %d \n', round(min_SD_meniscus,5));
% figure; histogram(SD_T2_meniscus_valid(:), 50);
% title('Distribution of SD(T2) Meniscus');
% 
% figure; histogram(SD_T2_cartilage_valid(:), 50);
% title('Distribution of SD(T2) Cartilage');
% 
% figure; histogram(SD_T2_combined(:), 50);
% title('Distribution of Combined SD(T2)');

%% Plot Heatmap for Combined SD
% Define grid
[TE_grid, TR_grid] = meshgrid(TE_range, TR_range);
fa_slice_idx = round(length(FA_range) / 2); 
TE_max_grid = TR_grid / 2 - 4.31; % Maximum allowable TE from scanner hardware
invalid_region = TE_grid > TE_max_grid;

% Plot SD(T2) for Meniscus
figure;
hold on;
imagesc(TE_range, TR_range, SD_T2_meniscus_valid(:,:,fa_slice_idx)); 
set(gca, 'YDir', 'normal'); 
colormap('parula');
colorbar;
clim([0,50]);
xlabel('TE (ms)');
ylabel('TR (ms)');
title('SD of T2 for Meniscus');
contour(TE_grid, TR_grid, invalid_region, [0.5, 0.5], 'r', 'LineWidth', 2);
legend('Forbidden Region');
hold off;

% Plot SD(T2) for Cartilage
figure;
hold on;
imagesc(TE_range, TR_range, SD_T2_cartilage_valid(:,:,fa_slice_idx)); 
set(gca, 'YDir', 'normal'); 
colormap('parula');
colorbar;
clim([0,50]);
xlabel('TE (ms)');
ylabel('TR (ms)');
title('SD of T2 for Cartilage');
contour(TE_grid, TR_grid, invalid_region, [0.5, 0.5], 'r', 'LineWidth', 2);
legend('Forbidden Region');
hold off;

% Compute and Plot Combined SD(T2)
SD_T2_combined = sqrt(SD_T2_meniscus_valid.^2 + SD_T2_cartilage_valid.^2);
figure;
hold on;
imagesc(TE_range, TR_range, SD_T2_combined(:,:,fa_slice_idx)); 
set(gca, 'YDir', 'normal'); 
colormap('parula');
colorbar;
clim([0,20]);
xlabel('TE (ms)');
ylabel('TR (ms)');
pbaspect([1 1.2 1]);
title('Combined SD of T2 Estimation');
contour(TE_grid, TR_grid, invalid_region, [0.5, 0.5], 'r', 'LineWidth', 2);
plot(TE_opt_combined, TR_opt_combined, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(TE_opt_combined + 0.2, TR_opt_combined, 'Optimised Protocol', 'Color', 'red', 'FontSize', 9);
legend({'Forbidden Region (TE > TR/2 - 4.31)', 'Optimised Protocol'}, 'Location', 'northwest');
hold off;
set(gcf, 'MenuBar', 'none');
set(gcf, 'ToolBar', 'none');
exportgraphics(gcf, 'combined_sd_heatmap_final.png', 'Resolution', 300);

%% Number of Scans
SE = 1; % ms
N = (min_SD_combined / SE)^2;
fprintf('We need %d scans\n', round(N));


% Define a range of SE values (1ms to 4ms)
SE_values = linspace(1, 4, 20);

% Compute the required number of scans for each SE
N_values = (min_SD_combined ./ SE_values) .^ 2;

% Plot the results
figure;
plot(SE_values, N_values, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Standard Error (SE) [ms]');
ylabel('Required Number of Scans (N)');
title('Number of Scans Required vs. SE');
grid on;
set(gca, 'YScale', 'log'); % Log scale for better visualization
