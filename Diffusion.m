clear; clc; close all;

%% Define Tissue Properties
T1_meniscus = 960;   % T1 for meniscus (ms)
T2_meniscus = 26.7;  % T2 for meniscus (ms)
D_meniscus  = 5e-2; % Diffusion coefficient (mm²/s)
Tau = 3.62; % ms
Gmax = 45mT
T1_cartilage = 900;  % T1 for cartilage (ms)
T2_cartilage = 39;   % T2 for cartilage (ms)
D_cartilage  = 0.75; % Diffusion coefficient (mm²/s)

%% Define MRI Parameters
TE_range = linspace(3, 6, 10);
TR_range = linspace(20, 34, 10);
FA_range = linspace(10, 50, 10);
b_values = linspace(0, 1000, 5); % b-values for diffusion weighting

num_iterations = 1000;

% Initialize matrices
SD_T2_meniscus = NaN(length(TE_range), length(TR_range), length(FA_range));
SD_T2_cartilage = NaN(length(TE_range), length(TR_range), length(FA_range));

%% Monte Carlo Simulation
for te_idx = 1:length(TE_range)
    for tr_idx = 1:length(TR_range)
        for fa_idx = 1:length(FA_range)
            for b_idx = 1:length(b_values)
                % MRI parameters
                TE_DESS = TE_range(te_idx);
                TR_DESS = TR_range(tr_idx);
                FA = deg2rad(FA_range(fa_idx));
                b = b_values(b_idx);
                
                % Simulate for both meniscus and cartilage
                for tissue = ["meniscus", "cartilage"]
                    if tissue == "meniscus"
                        T1 = T1_meniscus;
                        T2 = T2_meniscus;
                        D  = D_meniscus;
                    else
                        T1 = T1_cartilage;
                        T2 = T2_cartilage;
                        D  = D_cartilage;
                    end
                    
                    % Monte Carlo Simulation
                    T2_estimates = zeros(num_iterations, 1);
                    
                    for iter = 1:num_iterations
                        % Compute relaxation & diffusion factors
                        E1 = exp(-TR_DESS / T1);
                        E2 = exp(-TE_DESS / T2);
                        diffusion_term = exp(-b * D);
                        
                        % Compute r factor and S2/S1 ratio with diffusion
                        r = sqrt((1 - E2^2) / ((1 - E1 * cos(FA))^2 - E2^2 * (E1 - cos(FA))^2));
                        S2_S1_ratio = exp(-2 * (TR_DESS - TE_DESS) / T2 - b * D) * sin(FA/2)^2 * ...
                            ((1 + exp(-TR_DESS / T1 - TR_DESS * b * D)) / (1 - cos(FA) * exp(-TR_DESS / T1 - TR_DESS * b * D)));
                        
                        % Generate noisy signals
                        signal_FID = 100 * exp(-TE_DESS/T2) * diffusion_term; % Add diffusion decay
                        signal_SE  = signal_FID * S2_S1_ratio;
                        
                        noisy_FID = signal_FID + randn;
                        noisy_SE  = signal_SE + randn;
                        
                        % Estimate T2 using the ratio
                        R_meas = noisy_SE / noisy_FID;
                        if R_meas > 0
                            T2_estimates(iter) = -2 * (TR_DESS - TE_DESS) / log(R_meas);
                        else
                            T2_estimates(iter) = NaN;
                        end
                    end
                    
                    % Compute SD(T2)
                    if tissue == "meniscus"
                        SD_T2_meniscus(te_idx, tr_idx, fa_idx, b_idx) = nanstd(T2_estimates);
                    else
                        SD_T2_cartilage(te_idx, tr_idx, fa_idx, b_idx) = nanstd(T2_estimates);
                    end
                end
            end
        end
    end
end

%% Visualization: Effect of Diffusion on T2 Estimation
figure;
hold on;
plot(b_values, squeeze(SD_T2_meniscus(5,5,5,:)), 'r-o', 'DisplayName', 'Meniscus');
plot(b_values, squeeze(SD_T2_cartilage(5,5,5,:)), 'b-o', 'DisplayName', 'Cartilage');
xlabel('b-value (s/mm²)');
ylabel('SD of T2 Estimation (ms)');
title('Effect of Diffusion on T2 Estimation in qDESS');
legend;
hold off;
