clear; clc; close all;

% Define parameters
T1_bone = 246;
T1_cartilage = 900;
T1_meniscus = 960;

T2_bone = 0.04 * 1000;       % convert to ms
T2_cartilage = 0.039 * 1000;
T2_meniscus = 0.0267 * 1000;

TR = linspace(0, 1000, 200);  % ms
TE = linspace(0, 200, 200);   % ms
M0 = 100;

% T1 Recovery: Mz = M0 * (1 - exp(-TR/T1))
Mz_bone = M0 * (1 - exp(-TR / T1_bone));
Mz_cartilage = M0 * (1 - exp(-TR / T1_cartilage));
Mz_meniscus = M0 * (1 - exp(-TR / T1_meniscus));

% T2 Decay: Mxy = M0 * exp(-TE/T2)
Mxy_bone = M0 * exp(-TE / T2_bone);
Mxy_cartilage = M0 * exp(-TE / T2_cartilage);
Mxy_meniscus = M0 * exp(-TE / T2_meniscus);

% Plot
figure('Position', [100 100 1000 400]);

subplot(1, 2, 1);
plot(TR, Mz_bone, '--k', 'LineWidth', 1.5); hold on;
plot(TR, Mz_cartilage, '-b', 'LineWidth', 1.5);
plot(TR, Mz_meniscus, '-.g', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Longitudinal Magnetization (Mz)');
title('T1 Recovery');
legend('Bone', 'Cartilage', 'Meniscus');
grid on;

subplot(1, 2, 2);
plot(TE, Mxy_bone, '--k', 'LineWidth', 1.5); hold on;
plot(TE, Mxy_cartilage, '-b', 'LineWidth', 1.5);
plot(TE, Mxy_meniscus, '-.g', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Transverse Magnetization (Mxy)');
title('T2 Decay');
legend('Bone', 'Cartilage', 'Meniscus');
grid on;

sgtitle('T1 Recovery and T2 Decay Curves for Different Knee Tissues');
saveas(gcf, 't1_t2_recovery_decay.png');
