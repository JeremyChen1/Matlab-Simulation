data = readtable('Merged_T2_Stats.csv', 'VariableNamingRule', 'preserve');

subjects = unique(data.Subject);
tissues = {'Femoral cartilage', 'Tibial cartilage', 'Meniscus'};
colors = lines(length(subjects));  % 每个 HV 一个颜色

figure;
hold on;

for s = 1:length(subjects)
    for t = 1:length(tissues)
        idx_ori = strcmp(data.Subject, subjects{s}) & strcmp(data.Tissue, tissues{t}) & strcmp(data.Protocol, 'Original');
        idx_opt = strcmp(data.Subject, subjects{s}) & strcmp(data.Tissue, tissues{t}) & strcmp(data.Protocol, 'Optimised');

        orig_val = data{idx_ori, 'SD(T2) (ms)'};
        opt_val = data{idx_opt, 'SD(T2) (ms)'};

        % x 轴位置做 slight offset 避免重叠
        x = [t - 0.1, t + 0.1];

        plot(x, [orig_val, opt_val], 'o-', 'Color', colors(s,:), ...
            'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', subjects{s});
    end
end

xlim([0.5, length(tissues)+0.5]);
xticks(1:length(tissues));
xticklabels(tissues);
ylabel('SD(T2) (ms)');
legend('Location', 'northwest');
title('T2 SD paired comparison (Original vs Optimised)');

saveas(gcf, 'T2_SD_perSubject.png');
