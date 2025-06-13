% Load CNR data
cnr_data = readtable('Merged_CNR_Summary.csv');

% Only use "Femoral cartilage vs Meniscus"
cnr_filtered = cnr_data(strcmp(cnr_data.Comparison, 'Femoral cartilage vs Meniscus'), :);
subjects = unique(cnr_filtered.Subject);

figure;
hold on;
for i = 1:length(subjects)
    orig = cnr_filtered.CNR(strcmp(cnr_filtered.Subject, subjects{i}) & strcmp(cnr_filtered.Protocol, 'Original'));
    opt  = cnr_filtered.CNR(strcmp(cnr_filtered.Subject, subjects{i}) & strcmp(cnr_filtered.Protocol, 'Optimised'));
    plot([1 2], [orig opt], '-o', 'LineWidth', 2);
end
xlim([0.8 2.2]);
xticks([1 2]);
xticklabels({'Original','Optimised'});
ylabel('CNR');
title('CNR: Femoral cartilage vs Meniscus (Paired Comparison)');
saveas(gcf, 'cnr_comparison.png');