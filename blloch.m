% T1 and T2 Relaxation Bloch Vectors Visualization
figure('Position', [100, 100, 1200, 600]);

% T1 relaxation vectors (along z-axis)
subplot(1,2,1);
hold on; grid on; axis equal;
title('T1 Relaxation');
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1,1]); ylim([-1,1]); zlim([0,1]);
view([-45 20]);
for z = linspace(0, 1, 5)
    quiver3(0, 0, 0, 0, 0, z, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
end

% T2 relaxation vectors (in xy-plane, decaying x)
subplot(1,2,2);
hold on; grid on; axis equal;
title('T2 Relaxation');
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1,1]); ylim([-1,1]); zlim([0,1]);
view([-45 20]);
for x = linspace(1, 0, 5)
    quiver3(0, 0, 0, x, 0, 1, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
end
saveas(gcf, 't1_t2_bloch_vectors.png');
