%%
% a = [0, 5, 10, 15, 20];

% Long axis data
values_405_long = [0.828571429, 0.914285714, 1, 0.957142857, 0.785714286];
values_488_long = [0.788461538, 0.826923077, 1, 0.846153846, 0.711538462];
values_633_long = [0.865671642, 0.820895522, 1, 0.671641791, 0.567164179];

% Short axis data
values_405_short = [0.835164835, 0.945054945, 1, 0.956043956, 0.725274725];
values_488_short = [0.787878788, 0.893939394, 1, 0.787878788, 0.772727273];
values_633_short = [0.850574713, 0.781609195, 1, 0.67816092, 0.597701149];

% Plotting
figure;
% 405 nm
plot(a, values_405_long, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);  % Blue circles with solid line
hold on;
plot(a, values_405_short, 'b^--', 'LineWidth', 1.5, 'MarkerSize', 8);  % Blue triangles with dashed line

% 488 nm
plot(a, values_488_long, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);  % Red squares with solid line
plot(a, values_488_short, 'rv--', 'LineWidth', 1.5, 'MarkerSize', 8);  % Red triangles with dashed line

% 633 nm
plot(a, values_633_long, 'g^-', 'LineWidth', 1.5, 'MarkerSize', 8);  % Green triangles with solid line
plot(a, values_633_short, 'gd--', 'LineWidth', 1.5, 'MarkerSize', 8);  % Green diamonds with dashed line

% Add labels and title
xlabel('Position');
ylabel('Signal intensity');
title('Plot of signal intensity depending on z-position');

% Add legend
legend('405 Long', '405 Short', '488 Long', '488 Short', '633 Long', '633 Short', 'Location', 'best');

% Adjust axis properties
grid on;  % Show gridlines
ylim([0, 1.1]);  % Adjust y-axis limits if necessary

% Hold off to end plotting
hold off;
%%
% % Simple test plot
% pos = [0, 5, 10, 15, 20];
% values = [0.828571429, 0.914285714, 1, 0.957142857, 0.785714286];
% 
% figure;
% plot(pos, values, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
% xlabel('Position');
% ylabel('Signal intensity');
% title('Test Plot');
% grid on;
% ylim([0, 1.1]);