% Homework 01. Problem 2. Analyses of cancellation
clear;
clc;

%% calculation of expression (a) and (b)

x = 10.^(-(1:14));

% calculate (a)
sec_x = 1.0 ./ cos(x);
tan_x = tan(x);
a_orig = (1.0 - sec_x) ./ (tan_x.^2);
a_stab = -cos(x) ./ (1.0 + cos(x));

% calculate (b)
b_orig = (1.0 - (1.0 - x).^3) ./ x;
b_stab = 3.0 - 3.0*x + x.^2;

% calculate relative errors
err_a = abs(a_orig - a_stab) ./ abs(a_stab);
err_b = abs(b_orig - b_stab) ./ abs(b_stab);


%% Plot of (a)

fig = figure('Position', [100, 100, 500, 180], 'Color', 'w');

subplot(1,2,1);
semilogx(x, a_orig, 'b.-', 'DisplayName', 'Original (a)', 'LineWidth', 1.5);
hold on;
box on;
semilogx(x, a_stab, 'r.-', 'DisplayName', 'Eq. (2.1)', 'LineWidth', 1.5);
ylim([-0.6, 0.1]);

set(gca, 'XDir', 'reverse');
xlabel('\it{x}', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('\it{y}', 'FontName', 'Times New Roman', 'FontSize', 12);
legend('Location','northwest', 'FontSize', 10);
legend('boxoff');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

subplot(1,2,2);
loglog(x, err_a, 'g.-', 'LineWidth', 1.5);

set(gca, 'XDir', 'reverse');
xlabel('\it{x}', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Relative error', 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

hold off;

% save the figure
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(fig, 'Cancellation_analysis_a.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
fprintf('Saved figure Cancellation_analysis_a.pdf\n');


%% Plot of (b)

fig = figure('Position', [100, 100, 500, 180], 'Color', 'w');

subplot(1,2,1);
semilogx(x, b_orig, 'b.-', 'DisplayName', 'Original (b)', 'LineWidth', 1.5);
hold on;
box on;
semilogx(x, b_stab, 'r.-', 'DisplayName', 'Eq. (2.2)', 'LineWidth', 1.5);
ylim([2.7, 3.1]);

set(gca, 'XDir', 'reverse');
xlabel('\it{x}', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('\it{y}', 'FontName', 'Times New Roman', 'FontSize', 12);
legend('Location','southeast', 'FontSize', 10);
legend('boxoff');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

subplot(1,2,2);
loglog(x, err_b, 'g.-', 'LineWidth', 1.5);
set(gca, 'YTick', logspace(-15, 0, 4));

set(gca, 'XDir', 'reverse');
xlabel('\it{x}', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Relative error', 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

hold off;

% save the figure
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(fig, 'Cancellation_analysis_b.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
fprintf('Saved figure Cancellation_analysis_b.pdf\n');