%% Postprocess for the calculation results
% Note: run `main.m` before you run this program
clear; 
clc; 
close all;

outdir = 'output';

% Load the data
load(fullfile(outdir, 'hilbert_workspace.mat'));

% Set figure size
fig_width = 320;
fig_height = 240;

%% Plot 1: condition number
fig1 = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);
semilogy(n_list, cond_inf, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on;
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\mathrm{cond}_{\infty}(H_n)$', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);
save_figure(fig1, outdir, 'fig_cond');

%% Plot 1b: log(cond_inf) and linear fit for n <= 13
fig1b = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);

logc = log(cond_inf);
plot(n_list, logc, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;

fit_mask = n_list <= min(13, n_list(end));
if nnz(fit_mask) >= 2
    % Linear fit and plot
    mdl = fitlm(n_list(fit_mask), logc(fit_mask));
    RSquare = mdl.Rsquared.Ordinary;    % R^2
    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);

    plot(n_list, slope*n_list + intercept, '--', 'LineWidth', 1.5);
    text('Interpreter', 'latex', 'String', ...
        sprintf('$f(n) = %.3f n %+ .3f$', slope, intercept), ...
         'Position', [2.5 50], 'FontSize', 12);
    text('Interpreter', 'latex', 'String', sprintf('$R^2 = %.5f$', RSquare), ...
         'Position', [2.5 43], 'FontSize', 12);
end
grid on; box on;
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\ln(\mathrm{cond}_{\infty}(H_n))$', 'FontSize', 12, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);
save_figure(fig1b, outdir, 'fig_logcond_fit');

%% Plot 2: relative residual and relative error
fig2 = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);
semilogy(n_list, res_inf, 's-', 'LineWidth', 1.8, 'MarkerSize', 6); 
hold on;
semilogy(n_list, err_inf, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on;
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\|r_n\|$ or $\|e_n\|$', 'FontSize', 12, 'Interpreter', 'latex');
legend({'$\|r_n\|$', '$\|e_n\|$'}, 'Location', 'northwest', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);
save_figure(fig2, outdir, 'fig_res_err');

%% Plot 3: heatmap of significant digits
max_n = max(n_list);
digit_mat = NaN(max_n, m);
for k = 1:m
    n = n_list(k);
    digit_mat(1:n, k) = all_digits{k};
end

fig3 = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);
h_img = imagesc(n_list, 1:max_n, digit_mat);
set(h_img, 'AlphaData', ~isnan(digit_mat));
set(gca, 'Color', 'w');

axis xy;
cb = colorbar;
cb.Label.String = 'Effective significant digits';
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('index', 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);
save_figure(fig3, outdir, 'fig_digits_heatmap');

%% Plot 3b: min/max/avg. significant digits
fig3b = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);
plot(n_list, min_sig_digit, 's-', 'LineWidth', 1.8, 'MarkerSize', 6); 
hold on;
plot(n_list, max_sig_digit, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);
plot(n_list, avg_sig_digit, 'x-', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on;
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Effective significant digits', 'FontSize', 12);
legend({'Min. significant digits', 'Max. significant digits', ...
    'Avg. significant digits'}, 'Location', 'best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);
save_figure(fig3b, outdir, 'fig_digits_min_max_avg');

%% Plot 4: \|(H^{-1})_j\|_1 for n = 8,10,12
fig4 = figure('Color', 'w', 'Position', [100,100,2.2*fig_width,fig_height]);

ax1 = subplot(1,2,1);   % left
ax2 = subplot(1,2,2);   % right

selected_n_list = [8 10 12];
for n=selected_n_list
    H = hilb(n);
    G = inv(H);     % the inverse of matrix H
    gamma = sum(abs(G), 2); % \|(H^{-1})_j\|_1 = sum_k |(H^{-1})_{jk}|

    e = all_xbar{n-1}- ones(n,1);  % the error
    r = H*e;
    e_bound = gamma * norm(r, inf);     % upper bound of the error

    % Plot
    c = get(ax1, 'ColorOrder');
    idx = find(selected_n_list == n);
    this_color = c(mod(idx-1, size(c,1)) + 1, :);

    semilogy(ax1, 1:n, gamma, 'o-', 'Color', this_color, ...
             'LineWidth', 1.5, 'MarkerSize', 5, ...
             'DisplayName', sprintf('$n = %d$', n));
    if (n == selected_n_list(1))
        hold(ax1, 'on');
    end

    semilogy(ax2, 1:n, abs(e), 'o-', 'Color', this_color, ...
             'LineWidth', 1.5, 'MarkerSize', 5, ...
             'DisplayName', sprintf('$n=%d$', n));
    if (n == selected_n_list(1))
        hold(ax2, 'on');
    end

    semilogy(ax2, 1:n, e_bound, '--', 'Color', this_color, ...
             'LineWidth', 1.5, 'DisplayName', sprintf('$n=%d$', n));
end

legend(ax1, 'Location', 'best', 'Interpreter', 'latex');
grid(ax1, 'on'); box(ax1, 'on');
xlabel(ax1, '$j$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel(ax1, '$\gamma_j=\|(H^{-1})_j\|_1$', 'FontSize', 12, 'Interpreter', 'latex');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

text(ax1, 0.06, 0.95, '(a)', 'Units', 'normalized', ...
     'FontName', 'Times New Roman', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

legend(ax2, 'Location', 'best', 'Interpreter', 'latex');
grid(ax2, 'on'); box(ax2, 'on');
xlabel(ax2, '$j$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel(ax2, '$|e_j| \ \mathrm{or}\ \gamma_j\|r\|_\infty$', 'FontSize', 12, ...
    'Interpreter', 'latex');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

text(ax2, 0.06, 0.95, '(b)', 'Units', 'normalized', ...
     'FontName', 'Times New Roman', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

save_figure(fig4, outdir, 'fig_H_inverse_norm-and-error_bound_components');
