%% Higher precision version
%% Hilbert matrix experiments, main program
% Note: This program may take longer time to finish running

clear; 
clc; 
close all;

digits(50);   % Set higher precision

%% Set output directory
outdir = 'output_higherPrecision';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% Define necessary variables
m = 49;
n_list = 2:m+1;  % n=2,3,4,...,20

cond_inf      = vpa(zeros(m,1));  % cond_{\infty}(H_n)
res_inf       = vpa(zeros(m,1));  % \|r_n\|_{\infty}
rel_res_inf   = vpa(zeros(m,1));  % \|r_n\|_{\infty} / \|b\|_{\infty}
err_inf       = vpa(zeros(m,1));  % \|\bar{x} - x\|_{\infty}
rel_err_inf   = vpa(zeros(m,1));  % \|\bar{x} - x\|_{\infty} / \|x\|_{\infty}
min_sig_digit = vpa(zeros(m,1));  % minimum significant figures of \bar{x}
max_sig_digit = vpa(zeros(m,1));  % maximum significant figures of \bar{x}
avg_sig_digit = vpa(zeros(m,1));  % average significant figures of \bar{x}
threshold_n   = [];               % min. n for \|\bar{x} - x\|_{\infty}>=\|x\|_{\infty}

all_digits = cell(m,1);
all_xbar   = cell(m,1);

%% Calculation for Hilbert Matrix Problems
for k = 1:m
    n = n_list(k);

    H = 1./(sym(1:n) + sym(1:n).' - 1);
    x = vpa(ones(n,1));
    b = H * x;

    % Solve H*x = b
    xbar = H \ b;

    % Calculate residual and error
    r = b - H * xbar;
    e = xbar - x;

    % Calculate condition number and norms
    cond_inf(k)    = norm(H, inf) * norm(inv(H), inf);
    res_inf(k)     = norm(r, inf);
    rel_res_inf(k) = norm(r, inf) / norm(b, inf);
    err_inf(k)     = norm(e, inf);
    rel_err_inf(k) = norm(e, inf) / norm(x, inf);

    % Calculate significant figures
    sig_digits = effective_digits(xbar, x);
    min_sig_digit(k) = min(sig_digits);
    max_sig_digit(k) = max(sig_digits);
    avg_sig_digit(k) = mean(sig_digits);

    all_digits{k} = sig_digits;
    all_xbar{k}   = xbar;

    % Find possible threshold_n
    if isempty(threshold_n) && any(abs(e) >= abs(x))
        threshold_n = n;
    end
end

%% Save Data
save(fullfile(outdir, 'hilbert_workspace_higherPrecision.mat'));

%% Visualization

fig_width = 320;
fig_height = 240;

%% Plot 1: condition number
fig1 = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);
semilogy(n_list, cond_inf, 's-', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; box on;
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\mathrm{cond}_{\infty}(H_n)$', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);
save_figure(fig1, outdir, 'fig_cond');

%% Plot 1b: linear fit
fig1b = figure('Color', 'w', 'Position', [100,100,fig_width,fig_height]);

semilogy(n_list, cond_inf, 's-', 'LineWidth', 1.8, 'MarkerSize', 6); 
hold on;

% Convert sym to numbers
cond_inf_num = zeros(m, 1);
for i = 1:m
    cond_inf_num(i) = vpa(cond_inf(i), 50);
end

% Linear fit and plot
mdl = fitlm(n_list, log(cond_inf_num));

RSquare = mdl.Rsquared.Ordinary;    % R^2
intercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);

semilogy(n_list, exp(slope*n_list + intercept), '--', 'LineWidth', 1.5);

grid on; box on;
xlabel('$n$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\ln(\mathrm{cond}_{\infty}(H_n))$', 'FontSize', 12, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

text(0.1, 0.65, ['$y=\exp(',sprintf('%.3f n %+ .3f', slope, intercept),')$'], ...
    'Units', 'normalized', ...
     'FontName', 'Times New Roman', 'FontSize', 12, ...
     'Interpreter', 'latex',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(0.1, 0.55, sprintf('$R^2 = %.7f$', RSquare), 'Units', 'normalized', ...
     'FontName', 'Times New Roman', 'FontSize', 12, ...
     'Interpreter', 'latex',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

save_figure(fig1b, outdir, 'fig_logcond_fit');

%% Plot 2: relative residual vs relative error
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

%% Plot 3b: min/max/avg. significant digits changing with n
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
fig4 = figure('Color', 'w', 'Position', [100,100,fig_width,2.2*fig_height]);

ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

selected_n_list = [10 20 30];
for n=selected_n_list
    H = 1./(sym(1:n) + sym(1:n).' - 1);
    G = inv(H);     % the inverse of matrix H
    gamma = sum(abs(G), 2); % \|(H^{-1})_j\|_1 = sum_k |(H^{-1})_{jk}|

    e = vpa(all_xbar{n-1}) - vpa(ones(n,1));  % calculation error
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


%% Local functions

% Effective significant digits
function sig_digits = effective_digits(xbar, x)
    error_abs = abs(xbar - x);
    sig_digits = vpa(zeros(size(error_abs)));

    k = log10(abs(xbar)) + 1;

    % If absolute error is zero
    zero_mask = (error_abs == 0);
    current_digits = digits;   
    sig_digits(zero_mask) = current_digits;  

    % If absolute error is not zero
    nz = ~zero_mask;
    sig_digits(nz) = max(vpa(0), floor(k - log10(2 * error_abs(nz))));
end

% Save figure in multiple formats
function save_figure(fig, outdir, basename)
    pngfile = fullfile(outdir, [basename, '.png']);
    pdffile = fullfile(outdir, [basename, '.pdf']);
    figfile = fullfile(outdir, [basename, '.fig']);

    try
        exportgraphics(fig, pngfile, 'Resolution', 300);
        exportgraphics(fig, pdffile, 'ContentType', 'vector');
    catch
        print(fig, pngfile, '-dpng', '-r300');
        print(fig, pdffile, '-dpdf', '-vector');
    end

    savefig(fig, figfile);
end