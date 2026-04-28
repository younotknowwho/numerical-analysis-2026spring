%% Hilbert matrix experiments, main program
% Firstly, run this program; and then, run `postprocess.m`

clear; 
clc; 
close all;

%% Set output directory
outdir = 'output';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% Define necessary variables
m = 19;
n_list = 2:m+1;  % n=2,3,4,...,20

cond_inf      = zeros(m,1);  % cond_{\infty}(H_n)
res_inf       = zeros(m,1);  % \|r_n\|_{\infty}
rel_res_inf   = zeros(m,1);  % \|r_n\|_{\infty} / \|b\|_{\infty}
err_inf       = zeros(m,1);  % \|\bar{x} - x\|_{\infty}
rel_err_inf   = zeros(m,1);  % \|\bar{x} - x\|_{\infty} / \|x\|_{\infty}
min_sig_digit = zeros(m,1);  % minimum significant figures of \bar{x}
max_sig_digit = zeros(m,1);  % maximum significant figures of \bar{x}
avg_sig_digit = zeros(m,1);  % average significant figures of \bar{x}
threshold_n   = NaN; % min. n for \|\bar{x} - x\|_{\infty}>=\|x\|_{\infty}

all_digits = cell(m,1);
all_xbar   = cell(m,1);

%% Calculation for Hilbert Matrix Problems
for k = 1:m
    n = n_list(k);
    H = hilb(n);
    x = ones(n,1);
    b = H * x;

    % Solve H*x = b
    xbar = H \ b;

    % Calculate residual and error
    r = b - H * xbar;
    e = xbar - x;

    % Calculate condition number and norms
    cond_inf(k)    = cond(H, inf);
    res_inf(k)     = norm(r, inf);
    rel_res_inf(k) = norm(r, inf) / norm(b, inf);
    err_inf(k)     = norm(e, inf);
    rel_err_inf(k) = norm(e, inf) / norm(x, inf);

    % Calculate significant figures
    digits = effective_digits(xbar, x);
    min_sig_digit(k) = min(digits);
    max_sig_digit(k) = max(digits);
    avg_sig_digit(k) = mean(digits);

    all_digits{k} = digits;
    all_xbar{k}   = xbar;

    % Find possible threshold_n
    if isnan(threshold_n) && any(abs(e) >= abs(x))
        threshold_n = n;
    end
end

%% Save Data
save(fullfile(outdir, 'hilbert_workspace.mat'));
