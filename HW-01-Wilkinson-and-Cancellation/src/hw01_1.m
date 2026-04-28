% Homework 01. Problem 1. Analyses of Wilkinson polynomial
clear;
clc;

%% 1.1 generate the roots and the corresponding polynomial
roots_exact = 1:20;
a = poly(roots_exact);  % a = [a_{20}, a_{19}, ..., a_{0}]

writematrix(a', 'a.txt');

%% 1.2 perturb the coefficients and plot the roots

% set the parameters
epsilons = [1e-10, 1e-6, 1e-3];
num_trials = 50;

rng(123456);    % seed

% create a figure
fig = figure('Position', [100, 100, 600, 180], 'Color', 'w');

for idx = 1:length(epsilons)
    eps = epsilons(idx);
    subplot(1, 3, idx);

    hold on; 
    box on;
    
    % plot the exact roots
    plot(roots_exact, zeros(1, 20), 'rx', 'MarkerSize', 6, ... 
        'LineWidth', 1.2, 'DisplayName', 'Exact');
    
    % plot the perturbed roots
    for t = 1:num_trials
        % generate the random numbers
        r = randn(1, 20);

        % perturb the coefficients
        a_perturbed = a;
        a_perturbed(1:end-1) = a(1:end-1) .* (1 + eps * r);
        
        % solve the perturbed roots
        roots_perturbed = roots(a_perturbed);

        % plot
        if t == 1
            plot(real(roots_perturbed), imag(roots_perturbed), ...
                'b.', 'MarkerSize', 4, 'DisplayName', 'Perturbed');
        else
            plot(real(roots_perturbed), imag(roots_perturbed), ...
                'b.', 'MarkerSize', 4, 'HandleVisibility', 'off');
        end
    end
    
    % figure style settings
    title(['\epsilon = 10^{', num2str(log10(eps)), '}'], 'FontName', ...
        'Times New Roman', 'FontSize', 14);
    xlabel('Real Part', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel('Imaginary Part', 'FontName', 'Times New Roman', 'FontSize', 12);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', ...
        1.0, 'TickDir', 'in');

    hold off;
end

% save the figure
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(fig, 'Wilkinson_Perturbation.pdf', 'ContentType', ...
    'vector', 'BackgroundColor', 'none');
fprintf('The figure has been saved.\n');
