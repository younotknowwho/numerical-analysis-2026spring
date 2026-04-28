function digits = effective_digits(xbar, x)
%Calculate effective significant digits
    error_abs = abs(xbar - x);
    digits = zeros(size(error_abs));

    k = log10(abs(xbar)) + 1;

    % If absolute error is zero
    zero_mask = (error_abs == 0);
    digits(zero_mask) = 15;   

    % If absolute error is not zero
    nz = ~zero_mask;
    digits(nz) = max(0, floor(k - log10(2 * error_abs(nz))));
end