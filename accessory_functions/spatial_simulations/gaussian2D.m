function G = gaussian2D(rows, cols, center, sigma, height)
    % Create coordinate grids
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Gaussian function formula
    G = height * exp(-((X - center(1)).^2 + (Y - center(2)).^2) / (2 * sigma^2));
end
