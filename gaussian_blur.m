function y_blur = gaussian_blur(y, window_size, sigma)
    x = linspace(-window_size/2, window_size/2, window_size);
    ker = exp(-x.^2 / (2*sigma^2));
    ker = ker / sum(ker);
    y_blur = conv(y, ker, 'same');
end