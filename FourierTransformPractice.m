fs = 100;     %sampling frequency of 100 Hz
dx = 1 / fs;  %change in x
N = 100;      %number of samples
x = -N:dx:N;  %x vector(specific value of x)
G = 1/fs * fft(dist(x));    %discrete fourier transform
freq = (-N/dx:N/dx)/(length(G) * dx);   %frequency vector
Gshift = fftshift(G);   %shifting right end of fft to left

%plots
hold on

plot(freq, abs(Gshift),'b') %fft plot
plot(freq, FT(freq),'r')    %analytical plot
legend('fft','analytic')

hold off

% normalized gaussian function
function f = dist(x)
    a = 0.4;
    c = 1/(a * sqrt(2 * pi)); %normalization factor
    f = c * exp(-x.^2 ./ (2 * a.^2));
end

% analytic fourier transform
function F = FT(freq)
    a = 0.4;
    A = 1/(a * sqrt(2 * pi));
    F = A * exp(-a.^2 * freq.^2 / 2);
end
    