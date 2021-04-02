%% parameters, functions at end of code
fs = 100;     %sampling frequency
dx = 1 / fs;  %change in x
N = 1000;      %number of samples
x = (-N/2:N/2)*dx;  %x vector for a gaussian curve
fmax = fs / 2;      %nyquist limit
G = fft(dist(x)) * dx;    %discrete fourier transform, *dx for scaling
freq = (0:N)/(N * dx) - fmax;  %frequency resolution vector
Gshift = fftshift(G);   %shifting right end of fft to left

%% plotting fft
f1 = figure;
hold on

plot(freq, abs(Gshift),'b') %fft plot
plot(freq, FT(freq),'r')    %analytical plot
xlabel('frequency')
ylabel('fft')
title('Comparison between discrete and analytical FT')
legend('fft','analytic')

hold off

%% log-linear plot
f2 = figure;
%both plots still overlap
semilogy(freq, abs(Gshift), freq, FT(freq))
xlabel('x')
ylabel('G (log scale)')
title('log-linear scale')
grid on

%% hanning window
% reduces spectral leakage for uneven cycles.
Ghann = fft(hann(length(G)).*dist(x));
GhannShift = fftshift(Ghann);

f3 = figure;
hold on
plot(freq, abs(GhannShift)/N) % 1/N for scaling
title('fft with hanning window')
xlabel('frequencies')
hold off

%% functions
% normalized gaussian function
function f = dist(x)
    a = 0.4;
    c = 1/(a * sqrt(2 * pi)); %normalization factor
    f = c * exp(-x.^2 ./ (2 * a.^2));
end

% analytic fourier transform
function F = FT(freq)
    a = 0.4;
    F = exp(-2*pi.^2*freq.^2*a.^2);
end
    
