% A file downloaded from the net. Credit: Unknown
%A common use of Fourier transforms is to find the frequency components of 
% a signal buried in a noisy time domain signal. Consider data sampled at 
% 1000 Hz. Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and 
% 120 Hz sinusoid of amplitude 1 and corrupt it with some zero-mean random 
% noise:
clear all;close all
Fs = 1000;                    % Sampling frequency, per sec
Ts = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*Ts;                % Time vector
display('Consider the discrete time sum of two sinusoids corrupted by random noise')
fprintf('\n y [t] = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t) + 2*randn\n')
fprintf('\nSampling frequency Fs = %g\n',Fs)
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
y = x + 2*randn(size(t))*1;     % Sinusoids plus noise
plot(Fs*t(1:50),y(1:50)) %time is multiplied by Fs to conver to ms
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')
display('Here is the time-domain plot. It is difficult to identify the frequency components by looking at the')
display('original signal. Converting to the frequency domain, the discrete Fourier')
display('transform of the noisy signal y is found by taking the fast Fourier transform (fft)')
display('')
display('Paused...please press any key to continue')
pause

% transform of the noisy signal y is found by taking the fast Fourier 
% transform (FFT):
% original signal. Converting to the frequency domain, the discrete Fourier 
% transform of the noisy signal y is found by taking the fast Fourier 
% transform (FFT)')

% It is difficult to identify the frequency components by looking at the 
% original signal. Converting to the frequency domain, the discrete Fourier 
% transform of the noisy signal y is found by taking the fast Fourier 
% transform (FFT):

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;

% the frequency axis is chosen as k/(Nfft*h) where k goes from 0 to NFFT/2.
% The frequency response is periodic after that. Here h is the sampling
% period, thus h = 1/Fs. Thus f = (0:512)/(Nfft*h) = Fs*(0:512)/Nfft
figure
f2 = Fs*(-0.5:1/NFFT:0.5-1/NFFT);
plot(f2,abs(fftshift(Y))); xlabel('Frequency (Hertz)'); 
ylabel('|Y(F)|'); title('Two sided spectrum of 0.7*sin(2*pi*50*t) + sin(2*pi*120*t)')

f = Fs/2*linspace(0,1,NFFT/2+1); %
display('Two sided spectrum: symmetric for real signals')
display('the frequency axis is chosen as k/(Nfft*h) where k goes from 0 to NFFT/2.')
display('The frequency response is periodic after that. Here h is the sampling')
display('period, thus h = 1/Fs. Thus f = (0:512)/(Nfft*h) = Fs*(0:512)/Nfft')
display(' ')


display('The main reason the amplitudes are not exactly at 0.7 and 1 is because of') 
display('the noise. Several executions of this code (including recomputation of y')
display('will produce different approximations to 0.7 and 1. The other reason is') 
display('that you have a finite length signal. Increasing L from 1000 to 10000 in') 
display('the example above will produce much better approximations on average.')
display('')
display('Next lets look at the one-sided spectrum')
display('Paused...please press any key to continue')
pause
% Plot single-sided amplitude spectrum.
figure
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t) of 0.7*sin(2*pi*50*t) + sin(2*pi*120*t)')
xlabel('Frequency (Hz)')
ylabel('|Y(F)|')
display('One sided spectrum:')


