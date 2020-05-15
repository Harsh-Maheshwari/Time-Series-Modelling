% Plot periodogram of yearly sunspots average number. The data reads from
% 1700 to 2016.
clear all;close all
load Sunspot_yearly_avg.txt -ascii
year = Sunspot_yearly_avg(:,1)
x = Sunspot_yearly_avg(:,2);
Fs = 1; %Frequency is one measurement per year
plot(year,x,'rx'); hold on, plot(year,x,'--')
% Let me first make the data zero mean 
x = x-mean(x)
N = length(x);

% Now take (unitary DFT) fft of x (do we need to do zero padding?)
X = fft(x)/sqrt(N);

X_Periodogram = X.*conj(X);  % same as (abs(X))^2

%Lets make the frequency axis ready. 
f = 0:1/N:0.5; % This is basically k/N where k goes from 0 to N-1 
F = f*Fs; %This is the frequency in Hertz
w = 2*pi*f;
figure; subplot(3,1,1); title('One sided Periodogram, Total Power')
plot(f,2*X_Periodogram(1:length(f)))
xlabel('Relative Frequency, f = F/Fs = k/N')
ylabel('|X_N(f)|^2')

subplot(3,1,2)
plot(F,2*X_Periodogram(1:length(f)))
xlabel('Digital Frequency, F = Fs*f, samples/year, Hertz')
ylabel('|X_N(F)|^2')

subplot(3,1,3); 
% In order to plot wrt angular frequency we need to divide the fft command
% based values further with 2*pi. 
plot(w,2*X_Periodogram(1:length(f))./(2*pi));
xlabel('Angular Digital Frequency, w =  2*pi*f, rad')
ylabel('|X_N(w)|^2')



% Lets use Matlab's PERIODOGRAM command
% Below plots from 0 to Fs/2
[pxx2,f1_mat] = periodogram(x,ones(size(x)),N,Fs); subplot(3,1,2), hold on, plot(f1_mat,pxx2,'r--')

% This one plots from 0 to pi
[pxx1,w1_mat] = periodogram(x,ones(size(x)),N); subplot(3,1,3), hold on, plot(w1_mat,pxx1,'r--')


