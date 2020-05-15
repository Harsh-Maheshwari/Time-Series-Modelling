% This file helps visualize white noise and colored noise
clear all;close all
display('Lets plot a white noise sequence of 1000 samples')
Ts =1; %sampling period
Fs = 1/Ts; %sampling frequency
fprintf(1,'\nLet us assume that sampling time is %g sec\n',Ts);
e = 0.1*randn(2048,1);
plot(e)
xlabel('Sample, t'); 
ylabel('e[t]');
title('White Noise sequence');
fprintf(1,'\n See very carefully, the sequence is uncorrelated in time\n');
display(' ')
display('Paused...press any key to continue...')
pause

% Assume a slowly decaying TF
fprintf(1,'\n We can get colured noise by filtering white noise\n')
num = 1/10;den = [1 -0.9];
display('Consider the following linear system')
mytf = tf(num,den,'Ts',1,'Variable','q^-1')
display(' ')
display('Paused...press any key to continue...')
pause

colored_noise = filter(num,den,e);
hold on
display('The coloured noise looks as follows')
display('Note the time correlation in the signal')
plot(colored_noise,'r--')
title('White noise, colored noise')
legend('white','colored'); ylabel('e[t], v[t]');

display(' ')
display('Paused...press any key to continue...')
pause

%Now lets plot autocorrelations of the two sequences
display('  ')
display('Now lets plot autocorrelations of the two sequences')
maxlag = 11;
Re = xcorr(e,maxlag,'unbiased');
Rcolored = xcorr(colored_noise,maxlag,'unbiased');
figure(2)
subplot(2,1,1)
stem(-maxlag:maxlag,Re),title('autocorrelation sequence')
legend('white'), xlabel('\tau'); ylabel('R_e(\tau)')
subplot(2,1,2)
stem(-maxlag:maxlag,Rcolored)
legend('colored')
 xlabel('\tau'); ylabel('R_v(\tau)')
display(' ')
display('Autocorrelation of white noise is theoretically a dirac delta')
display('Autocorrelation of colored noise is that of an AR(1) process')
display('Now lets plot the power spectral density')
display('We will use the SPA command')
display(' ')
display('Paused...press any key to continue...')
pause

% Now find the signal spectra using spa command
Hanning_window = 30;
fprintf(1,'\nThere is an interesting parameter called Hannings Window size used by SPA\n')
fprintf(1,'\nIn this work we set it to M = %g\n',Hanning_window)
iddatae = iddata([],e,Ts,'Domain','Time','InterSample','zoh',...
                        'Name','e','Notes','only white noise');
idfrdspa_e = spa(iddatae,Hanning_window);
phi_e = get(idfrdspa_e,'SpectrumData'); % this is a 3-d array
phi_e = squeeze(phi_e); % squeeze out the non-singleton elements
freqE = get(idfrdspa_e,'Frequency'); 
figure(3)
hold on
loglog(freqE,phi_e); title('Power Spectral Density of white noise')

iddatav = iddata([],colored_noise,Ts,'Domain','Time','InterSample','zoh',...
                        'Name','v','Notes','only coloured noise');
idfrdspa_v = spa(iddatav,Hanning_window);
phi_v = get(idfrdspa_v,'SpectrumData'); % this is a 3-d array
phi_v = squeeze(phi_v); % squeeze out the non-singleton elements
freqV = get(idfrdspa_v,'Frequency'); 
hold on
loglog(freqV,phi_v,'r'); title('One-sided Power Spectral Density')
legend('white','colored'); , xlabel('\omega'); ylabel('\Phi_e(\omega), \Phi_v(\omega)')
%set(gca,'YLIM',[-0.1 1]);  grid on
display('PSD of white noise is nearly flat (equal power at all frequencies)')
display('Moreover the flat power equals variance of the white noise')
display('PSD of coloured noise is concentrated at low frequencies')
display('This suggests the filter was a low pass filter')