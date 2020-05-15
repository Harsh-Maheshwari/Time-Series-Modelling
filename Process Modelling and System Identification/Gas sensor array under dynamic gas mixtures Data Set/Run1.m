T = table2array(readtable("csv/CO_per25Sec.csv"));
y = T(3:1684,2:3); %% Removing first few observations because they tend to be highly unsteady
u = T(3:1684,4:19);
t = T(3:1684,1);
clearvars T
Ts = 25;

plot(y(:,1))
title( 'y1 vs time')
pause

plot(y(:,2))
title( 'y2 vs time')
pause

plot(u)
title('U vs time')
pause

%EDA%

% Corr
plot(xcorr(y(:,1)))
title( 'Auto-correlation of y1')
legend('y1')
pause

plot(xcorr(y(:,2)))
title( 'Auto-correlation of y2')
legend('y2')
pause

plot(xcorr(y(:,1),y(:,2)))
title( 'Cross-correlation of y1 and y2')
legend('y1 x y2')
pause

% Mean Variables
y_mean = y-mean(y);
plot(y_mean)
title('y-Mean\_y vs time')
pause

u_mean = u-mean(u);
plot(u_mean)
title('u-Mean\_u vs time')
pause

% Corr using means :
plot(xcorr(y_mean(:,1)))
title( 'Auto-correlation of y1 mean')
legend('y1 mean')
pause

plot(xcorr(y_mean(:,2)))
title( 'Auto-correlation of y2 mean')
legend('y2 mean')
pause

plot(xcorr(y_mean(:,1),y_mean(:,2)))
title( 'Cross-correlation of y1 mean and y2 mean')
legend('y1 mean x y2 mean')
pause

% PSDs  Power-spectral-density Analysis
py1 = periodogram(y_mean(:,1)); 
plot(py1)
title( 'Periodogram of y1 Mean deviation')
pause

py2 = periodogram(y_mean(:,2));
plot(py2)
title( 'Periodogram of y2 Mean deviation')
pause

pu = periodogram(u_mean);
plot(pu)
title( 'Periodogram of U(input) Mean deviation')
pause

%%% Now we will always use the mean deviation variables
y = y_mean;
u = u_mean;
clear y_mean u_mean
variance_u = var(u);
variance_y = var(y);

%DFT
% y_diff1 = y(2:840,:)-y(1:839,:);              % Diff -1
% y_diff2 = y_diff1(2:839,:)-y_diff1(1:838,:);  % Diff -2
% y_diff3 = y_diff2(2:838,:)-y_diff2(1:837,:);  % Diff -3

% DFT y1
dft_y1 = fft(y(:,1));                               % Compute DFT y1
m = abs(dft_y1);                                    % Magnitude
dft_y1(m<1e-6) = 0;
p = unwrap(angle(dft_y1));                          % Phase
f = (0:length(dft_y1)-1)*100/length(dft_y1);        % Frequency vector
subplot(2,1,1)
plot(f,m)
title('Magnitude from DFT of y1')
subplot(2,1,2)
plot(f,p*180/pi)
title('Phase from DFT of y1')
pause

% DFT y2
dft_y2 = fft(y(:,2));                               % Compute DFT y2
m = abs(dft_y2);                                    % Magnitude
dft_y2(m<1e-6) = 0;
p = unwrap(angle(dft_y2));                          % Phase
f = (0:length(dft_y2)-1)*100/length(dft_y2);        % Frequency vector
subplot(2,1,1)
plot(f,m)
title('Magnitude from DFT of y2')
subplot(2,1,2)
plot(f,p*180/pi)
title('Phase from DFT of y2')

clear f m p pu py1 py2 dft_y1 dft_y2

