% CL 625: Process modeling and identification
% Sharad Bhartiya, Instructor
%
% The intention of this file is to show that if we have colored noise and we
% use simply y=Gu + e that is think of v as white then estimate of G is
% biased
%
% This file shows how to do the following
% 1. How to define an idmodel using idpoly
% 2. How to generate an input signal for identification using idinput and iddata. 
% 3. Use of sim command
% 4. Estimate parameters of different model structures

% First the plant to generate data. This is created using the idpoly
% command as follows for a general model Ay = (B/F)u + (C/D)e
% M = IDPOLY(A,B,C,D,F,NoiseVariance,Ts). 
% We will assume only an ARMAX model Ay = Bu + Ce and default value of Ts=1
clear all;close all
% 2. How to define an idmodel using idpoly
sampleT = 0.1;
Apoly = [1 -1.5 0.7]; Bpoly = [0 0.8 -0.6]; Cpoly = [1 -0.7];
% Apoly = [1 -0.8]; Bpoly = [0 0.8]; Cpoly = [1 -0.7];
%Assume standard deviation to be 10% of gain (variance = std^2)
NoiseVariance = (0.1)^2*abs(sum(Bpoly)/sum(Apoly));
idm = idpoly(Apoly,Bpoly,Cpoly,[],[],NoiseVariance,sampleT);
figure(1);subplot(2,1,1);step(idm);[ystep,tstep,ysdstep] = step(idm);

%3. Generate step response 
display('Let us look at the step response of our system using step command')
figure(1);step(idm);[ystep,tstep,ysdstep] = step(idm);

%4. How to implement idinput
%Lets also plot the bode plot. Note that for discrete time systems, bode
%plots only up to frequency pi (for digital angular frequency in radians)
%which it converts to some analog kind of frequency by dividing by Ts. So
%the x-axis is w/Ts. So to get the bandwidth, multiply the x-axis value by
%Ts to get bandwidth in radians
subplot(2,1,2);bode(idm)

BW = bandwidth(idm); %in w/Ts units

% So convert to digital frequency (F/Fs) as follows: w = 2*pi*F/Fs
BW_digital_relative = BW*sampleT/(2*pi) %This is F/Fs,

% As given in Ljung, choose B = 2.5 time BW_digital_relative
% How to generate an input signal for identification using idinput. 
% In order to generate data, we need an input to drive the model. The input
% u is defined using idinput
% U = IDINPUT(N,TYPE,BAND,LEVELS) where N = [P Nu M] gives a M*P-by-Nu 
% input, periodic with period P and with M periods.
%B = min(1,2.5*BW_digital_relative); % this is too high
B = BW_digital_relative;
% if tstep(end)/sampleT > 25
% B = 25*sampleT/tstep(end)
% else 
%     B = 1 %Sampling itself is slow then choose B high
% end
% %B=0.0625
u = idinput(1024,'PRBS',[0 B],[-1 1]); 
%Scale u
figure(2); plot(u);title('Input to plant');

%5. How to save u (or any time series) or iddata that is u and y as an
%iddata object
% Let us define the iddata object using the iddata command
% DAT = IDDATA(Y,U,Ts)
iddatau = iddata([],u,sampleT,'Domain','Time','InterSample','zoh',...
                        'Name','u1','Notes','only input data');

% 6. Use of sim command                    
% Let us simulate the model idm using input data iddatau using sim. I need
% to check this but note that for input starting from u(0) the first
% element of the output is y(1) and not y(0)
iddatay = sim(idm,iddatau,'Noise'); %Note sim is used by various toolboxes differently
iddatay1 = sim(idm,iddatau); % No noise here...

% To plot output first extract it from the iddata object iddatay
y = get(iddatay,'Outputdata'); %or just y=iddata.Outputdata
y1 = get(iddatay1,'Outputdata'); %or just y=iddata.Outputdata
figure(3), plot(y); title('Plant output with and without noise')
hold on;plot(y1,'r--');legend('with noise','without noise')

% 7. Now generate iddata object using both input and output data
iddatayu = iddata(y,u,sampleT,'Domain','Time','InterSample','zoh',...
                  'OutputName','y','InputName','u','Notes','input-output data');
% 8. Use iddata to generate model parameters (ARX, ARMAX)              
M_ARX = arx(iddatayu,'na',2,'nb',2,'nk',1);
M_ARX_high_order = arx(iddatayu,'na',10,'nb',10,'nk',1);
M_ARMAX = armax(iddatayu,'na',2,'nb',2,'nc',1,'nk',1);

%9. k-step ahead predictions using predict
k=2;
figure(4)
predict(M_ARX,'r',M_ARX_high_order,'k--',M_ARMAX,'b--',iddatayu,k)

legend('M_ARX','M_ARX_high_order', 'M_ARMAX')

%10. Residual analysis using resid command
figure(5)
resid(M_ARX,iddatayu,'corr')
figure(6)
resid(M_ARX_high_order,iddatayu,'corr')
figure(7)
resid(M_ARMAX,iddatayu,'corr')

% 11. Step response;
% figure;step(M_ARX); hold on; step(M_ARX_high_order), step(M_ARMAX); 
% legend('M_ARX','M_ARX_high_order', 'M_ARMAX')

%12 Frequency response comparison using polar corordinates (Bode plot)
figure(8)
[mag,phase,w] = bode(idm);
loglog(squeeze(w),squeeze(mag))
hold on
[mag,phase,w] = bode(M_ARX);
loglog(squeeze(w),squeeze(mag),'r--')
[mag,phase,w] = bode(M_ARX_high_order);
loglog(squeeze(w),squeeze(mag),'k-.')
hold on
[mag,phase,w] = bode(M_ARMAX);
loglog(squeeze(w),squeeze(mag),'g:')
legend('plant_armax(2,2,1)','arx(2,2)','arx(10,10)','armax(2,2,1)')

figure(9)
[y,t] = step(idm); plot(t,y);hold on; 
[y,t] = step(M_ARX); plot(t,y,'r--');
[y,t] = step(M_ARX_high_order); plot(t,y,'k-.')
[y,t] = step(M_ARMAX);plot(t,y,'g:')
legend('plant_armax(2,2,1)','arx(2,2)','arx(10,10)','armax(2,2,1)')

%To evaluate the best model using %fit
figure(10)
compare(iddatayu,M_ARX,'rx',M_ARX_high_order,'g--',M_ARMAX,'ks',1)

%To evaluate the impulse response coeffs using correlation analysis using
%data
figure(11)
cra(iddatayu)