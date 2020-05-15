% CL 625: Process modeling and identification
% Sharad Bhartiya, Instructor
%
% This file shows how to do the following
% 1. How to define an idmodel using idpoly
% 2. How to generate an input signal for identification using idinput and iddata. 
% 3. Use of sim command for model simulation (not prediction!)
% 4. How to generate empirical transfer function estimate using etfe
% 5. How to generate etfe using first principles definition using fft command, etfe =
%                                                             Yn(w)./Un(w)
% 6. Ghat using power spectrum and cross-spectrum phi_yu = Ghat*phi_u
fprintf(1,'\nThis file shows how to do the following: \n')
fprintf(1,'1. How to define an idmodel using idpoly\n')
fprintf(1,'2. How to generate an input signal for identification using idinput and iddata\n') 
fprintf(1,'3. Use of sim command for model simulation (not prediction!)\n')
fprintf(1,'4. How to generate empirical transfer function estimate using etfe\n')
fprintf(1,'5. How to generate etfe using first principles definition using fft command, etfe = Yn(w)./Un(w)\n')
fprintf(1,'6. Ghat using power spectrum and cross-spectrum phi_yu = Ghat*phi_u\n')
fprintf(1,'\n press any key to continue ....\n')
pause


% First the plant to generate data. This is created using the idpoly
% command as follows for a Box-Jenkins model Ay = (B/F)u + (C/D)e
% M = IDPOLY(A,B,C,D,F,NoiseVariance,Ts). 
% We will assume only an ARMAX model Ay = Bu + Ce and value of Ts=0.1
display(' ')
display(' First the plant to generate data. This is created using the idpoly')
display('command as follows for a Box-Jenkins model Ay = (B/F)u + (C/D)e')
display('M = IDPOLY(A,B,C,D,F,NoiseVariance,Ts)') 
display('We will assume only an ARMAX model Ay = Bu + Ce and value of Ts=0.1')
clear all;close all
% 1. How to define an idmodel using idpoly
echo on
sampleT = 0.1;
Apoly = [1 -1.5 0.7]; Bpoly = [0 0.8 -0.6]; Cpoly = [1 -0.7];
% Apoly = [1 -0.8]; Bpoly = [0 0.8]; Cpoly = [1 -0.7];
%Assume standard deviation to be 5% of gain (variance = std^2)
display('Assume standard deviation to be 5% of gain (variance = std^2)')
NoiseVariance = (0.05)^2*abs(sum(Bpoly)/sum(Apoly));
idm = idpoly(Apoly,Bpoly,Cpoly,[],[],NoiseVariance,sampleT);
echo off
display('Let us look at the step response of our system using step command')
figure(1);step(idm);[ystep,tstep,ysdstep] = step(idm);
fprintf(1,'\nSo this is how the step response looks like\n press any key to continue ....\n')
pause

% 2. How to generate an input signal for identification using idinput. 
% In order to generate data, we need an input to drive the model. The input
% u is defined using idinput
% U = IDINPUT(N,TYPE,BAND,LEVELS) where N = [P Nu M] gives a M*P-by-Nu 
% input, periodic with period P and with M periods.
display('To generate output from our model we need an input signal u')
display('How to generate an input signal for identification using idinput') 
display('In order to generate data, we need an input to drive the model. The input')
display('u is defined using idinput')
display('U = IDINPUT(N,TYPE,BAND,LEVELS) where N = [P Nu M] gives a M*P-by-Nu ')
display('input, periodic with period P and with M periods.')
if tstep(end)/sampleT > 25
B = 25*sampleT/tstep(end)
else 
    B = 1 %Sampling itself is slow then choose B high
end
fprintf(1,'\n The type of input signal critically depends on value of B.  Here we choose B = %g\n', B)
display('pause....press any key to continue')
pause
echo on
u = idinput(1000,'PRBS',[0 B],[-1 1]); 
%Scale u
figure(2); plot(u);title('Input to plant, a PRBS signal');
echo off
fprintf(1,'\nThis is a PRBS input which is very popular \n press any key to continue ....\n')
pause
% Let us define the iddata object using the iddata command
display('Let us define the iddata object using the iddata command')
% DAT = IDDATA(Y,U,Ts)
echo on
iddatau = iddata([],u,sampleT,'Domain','Time','InterSample','zoh',...
                        'Name','u1','Notes','only input data');
echo off
% 3. Use of sim command                    
% Let us simulate the model idm using input data iddatau using SIM. I need
% to check this but note that for input starting from u(0) the first
% element of the output is y(1) and not y(0)
display('Now we are ready to simulate our idmodel object IDM.')
display('Let us simulate the model idm using input data iddatau using SIM command')
echo on
iddatay = sim(idm,iddatau,'Noise'); %Note sim is used by various toolboxes differently
iddatay1 = sim(idm,iddatau); % No noise here...
% To plot output first extract it from the iddata object iddatay
y = get(iddatay,'Outputdata'); %or just y=iddata.Outputdata
y1 = get(iddatay1,'Outputdata'); %or just y=iddata.Outputdata
echo off
figure(3), plot(y); title('Plant output with and without noise')
hold on;plot(y1,'r--');legend('with noise','without noise')
fprintf(1,'\nThis is response to the PRBS input with noise (Gu + v) and without noise (Gu)\n press any key to continue ....\n')
pause

% 4. How to generate empirical transfer function estimate using etfe
%Object containing both input as well as output data. 
display('*************NON-PARAMETRIC SYSTEM IDENTIFICATION**************')
display('Let us say you had only the PRBS input and its response and you wanted to find out transfer function G')
display('One way is using the ETFE command: G = Yn(w)./Un(w)')
display('We first need to store input-output data in the iddata object and then use ETFE command..press any key to continue')
pause
echo on
iddatayu = iddata(y,u,sampleT,'Domain','Time','InterSample','zoh',...
                  'OutputName','y','InputName','u','Notes','input-output data');
idfrdme = etfe(iddatayu); % This generates G in frequency domain

%we could also use a Hamming window to smoothen the etfe
display('we could also use a Hamming window to smoothen the etfe')
M = min(length(y)/10,30)
idfrdme_smooth = etfe(iddatayu,M)
% Now we should extract out the etfe Ghat to plot it. To see properties
% of the idfrd object, 
% idprops idfrd
% Note that G is stored in ResponseData property
etfe_G = get(idfrdme,'ResponseData'); % this is a 3-d array
etfe_G = squeeze(etfe_G); % squeeze out the non-singleton elements
freqG = get(idfrdme,'Frequency'); % in rad/sec
echo off
% freqG = freqG./(2*pi); %converted to hertz
figure(4)
 bode(idfrdme), hold on, loglog(freqG,abs(etfe_G),'r--'),title('Raw ETFE')
 figure(5),bode(idfrdme_smooth,'r--') ; title('Smoothened ETFE') %compare non-smooth v/s smooth idfrd
 display('This is the raw ETFE and smooth ETFE..press any key to continue')
pause

%5. Now etfe using first principles definition using fft command 
 %Okay lets check if this is really based on first principles definition of
 % etfe being finite DFT of y divided by finite DFT of u
 fprintf(1,'It is great to know some software, but in class we used first principles definitions\n')
 display('It is worthwhile to see if these first definitions work')
 display('So let us use the fft command and find G = (DFT of y)/(DFT of u)')
% Now lets try to calculate the dft using the definition
% First lets take FFT of output y
echo on
NFFT = 2^nextpow2(length(y)); % Next power of 2 from length of y
ffty = fft(y,NFFT)./sqrt(length(y)); % using Ljung's definition

%Now lets take fft of input u
fftu = fft(u,NFFT)./sqrt(length(u));% using Ljung's definition
Ghat = ffty./fftu;

% corresponding frequency (in rad/sampling interval or rad/sec)
freqYdef = (2*pi/sampleT)*(1/length(ffty):1/length(ffty):0.5);
echo off
display('Note that we convert magnitude to dB by 20 log10 (abs(G))')
figure()

subplot(2,1,1);hold on; semilogx([freqYdef],20.*log10(abs(Ghat(1:512)))), xlabel('\Omega, rad/s')
ylabel('|G(j\Omega)|')
title('First principles based ETFE using fft')
fprintf(1,'\ncompare this first principles calculation with etfes obtained using etfe command..paused\n')
pause
% Can also use ffplot to draw frequency response from IDFRD object
% figure(6);ffplot(idfrdme); ; %Note this is in Hertz
% title('using ffplot command');hold on; semilogy(freqY/(2*pi),periodY,'r--')
% legend('using ffplot command','using semilogy')

echo on
% 6. Ghat using power spectrum and cross-spectrum phi_yu = Ghat*phi_u
% Transfer function estimation using power spectrum, phi_yu = Ghat*phi_u
fprintf(1,'\n Yet another way is to use the PSD which can be obtained using SPA command')
idfrdspa = spa(iddatayu);
phi_yu = get(idfrdspa,'ResponseData'); % this is a 3-d array
phi_yu = squeeze(phi_yu); % squeeze out the non-singleton elements
freqYU = get(idfrdspa,'Frequency');
idfrdspa_u = spa(iddatau);
phi_u = get(idfrdspa_u,'SpectrumData'); % this is a 3-d array
phi_u = squeeze(phi_u); % squeeze out the non-singleton elements
freqU = get(idfrdspa_u,'Frequency');
if sum(freqYU == freqU) == length(freqU)
    okay = 1
else
    okay = 0
    error('Freq of phi_yu and phi_u do not match')
end

Ghat = phi_yu./phi_u;

figure();subplot(2,1,1);hold on
loglog(freqYU,abs(Ghat),'-.');hold on;
legend('smooth etfe based','power spectrum based')

figure()
semilogx(freqYU,20*log10(abs(Ghat))); grid on

figure()
loglog(freqU,phi_u), title('Power spectrum of PRBS')



