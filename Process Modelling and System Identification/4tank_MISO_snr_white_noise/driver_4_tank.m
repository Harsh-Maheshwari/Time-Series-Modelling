% 1) This file will simulate the 4-tank system for a PRBS input
% 2) Linearize the 4-tank system and simulate the linearized system
% 3) Identify time-series models and compare performances
 
clear all;close all
global A1 A2 A3 A4 a1 a2 a3 a4 kc g gamma1 gamma2 k1 k2
global A B L C
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextFontSize',16);

%Model parameters
A1 = 28; A2 = 32; A3= 28; A4 = 32; a1 = 0.071; a2 = 0.057; a3 =0.071;
a4 = 0.057; kc =0.5; g= 981; gamma1 = 0.7; gamma2 = 0.6; k1 = 3.33;
k2 = 3.35;

%Simulation parameters
T0 = 0; % initial time
Ts = 1;
tspan = [T0:Ts];
% xhat0 = [2 2 2 2]'; %initial conditions for estimates
% h0 = [h0;xhat0]; %initial conditions for combined plant and estimated states
% v =[3;3];
Npoints = 1023;
u = idinput([Npoints 2 1],'PRBS',[0 1/20],[2 4]); %Full length PRBS
myoptions = odeset('NonNegative',[1 2 3 4]);
h0 = [12.2630   12.7832    1.6339    1.4090]
% [tout,hout] = ode45(@(t,h) four_tank(t,h,[3 3]),[0 400],h0,myoptions);

v = u;
v(:,1) = u(:,1);
v(:,2) = u(:,2);

vss_at_3volts = 3.*ones(size(v));
hss_at_3volts = [12.2630   12.7832    1.6339    1.4090];
h0 = hss_at_3volts;[0;0;0;0]';
hout = h0; 
tout = 0;

% We would like to corrupt measurements with noise at a given SNR
% So first we will need to find the variance of the true signal. So lets
% first simulate the noise free system response
for i = 1:Npoints
    tspan = [(i-1)*Ts i*Ts];
    [tout1,hout1] = ode15s(@(t,h) four_tank(t,h,v(i,:)),tspan,hout(end,:),myoptions);
    tout = [tout;i*Ts];hout = [hout;hout1(end,:)];
end

subplot(2,1,1)
plot(tout(1:end-1),v(:,1));hold on;title('Input Voltage to Pump 1');
xlabel('Time, [seconds]');ylabel('V_1 [Volts]')
subplot(2,1,2)
plot(tout(1:end-1),v(:,2));hold on;title('Input Voltage to Pump 2');
xlabel('Time, [seconds]');ylabel('V_2 [Volts]')

figure(2)
subplot(2,1,1)
plot(tout,hout(:,1));hold on;
xlabel('Time, [seconds]');ylabel('h_1 [cm]')
subplot(2,1,2)
plot(tout,hout(:,2));hold on;
xlabel('Time, [seconds]');ylabel('h_2 [cm]')

figure(3)
subplot(2,1,1)
plot(tout,hout(:,3));hold on;
xlabel('Time, [seconds]');ylabel('h_3 [cm]')
subplot(2,1,2)
plot(tout,hout(:,4));hold on;
xlabel('Time, [seconds]');ylabel('h_4 [cm]')


%Initialize noise states
snr = 10; %signal to noise ratio
signal_var = var(hout(:,1:2)); % Signal variance
noise_var = signal_var./10; %Noise variance
n1=1; n2=1;
hmeas = hout(1,[1 2]) + sqrt(noise_var).*randn(size(noise_var));

% Let us corrupt signal with noise at the SNR level
vk = repmat(sqrt(noise_var),size(hout,1),1).*randn(size(hout(:,[1 2])));
hmeas = hout(:,[1 2]) + vk;

% sol = ode45(@(t,h) four_tank(t,h,v),tspan,h0,[]);
% yout = hout(:,[1 2]) + sqrt(0.01)*randn(size(hout(:,1:2)));

% Plot results
figure(2)
subplot(2,1,1)
plot(tout,hmeas(:,1),'r--')
subplot(2,1,2)
plot(tout,hmeas(:,2),'r--')
% fprintf(1,'\nPaused....\n')
% pause
% fprintf(1,'\nResumed....\n')

% Linearize 4-tank system
% Operating point of linearization
xss = [12.2630   12.7832    1.6339    1.4090];
uss = [3;3];
[A,B,C,D] = linmod('four_tank_sfun',xss,uss)

lin_4tank_sys = ss(A,B,C,D)

[hmeas_lin_dev,t,hlin_dev] = lsim(lin_4tank_sys,v-vss_at_3volts,tout(1:end-1),h0-hss_at_3volts);

%Convert deviation variables to engineering form
hlin = hlin_dev + repmat(hss_at_3volts,size(hlin_dev,1),1);

% Plot linear results and compare with nonlinear results
figure(2)
% subplot(2,1,1); plot(t,hlin(:,1),'g');
% subplot(2,1,2); plot(t,hlin(:,2),'g');
figure(3),title('Comparison of linear and nonlinear models')
subplot(2,1,1); plot(t,hlin(:,3),'g');
subplot(2,1,2); plot(t,hlin(:,4),'g');

% Lets use the input output data to generate time-series models
sampleT = 1;
hmeas_dev=hmeas-repmat(hss_at_3volts(1:2),size(hmeas,1),1);
iddatayu = iddata(hmeas_dev(1:end-1,:),v-vss_at_3volts,sampleT,'Domain','Time','InterSample',{'zoh';'zoh'},...
                  'OutputName',{'h1','h2'}','InputName',{'v1bar','v2bar'},'Notes','input-output data');
MARX_MV = arx(iddatayu,'na',[2 2;2 2],'nb',[2 2;2 2],'nk',[1 1;1 1]);
MARX_high_order_MV = arx(iddatayu,'na',[10 10;10 10],'nb',[10 10;10 10],'nk',[1 1;1 1]);

y_arx_dev = sim(MARX_MV,v-vss_at_3volts);
y_arx = y_arx_dev + repmat(hss_at_3volts(1:2),size(y_arx_dev,1),1);
figure(2)
subplot(2,1,1);plot(1:Npoints,y_arx(:,1),'m--')
subplot(2,1,2);plot(1:Npoints,y_arx(:,2),'m--')
y_arx_high_order_dev = sim(MARX_high_order_MV,v-vss_at_3volts);
y_arx_high_order = y_arx_high_order_dev + repmat(hss_at_3volts(1:2),size(y_arx_high_order_dev,1),1);
figure(2)
subplot(2,1,1);plot(1:Npoints,y_arx_high_order(:,1),'k--')
subplot(2,1,2);plot(1:Npoints,y_arx_high_order(:,2),'k--')
legend('Signal','Meas','ARX_MV(2,2)','ARX_MV(10,10')

%Unlike ARX, ARMAX cannot handle MIMO data. It is restricted to MISO. So
%lets first geneate two MISO data series
%y1 v/s u1 and u2
iddatay1u = iddata(hmeas_dev(1:end-1,1),v-vss_at_3volts,sampleT,'Domain','Time','InterSample',{'zoh';'zoh'},...
                  'OutputName',{'h1'}','InputName',{'v1bar','v2bar'},'Notes','input-output data');
%y2 v/s u1 and u2              
iddatay2u = iddata(hmeas_dev(1:end-1,2),v-vss_at_3volts,sampleT,'Domain','Time','InterSample',{'zoh';'zoh'},...
                  'OutputName',{'h2'}','InputName',{'v1bar','v2bar'},'Notes','input-output data');
 
% MISO ARX for y1              
MARXy1u = arx(iddatay1u,'na',[2],'nb',[2 2],'nk',[1 1]);
% Lets simulate (note this is not prediction) 
y1_arx_dev = sim(MARXy1u,v-vss_at_3volts);
y1_arx = y1_arx_dev + repmat(hss_at_3volts(1),size(y1_arx_dev,1),1);

% MISO ARX for y2
MARXy2u = arx(iddatay2u,'na',[2],'nb',[2 2],'nk',[1 1]);
% Lets simulate (note this is not prediction) 
y2_arx_dev = sim(MARXy2u,v-vss_at_3volts);
y2_arx = y2_arx_dev + repmat(hss_at_3volts(2),size(y2_arx_dev,1),1);

% MISO ARX for y1
MARXy1u_high_order = arx(iddatay1u,'na',[10],'nb',[10 10],'nk',[1 1]);
% Lets simulate (note this is not prediction) 
y1_arx_high_order_dev = sim(MARXy1u_high_order,v-vss_at_3volts);
y1_arx_high_order = y1_arx_high_order_dev + repmat(hss_at_3volts(1),size(y1_arx_high_order_dev,1),1);

% MISO ARX for y2
MARXy2u_high_order = arx(iddatay2u,'na',[10],'nb',[10 10],'nk',[1 1]);
% Lets simulate (note this is not prediction) 
y2_arx_high_order_dev = sim(MARXy2u_high_order,v-vss_at_3volts);
y2_arx_high_order = y2_arx_high_order_dev + repmat(hss_at_3volts(2),size(y2_arx_high_order_dev,1),1);

figure(2)
subplot(2,1,1);plot(1:Npoints,y1_arx,'m')
subplot(2,1,2);plot(1:Npoints,y2_arx,'m')
subplot(2,1,1);plot(1:Npoints,y1_arx_high_order,'k')
subplot(2,1,2);plot(1:Npoints,y2_arx_high_order,'k')

%MISO ARMAX for y1
MARMAXy1u = armax(iddatay1u,'na',[2],'nb',[2 2],'nc',[1],'nk',[1 1]);
% Lets simulate (note this is not prediction) 
y1_armax_dev = sim(MARMAXy1u,v-vss_at_3volts);
y1_armax = y1_armax_dev + repmat(hss_at_3volts(1),size(y1_armax_dev,1),1);

% MISO ARMAX for y2
MARMAXy2u = armax(iddatay2u,'na',[2],'nb',[2 2],'nc',[1],'nk',[1 1]);
% Lets simulate (note this is not prediction) 
y2_armax_dev = sim(MARMAXy2u,v-vss_at_3volts);
y2_armax = y2_armax_dev + repmat(hss_at_3volts(2),size(y2_armax_dev,1),1);

figure(2)
subplot(2,1,1);plot(1:Npoints,y1_armax(:,1),'c')
subplot(2,1,2);plot(1:Npoints,y2_armax(:,1),'c')
title('Data v/s Model Structure Simulations')

% State space model
M_SS = n4sid(iddatayu)
y_ss_dev = sim(M_SS,v-vss_at_3volts);
y_ss = y_ss_dev + repmat(hss_at_3volts(1:2),size(y_arx_dev,1),1);
figure(2)
subplot(2,1,1);plot(1:Npoints,y_ss(:,1),'y')
subplot(2,1,2);plot(1:Npoints,y_ss(:,2),'y')
legend('Signal','Meas','ARX_MV(2,2)','ARX_MV(10,10)','ARX(2,2)','ARX(10,10)','ARMAX(2,2,2)','n4sid')

figure
compare(iddatay1u,MARXy1u,MARXy1u_high_order,MARMAXy1u)

figure
compare(iddatay2u,MARXy2u,MARXy2u_high_order,MARMAXy2u)

figure
resid(MARXy1u,iddatay1u,'corr')
figure
resid(MARXy2u,iddatay2u,'corr')
figure
resid(MARXy1u_high_order,iddatay1u,'corr')
figure
resid(MARXy2u_high_order,iddatay2u,'corr')
figure
resid(MARMAXy1u,iddatay1u,'corr')
figure
resid(MARMAXy2u,iddatay2u,'corr')

fpe_MARMAXy1u = fpe(MARMAXy1u)
fpe_MARMAXy2u = fpe(MARMAXy2u)              
fpe_MARX_MV = fpe(MARX_MV)             
fpe_MARX_high_order_MV = fpe(MARX_high_order_MV)              
  fpe_MARXy1u = fpe(MARXy1u)             
  fpe_MARXy1u_high_order = fpe(MARXy1u_high_order)            
  fpe_MARXy2u = fpe(MARXy2u)           
  fpe_MARXy2u_high_order = fpe(MARXy2u_high_order)               
  fpe_M_SS = fpe(M_SS)           
% Implementation of the Luenberger Observer
% % the input sequence
% u_t.time = [0:Ts:length(u)-1]';
% u_t.signals.values = v; % the two voltages
% u_t.signals.dimensions = 2;
% 
% %check for observability
% OB = obsv(A,C)
% 
% fprintf(1,'\nRank of observability matrix = %g\n',rank(OB));
% L = place(A',C',[-2;-3;-1+2i;-1-2i]);
% L = L';
% xhatminus = [2 2 2 2];
% xhatplus = []
% tspan_e = tspan
% xhat0 = [2 2 2 2];
% for i = 2:size(tout,1)
%     xhatminus(i,:) = (A*xhatminus(i-1,:)' + Bd*[0;0] + L*((yout(i-1,:)'-[12.3;12.8])-...
%                      Cd*xhatminus(i-1,:)' - Dd*[0;0]))';
%     xhatplus(i,:) =  (xhatminus(i,:)' + M*((yout(i,:)'-[12.3;12.8])- Cd*xhatminus(i,:)'))';
% end
% xhatplus = hout(:,5:8);
% xhatplus = xhatplus + repmat([12.3;12.8;1.64;1.41]',size(tout,1),1);
% figure(1);hold on;
% plot(tout,xhatplus(:,1),'k:')
% figure(2);hold on
% plot(tout,xhatplus(:,2),'k:')
% figure(3);hold on;plot(tout,xhatplus(:,3),'k:')
% figure(4);hold on;plot(tout,xhatplus(:,4),'k:')
% 
