temp_y = y; temp_u = u;
temp_data = iddata(temp_y,temp_u,Ts);
y_test = y(1401:1682,:);
u_test = u(1401:1682,:);
data_test = iddata(y_test,u_test,Ts);
y = y(1:1400,:);
u = u(1:1400,:);
data = iddata(y,u,Ts);

%udata = iddata([],u_test,Ts);

l = size(y);
nrows = l(1);
clear l

%5.
%Impulse
impulse_sys = impulseest(data);
plot(impulse_sys);
pause

%6
%Noise 1% 10%
e = 0.01*variance_y.*randn(nrows,2);
y_noise_1 = y + e;
e = 0.1*variance_y.*randn(nrows,2);
y_noise_2 = y + e;
hold on
plot(y_noise_1(:,1))
plot(y_noise_2(:,1))
plot(y(:,1))
legend('y1\_1% noise','y1\_10% noise','y1\_original')
pause

%Noise 1%
data_noise_1 = iddata(y_noise_1,u,Ts);  
impulse_sys_noise_1 = impulseest(data_noise_1);
plot(impulse_sys_noise_1);
pause
%Noise 10%
data_noise_2 = iddata(y_noise_2,u,Ts);  
impulse_sys_noise_2 = impulseest(data_noise_2);
plot(impulse_sys_noise_2);
pause

%7
e = randn(nrows,2);
colored_noise = filter(1,[1 -0.95],e);
plot(colored_noise,'r--')
title( 'colored noise')
legend('colored'); ylabel('v[t]');
pause
close all 
y_coloured=y+colored_noise ;
hold on
plot(y_coloured(:,1))
plot(y(:,1))
legend('Coloured Output y1','Original Output y1')
pause
close all
hold on
plot(y_coloured(:,2))
plot(y(:,2))
legend('Coloured Output y2','Original Output y2')
pause
close all

data_coloured = iddata(y_coloured,u,Ts);
impulse_sys_coloured = impulseest(data_coloured);
plot(impulse_sys_coloured)
pause

%8
%corr for impulse point estimation
%d = iddata(y(:,1),u(:,1),50)
%[ir,R,cl] = cra(d,20,46,2) 

%10 
ETFE_white = etfe(data_noise_1);
linearSystemAnalyzer({'bode'},ETFE_white)

ETFE_colour = etfe(data_coloured);
linearSystemAnalyzer({'bode'},ETFE_colour)