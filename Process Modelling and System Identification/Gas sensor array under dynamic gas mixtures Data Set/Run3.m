% ARX
k = 36; %22;
arx_order = [k*[2,2;2,2],k*[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]];
arx_sys = arx(data,arx_order);
yp = predict(arx_sys,data_test,4);

hold on
plot(y_test(:,1))
plot(yp.y(:,1))
legend('original y1','predicted y1')
pause
close all

hold on
plot(y_test(:,2))
plot(yp.y(:,2))
legend('original y2','predicted y2')
pause
close all

%ARMAX
k = 18;%9;
armax_order = [k*[2,2;2,2],k*[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],k*[2;2],[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]];
armax_sys = armax(data,armax_order);
yp = predict(armax_sys,data_test,4);

hold on
plot(y_test(:,1))
plot(yp.y(:,1))
legend('original y1','predicted y1')
pause
close all

hold on
plot(y_test(:,2))
plot(yp.y(:,2))
legend('original y2','predicted y2')
pause
close all

%Non Linear ARX
k = 10; % Dosen't depend on order ???
nlarx_order = [k*[2,2;2,2],k*[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]];
nlarx_sys = nlarx(data,arx_order);

% State Space N4SID
% n4sid_sys = n4sid(data,'best');

% OE Model
k = 7;
oe_order = [k*[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],k*[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2;2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]];
oe_sys = oe(data,oe_order);
