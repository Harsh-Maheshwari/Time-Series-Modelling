function [hdot] = four_tank(t,h,v)
%v =[3;3]
global A1 A2 A3 A4 a1 a2 a3 a4 kc g gamma1 gamma2 k1 k2 
global A B L C
hdot = [-a1/A1*sqrt(2*g*h(1))+a3/A1*sqrt(2*g*h(3))+gamma1*k1/A1*v(1)
        -a2/A2*sqrt(2*g*h(2))+a4/A2*sqrt(2*g*h(4))+gamma2*k2/A2*v(2)
        -a3/A3*sqrt(2*g*h(3))+(1-gamma2)*k2/A3*v(2);
        -a4/A4*sqrt(2*g*h(4))+(1-gamma1)*k1/A4*v(1)];
%     dxhatbydt = A*h(5:8) + B*(v-[3;3]) + L*((h(1:2)+sqrt(0.01)*randn(2,1)-[12.3;12.8])-C*h(5:8));
    hdot = hdot ;%+ sqrt(0.1).*randn(4,1);
%     hdot = [hdot;dxhatbydt];
return    