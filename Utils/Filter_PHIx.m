function [ Yk ] = Filter_PHIx( Xk_1, Xk, Yk_1 )
%对水平误差角PHIx进行一阶低通滤波（Ref：房建成.一种新的惯导系统静基座快速初始对准方法[J].北京航空航天大学学报,1999(06):728-731.）

T = 10e-3;
Wdc = 0.001256;
Wac = tan(Wdc*T/2)*2/T;

Yk = Wac*T/2/(1 + Wac*T/2)*(Xk + Xk_1) + (1 - Wac*T/2)/(1 + Wac*T/2)*Yk_1;
end
