function [ g ] = G_H( H )
%由高度H计算重力加速度g

Re = 6378137;
g0 = 9.7803;

g = g0*(1 - 2*H/Re);
end
