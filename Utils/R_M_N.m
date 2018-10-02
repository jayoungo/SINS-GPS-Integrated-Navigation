function [ RM, RN ] = R_M_N( L )
%由纬度L计算地球子午面曲率半径RM和卯酉圈曲率半径RN

Re = 6378137;
f = 1/298.257;

RM = Re*(1 - 2*f + 3*f*(sin(L))^2);
RN = Re*(1 + f*(sin(L))^2);
end
