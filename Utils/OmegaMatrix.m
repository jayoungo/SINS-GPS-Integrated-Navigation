function [ Omega ] = OmegaMatrix( w )
%计算角速度向量w对应的反对称矩阵Omega

if(size(w) ~= 3)
    error('Input Error! In function OmegaMatrix.');
end
x = w(1);
y = w(2);
z = w(3);

Omega = [0 -z y;z 0 -x;-y x 0];
end
