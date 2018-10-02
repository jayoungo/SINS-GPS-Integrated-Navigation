function [ omega ] = omega_matrix( w )
%计算四元数微分方程中的"ω"矩阵

if(size(w) ~= 3)
    error('Input Error! In function omega_matrix.');
end
x = w(1);
y = w(2);
z = w(3);

omega = [0 -x -y -z;x 0 z -y;y -z 0 x;z y -x 0];
end
