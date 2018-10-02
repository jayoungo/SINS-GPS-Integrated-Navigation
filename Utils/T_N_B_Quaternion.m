function [ Tnb ] = T_N_B_Quaternion( Q )
%由四元数Q计算方向余弦矩阵Tnb

if(size(Q) ~= 4)
    error('Input Error! In function T_N_B_Quaternion.');
end
q_0 = Q(1)^2;
q_1 = Q(2)^2;
q_2 = Q(3)^2;
q_3 = Q(4)^2;
q_01 = Q(1)*Q(2);
q_02 = Q(1)*Q(3);
q_03 = Q(1)*Q(4);
q_12 = Q(2)*Q(3);
q_13 = Q(2)*Q(4);
q_23 = Q(3)*Q(4);

Tnb = [q_0 + q_1 - q_2 - q_3 2*(q_12 - q_03) 2*(q_13 + q_02);...
    2*(q_12 + q_03) q_0 - q_1 + q_2 - q_3 2*(q_23 - q_01);...
    2*(q_13 - q_02) 2*(q_23 + q_01) q_0 - q_1 - q_2 + q_3];
end
