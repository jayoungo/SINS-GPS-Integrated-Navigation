function [ F_new ] = R_K_2( deltaT, F, d_F_0, d_F )
%使用二级二阶龙格-库塔方法求解常微分方程初值问题的一步数值解（多处复用）
%参数：步长deltaT，初值F，微分初值d_F_0，微分d_F

F_new = F + deltaT/2*(d_F_0 + d_F);
end
