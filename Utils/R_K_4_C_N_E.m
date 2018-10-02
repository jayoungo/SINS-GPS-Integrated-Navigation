function [ Cne_new ] = R_K_4_C_N_E( deltaT, Cne, Wnen_0, Wnen_1, Wnen_2 )
%使用四级四阶龙格-库塔方法求解常微分方程初值问题的一步数值解（位置矩阵微分方程）
%参数：步长deltaT，初值Cne，角速度：Wnen_0,Wnen_1,Wnen_2

k_1 = (-1).*OmegaMatrix(Wnen_0)*Cne;
k_2 = (-1).*OmegaMatrix(Wnen_1)*(Cne + deltaT/2*k_1);
k_3 = (-1).*OmegaMatrix(Wnen_1)*(Cne + deltaT/2*k_2);
k_4 = (-1).*OmegaMatrix(Wnen_2)*(Cne + deltaT*k_3);
Cne_new = Cne + deltaT/6*(k_1 + 2*k_2 + 2*k_3 + k_4);
end
