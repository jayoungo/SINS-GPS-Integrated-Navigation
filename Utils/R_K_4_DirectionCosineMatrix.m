function [ Tnb_new ] = R_K_4_DirectionCosineMatrix( deltaT, Tnb, Wbnb_0, Wbnb_1, Wbnb_2 )
%使用四级四阶龙格-库塔方法求解常微分方程初值问题的一步数值解（方向余弦矩阵微分方程）
%参数：步长deltaT，初值Tnb，角速度：Wbnb_0,Wbnb_1,Wbnb_2

k_1 = Tnb*OmegaMatrix(Wbnb_0);
k_2 = (Tnb + deltaT/2*k_1)*OmegaMatrix(Wbnb_1);
k_3 = (Tnb + deltaT/2*k_2)*OmegaMatrix(Wbnb_1);
k_4 = (Tnb + deltaT*k_3)*OmegaMatrix(Wbnb_2);
Tnb_new = Tnb + deltaT/6*(k_1 + 2*k_2 + 2*k_3 + k_4);
end
