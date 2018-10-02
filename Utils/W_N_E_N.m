function [ Wnen ] = W_N_E_N( VnE, VnN, L, H )
%由东向速度VnE，北向速度VnN，纬度L和高度H计算Wnen（ENU坐标系）

[RM,RN] = R_M_N(L);

Wnen = [-VnN/(RM + H);VnE/(RN + H);VnE/(RN + H)*tan(L)];
end
