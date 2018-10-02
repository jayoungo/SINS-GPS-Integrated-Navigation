function [ d_Vn ] = d_V_N( Tnb, Fb, Wnen, Cne, Vn, g )
%由方向余弦矩阵Tnb，比力Fb，角速度Wnen，位置矩阵Cne，速度Vn和重力加速度g计算载体加速度d_Vn

Weie = [0;0;7.292115e-5];

d_Vn = Tnb*Fb - (2.*OmegaMatrix(Cne*Weie) + OmegaMatrix(Wnen))*Vn - [0;0;g];
end
