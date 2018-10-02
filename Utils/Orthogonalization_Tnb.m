function [ Tnb_new ] = Orthogonalization_Tnb( Tnb )
%方向余弦矩阵Tnb正交化

Tbn = Tnb';
for i=1:3
    Tbn = (Tbn + inv(Tbn'))/2;
end
Tnb_new = Tbn';
end
