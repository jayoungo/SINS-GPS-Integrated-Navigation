function [ psi, theta, gamma ] = AnttitudeAngle_Tnb( Tnb )
%由方向余弦矩阵Tnb计算载体姿态角ψ,θ,γ

T_12 = Tnb(1,2);
T_22 = Tnb(2,2);
T_31 = Tnb(3,1);
T_32 = Tnb(3,2);
T_33 = Tnb(3,3);
psi = atan((-T_12)/T_22);
if(T_22 < 0)
    psi = psi + pi;
elseif(psi < 0)
    psi = psi + 2*pi;
end
theta = atan(T_32/sqrt(T_12^2 + T_22^2));
gamma = atan((-T_31)/T_33);
if(T_33 < 0)
    if(gamma > 0)
        gamma = gamma - pi;
    else
        gamma = gamma + pi;
    end
end
end
