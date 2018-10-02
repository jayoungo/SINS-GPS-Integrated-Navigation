function [ Q ] = Q_AnttitudeAngle( psi, theta, gamma )
%由偏航角ψ，俯仰角θ,横滚角γ计算四元数Q

psi = psi/2;
theta = theta/2;
gamma = gamma/2;

S_psi = sin(psi);
S_theta = sin(theta);
S_gamma = sin(gamma);
C_psi = cos(psi);
C_theta = cos(theta);
C_gamma = cos(gamma);

Q = [C_psi*C_theta*C_gamma - S_psi*S_theta*S_gamma;...
    C_psi*S_theta*C_gamma - S_psi*C_theta*S_gamma;...
    C_psi*C_theta*S_gamma + S_psi*S_theta*C_gamma;...
    C_psi*S_theta*S_gamma + S_psi*C_theta*C_gamma];
end
