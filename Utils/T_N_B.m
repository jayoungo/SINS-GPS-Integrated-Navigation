function [ Tnb ] = T_N_B( psi, theta, gamma )
%ÓÉÆ«º½½Ç¦×£¬¸©Ñö½Ç¦ÈºÍºá¹ö½Ç¦Ã¼ÆËã·½ÏòÓàÏÒ¾ØÕóTnb

S_psi = sin(psi);
S_theta = sin(theta);
S_gamma = sin(gamma);
C_psi = cos(psi);
C_theta = cos(theta);
C_gamma = cos(gamma);

Tnb = [C_gamma*C_psi - S_gamma*S_theta*S_psi C_gamma*S_psi + S_gamma*S_theta*C_psi -S_gamma*C_theta;...
    -C_theta*S_psi C_theta*C_psi S_theta;...
    S_gamma*C_psi + C_gamma*S_theta*S_psi S_gamma*S_psi - C_gamma*S_theta*C_psi C_gamma*C_theta]';
end
