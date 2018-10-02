%SINS/GPS组合导航（ENU坐标系；卡尔曼滤波；反馈校正）
%视组合点处误差已完全补偿，较"M5_1_SINS_GPS"效果更佳

function [ ] = M5_2_SINS_GPS( )
clc;clear;
close all;
addpath(genpath('Utils'));

%常量赋值
Weie = [0;0;7.292115e-5];
g_unit = 9.7803;%重力加速度
deltaT = 10e-3;
deltaT_G = 100*deltaT;%GPS数据采样周期
k1 = 3.828;%高度通道二阶阻尼系统参数（Ref：惯性导航基础P123.）
k2 = 3.2804;

%导航解算初始化
% %INS读初值
% dataPath = '../Data/';
% fidOut = fopen([dataPath,'INS.txt'],'r');
% initOut = fscanf(fidOut,'%e',[13,1]);
% lambda_0 = degree2radian(initOut(2));
% L_0 = degree2radian(initOut(3));
% H_0 = initOut(4);
% Vn_0 = [initOut(5);initOut(6);initOut(7)];
% psi_0 = degree2radian(initOut(11));
% theta_0 = degree2radian(initOut(12));
% gamma_0 = degree2radian(initOut(13));
dataPath = '../Data/';
fidGPS = fopen([dataPath,'GPS.gps'],'r');%GPS数据文件
initGPS = fscanf(fidGPS,'%e',[7,1]);
lambda_0 = degree2radian(initGPS(5));
L_0 = degree2radian(initGPS(6));
H_0 = initGPS(7);
Vn_0 = [0;0;0];%静基座数据
%[psi_0,theta_0,gamma_0] = M4_InitAlign();
psi_0 = degree2radian(89.850004);
theta_0 = degree2radian(1.9113951);
gamma_0 = degree2radian(1.0572407);
num = 0;

%初始计算
Tnb_0 = T_N_B(psi_0,theta_0,gamma_0);
Q_0 = Q_AnttitudeAngle(psi_0,theta_0,gamma_0);
Cne_0 = C_N_E(lambda_0, L_0);

%计算Wnen_0
Wnen_0 = W_N_E_N(Vn_0(1),Vn_0(2),L_0,H_0);
Wnin_0 = Cne_0*Weie + Wnen_0;
g_0 = G_H(H_0);

%IMU读Wbib_0，Fb_0：计算Wbnb_0，d_Vn_0
fidIn = fopen([dataPath,'IMU.txt'],'r');
imu_0 = fscanf(fidIn,'%e',[7,1]);
Wbib_0 = [imu_0(5);imu_0(6);imu_0(7)];%陀螺仪数据
Wbnb_0 = Wbib_0 - Tnb_0'*Wnin_0;
Fb_0 = [imu_0(2);imu_0(3);imu_0(4)];%加速度计数据
d_Vn_0 = d_V_N(Tnb_0,Fb_0,Wnen_0,Cne_0,Vn_0,g_0);
%HIGH读Hc_0
fidHc = fopen([dataPath,'HIGH.txt'],'r');%高度数据文件
high_0 = fscanf(fidHc,'%e',[3,1]);
Hc_0 = high_0(3);

%Ref：惯性导航基础P263.
%[组合导航卡尔曼滤波]误差参数
V_delta = 1;
W_epsilon = degree2radian(0.5)/3600;
W_d = degree2radian(0.5)/3600;
%F_delta = 1e-5*g_unit;
F_d = 1e-5*g_unit;
PHI = degree2radian(1);
POS = degree2radian(1);
P_delta = 5;%GPS水平位置测量误差（2～5m）

%[组合导航卡尔曼滤波]数值初始化
%Xk = [0 0 0 0 0 0 0 0 0 0]';
Pk = diag([PHI PHI PHI V_delta V_delta POS POS W_epsilon W_epsilon W_epsilon].^2);
%Q_KF = diag([W_d W_d W_d F_d F_d 0 0 0 0 0].^2);
Q_KF = diag([W_d W_d W_d F_d F_d F_d].^2);
R = diag([P_delta P_delta].^2);

Epsilon_b = [0;0;0];

%H，VnE，VnN，VnU解算（d_Vn）：二级R-K方法
%lambda，L，psi，theta，gamma解算（d_Q,d_Cne）：四级R-K方法
lambda = lambda_0;
L = L_0;
H = H_0;
H_d = H_0;
psi = psi_0;
theta = theta_0;
gamma = gamma_0;
Tnb = Tnb_0;
Wnen = Wnen_0;
Cne = Cne_0;
Vn = Vn_0;
VnU_d = Vn_0(3);
Q = Q_0;
d_VnU_d = d_Vn_0(3);

[RM,RN] = R_M_N(L);

%系统矩阵
Fn = Tnb*Fb_0;
F = zeros(10,10);
F(1,5) = -1/(RM + H);
F(1,2) = Weie(3)*sin(L) + Vn(1)/(RN + H)*tan(L);
F(1,3) = -Weie(3)*cos(L) - Vn(1)/(RN + H);
F(2,4) = 1/(RN + H);
F(2,6) = -Weie(3)*sin(L);
F(2,1) = -Weie(3)*sin(L) - Vn(1)/(RN + H)*tan(L);
F(2,3) = -Vn(2)/(RM + H);
F(3,4) = tan(L)/(RN + H);
F(3,6) = Weie(3)*cos(L) + Vn(1)/(RN + H)/(cos(L))^2;
F(3,1) = Weie(3)*cos(L) + Vn(1)/(RN + H);
F(3,2) = Vn(2)/(RM + H);
F(4,3) = Fn(2);
F(4,2) = -Fn(3);
F(4,4) = Vn(2)/(RM + H)*tan(L) - Vn(3)/(RM + H);
F(4,5) = 2*Weie(3)*sin(L) + Vn(1)/(RN + H)*tan(L);
F(4,6) = 2*Weie(3)*cos(L)*Vn(2) + Vn(1)*Vn(2)/(RN + H)/(cos(L))^2 + 2*Weie(3)*sin(L)*Vn(3);
F(5,1) = Fn(3);
F(5,3) = -Fn(1);
F(5,4) = -2*Weie(3)*sin(L) - 2*Vn(1)/(RN + H)*tan(L);
F(5,5) = -Vn(3)/(RM + H);
F(5,6) = -2*Weie(3)*cos(L)*Vn(1) - Vn(1)^2/(RN + H)/(cos(L))^2;
F(6,5) = 1/(RM + H);
F(7,4) = 1/((RN + H)*cos(L));
F(7,6) = Vn(1)/(RN + H)/cos(L)*tan(L);

F(1:3,8:10) = Tnb;

%系统噪声矩阵
G = zeros(10,6);
G(1:3,1:3) = Tnb;
G(4:5,4:6) = Tnb(1:2,:);

%连续时间系统离散化
PHIk_k_1 = eye(10) + F.*deltaT_G + F^2.*(deltaT_G^2/2);%一步转移阵
GAMMAk_1 = deltaT_G.*(eye(10) + F.*(deltaT_G/2) + F^2.*(deltaT_G^2/6))*G;%系统噪声驱动阵

%量测阵
Hk = [zeros(2,5) diag([(RM + H) (RN + H)*cos(L)]) zeros(2,3)];

anttitude_res = [];
location_res = [];
high_res = [];
speed_res = [];
Xk_res = [];
while(fseek(fidIn,2,0) == 0)
    %解算未结束
    fseek(fidIn,-2,0);
    
    g = G_H(H);
    
    %IMU读Wbib，Fb
    imu = fscanf(fidIn,'%e',[7,1]);
    Wbib = [imu(5);imu(6);imu(7)] - Epsilon_b;
    Fb = [imu(2);imu(3);imu(4)];
    num = num + 1;
    
    d_Vn = d_V_N(Tnb,Fb,Wnen,Cne,Vn,g);
    Vn(1) = R_K_2(deltaT,Vn(1),d_Vn_0(1),d_Vn(1));
    Vn(2) = R_K_2(deltaT,Vn(2),d_Vn_0(2),d_Vn(2));
    Vn(3) = R_K_2(deltaT,Vn(3),d_Vn_0(3),d_Vn(3));
    H = R_K_2(deltaT,H,Vn_0(3),Vn(3));
    if(mod(num,4) == 0)
        %HIGH读Hc
        high = fscanf(fidHc,'%e',[3,1]);
        Hc = high(3);
        
        Vn(3) = R_K_2(deltaT*4,VnU_d,d_VnU_d - k2*(H_d - Hc_0),d_Vn(3) - k2*(H - Hc));
        H = R_K_2(deltaT*4,H_d,VnU_d - k1*(H_d - Hc_0),Vn(3) - k1*(H - Hc));
        
        Hc_0 = Hc;
        VnU_d = Vn(3);
        H_d = H;
        d_VnU_d = d_Vn(3);
    end
    
    d_Vn_0 = d_Vn;
    
    Wnen = W_N_E_N(Vn(1),Vn(2),L,H);
    Wnin = Cne*Weie + Wnen;
    Wbnb = Wbib - Tnb'*Wnin;
    
    if(mod(num,2) == 0)
        %已读取数据0，数据1，数据2
        %计算Q
        Q = R_K_4_Quaternion(deltaT*2,Q,Wbnb_0,Wbnb_1,Wbnb);
        Wbnb_0 = Wbnb;
        
        %归一化
        Q = Normalization_Q(Q);
        
        %计算Cne
        Cne = R_K_4_C_N_E(deltaT*2,Cne,Wnen_0,Wnen_1,Wnen);
        Wnen_0 = Wnen;
        
        %计算Tnb
        Tnb = T_N_B_Quaternion(Q);
        
        %计算psi，theta，gamma
        [psi,theta,gamma] = AnttitudeAngle_Tnb(Tnb);
        
        %计算lambda，L
        [lambda,L] = Location_Cne(Cne);
    else
        %已读取数据0，数据1
        Wnen_1 = Wnen;
        Wbnb_1 = Wbnb;
    end
    
    %SINS/GPS松组合卡尔曼滤波（反馈校正）
    if(mod(num,100) == 0)
        %数据更新
        [RM,RN] = R_M_N(L);
        
        %系统矩阵
        Fn = Tnb*Fb;
        F = zeros(10,10);
        F(1,5) = -1/(RM + H);
        F(1,2) = Weie(3)*sin(L) + Vn(1)/(RN + H)*tan(L);
        F(1,3) = -Weie(3)*cos(L) - Vn(1)/(RN + H);
        F(2,4) = 1/(RN + H);
        F(2,6) = -Weie(3)*sin(L);
        F(2,1) = -Weie(3)*sin(L) - Vn(1)/(RN + H)*tan(L);
        F(2,3) = -Vn(2)/(RM + H);
        F(3,4) = tan(L)/(RN + H);
        F(3,6) = Weie(3)*cos(L) + Vn(1)/(RN + H)/(cos(L))^2;
        F(3,1) = Weie(3)*cos(L) + Vn(1)/(RN + H);
        F(3,2) = Vn(2)/(RM + H);
        F(4,3) = Fn(2);
        F(4,2) = -Fn(3);
        F(4,4) = Vn(2)/(RM + H)*tan(L) - Vn(3)/(RM + H);
        F(4,5) = 2*Weie(3)*sin(L) + Vn(1)/(RN + H)*tan(L);
        F(4,6) = 2*Weie(3)*cos(L)*Vn(2) + Vn(1)*Vn(2)/(RN + H)/(cos(L))^2 + 2*Weie(3)*sin(L)*Vn(3);
        F(5,1) = Fn(3);
        F(5,3) = -Fn(1);
        F(5,4) = -2*Weie(3)*sin(L) - 2*Vn(1)/(RN + H)*tan(L);
        F(5,5) = -Vn(3)/(RM + H);
        F(5,6) = -2*Weie(3)*cos(L)*Vn(1) - Vn(1)^2/(RN + H)/(cos(L))^2;
        F(6,5) = 1/(RM + H);
        F(7,4) = 1/((RN + H)*cos(L));
        F(7,6) = Vn(1)/(RN + H)/cos(L)*tan(L);
        
        F(1:3,8:10) = Tnb;
        
        %系统噪声矩阵
        G = zeros(10,6);
        G(1:3,1:3) = Tnb;
        G(4:5,4:6) = Tnb(1:2,:);
        
        %连续时间系统离散化
        PHIk_k_1_N = eye(10) + F.*deltaT_G + F^2.*(deltaT_G^2/2);%一步转移阵
        GAMMAk_1_N = deltaT_G.*(eye(10) + F.*(deltaT_G/2) + F^2.*(deltaT_G^2/6))*G;%系统噪声驱动阵
        
        %量测阵
        Hk_N = [zeros(2,5) diag([(RM + H) (RN + H)*cos(L)]) zeros(2,3)];
        
        %GPS数据读取
        GPS = fscanf(fidGPS,'%e',[7,1]);
        lambda_G = degree2radian(GPS(5));
        L_G = degree2radian(GPS(6));
        
        %量测量
        Zk = [(L - L_G)*(RM + H);(lambda - lambda_G)*(RN + H)*cos(L)];
        
        %离散型卡尔曼滤波基本方程
        Pk_k_1 = PHIk_k_1*Pk*PHIk_k_1' + GAMMAk_1*Q_KF*GAMMAk_1';
        Kk = Pk_k_1*Hk'/(Hk*Pk_k_1*Hk' + R);
        %Pk = inv(inv(Pk_k_1) + Hk'/R*Hk);
        Pk = (eye(10) - Kk*Hk)*Pk_k_1*(eye(10) - Kk*Hk)' + Kk*R*Kk';
        %Xk_k_1 = PHIk_k_1*Xk;
        %Xk = Xk_k_1 + Kk*(Zk - Hk*Xk_k_1);
        Xk = Kk*Zk;%Xk_k_1 = 0;
        
        PHIk_k_1 = PHIk_k_1_N;
        GAMMAk_1 = GAMMAk_1_N;
        Hk = Hk_N;
        
        Xk_res = [Xk_res Xk];
        
        %惯导系统反馈修正：数学平台Tnb、速度Vn、位置L&lambda
        phi = Xk(1:3,:);
        Tnb = (eye(3,3) + OmegaMatrix(phi))*Tnb;
        Vn = Vn - [Xk(4:5);0];%Ref：惯性导航基础P150,152.
        L = L - Xk(6);
        lambda = lambda - Xk(7);
        
        [psi,theta,gamma] = AnttitudeAngle_Tnb(Tnb);
        Q = Q_AnttitudeAngle(psi,theta,gamma);
        Cne = C_N_E(lambda, L);
        
        Epsilon_b = Xk(8:10);
        
        disp(num);
    end
    
    Vn_0 = Vn;
    
    %解算结果存储
    anttitude_res = [anttitude_res;[psi theta gamma]];
    location_res = [location_res;[lambda L]];
    high_res = [high_res;H];
    speed_res = [speed_res;Vn']; %#ok<*AGROW>
end

%画图
fidOut = fopen([dataPath,'INS.txt'],'r');
ins = fscanf(fidOut,'%e');
%数据对齐
ins(1:13,:) = [];
ins = reshape(ins,13,size(ins,1)/13)';
anttitude_base = degree2radian(ins(:,11:13));
location_base = [degree2radian(ins(:,2:3)) ins(:,4)];
speed_base = ins(:,5:7);
anttitude_res = anttitude_res(100:100:59900,:);
location_res = location_res(100:100:59900,:);
high_res = high_res(100:100:59900,:);
speed_res = speed_res(100:100:59900,:);

index = (1:599).*(deltaT*100);

figure;
subplot(3,4,1)
plot(index,anttitude_base(:,1),'m-',index,anttitude_res(:,1),'b-')
title('偏航角解算结果')
xlabel('t/s')
ylabel('ψ/rad')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,2)
plot(index,anttitude_base(:,2),'m-',index,anttitude_res(:,2),'b-')
title('俯仰角解算结果')
xlabel('t/s')
ylabel('θ/rad')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,3)
plot(index,anttitude_base(:,3),'m-',index,anttitude_res(:,3),'b-')
title('横滚角解算结果')
xlabel('t/s')
ylabel('γ/rad')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,4)
plot(index,anttitude_base(:,1) - anttitude_res(:,1),'m-',index,anttitude_base(:,2) - anttitude_res(:,2),'b-',index,anttitude_base(:,3) - anttitude_res(:,3),'c-')
title('姿态解算误差')
xlabel('t/s')
ylabel('ε/rad')
legend({'ψ','θ','γ'},'FontSize',8,'Location','best');

subplot(3,4,5)
plot(index,speed_base(:,1),'m-',index,speed_res(:,1),'b-')
title('东向速度解算结果')
xlabel('t/s')
ylabel('VE/(m/s)')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,6)
plot(index,speed_base(:,2),'m-',index,speed_res(:,2),'b-')
title('北向速度解算结果')
xlabel('t/s')
ylabel('VN/(m/s)')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,7)
plot(index,speed_base(:,3),'m-',index,speed_res(:,3),'b-')
title('天向速度解算结果')
xlabel('t/s')
ylabel('VU/(m/s)')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,8)
plot(index,speed_base(:,1) - speed_res(:,1),'m-',index,speed_base(:,2) - speed_res(:,2),'b-',index,speed_base(:,3) - speed_res(:,3),'c-')
title('速度解算误差')
xlabel('t/s')
ylabel('ε/(m/s)')
legend({'VE','VN','VU'},'FontSize',8,'Location','best');

subplot(3,4,9)
plot(index,location_base(:,1),'m-',index,location_res(:,1),'b-')
title('经度解算结果')
xlabel('t/s')
ylabel('λ/rad')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,10)
plot(index,location_base(:,2),'m-',index,location_res(:,2),'b-')
title('纬度解算结果')
xlabel('t/s')
ylabel('L/rad')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,11)
plot(index,location_base(:,3),'m-',index,high_res,'b-')
title('高度解算结果')
xlabel('t/s')
ylabel('H/m')
legend({'标准','结果'},'FontSize',8,'Location','best')
hold on;
subplot(3,4,12)
yyaxis left
plot(index,location_base(:,1) - location_res(:,1),index,location_base(:,2) - location_res(:,2))
ylabel('ε/rad')
yyaxis right
plot(index,location_base(:,3) - high_res)
ylabel('ε/m')
title('位置解算误差')
xlabel('t/s')
legend({'λ','L','H'},'FontSize',8,'Location','best');

fclose(fidIn);
fclose(fidOut);
fclose(fidHc);
fclose(fidGPS);
disp('INS/GPS组合导航解算结束！');
end
