%基于四元数的载体姿态解算程序（ENU坐标系）

function [ ] = M2_Quaternion( )
clc;clear;
close all;
addpath(genpath('Utils'));

%不变量
Weie = [0;0;7.292115e-5];
deltaT = 10e-3*2;

%初始化
dataPath = '../Data/';
fidOut = fopen([dataPath,'INS.txt'],'r');
initOut = fscanf(fidOut,'%e',[13,1]);
psi_0 = degree2radian(initOut(11));
theta_0 = degree2radian(initOut(12));
gamma_0 = degree2radian(initOut(13));
lambda_0 = degree2radian(initOut(2));
L_0 = degree2radian(initOut(3));
num = 0;

%初始计算
Q_0 = Q_AnttitudeAngle(psi_0,theta_0,gamma_0);
Tnb_0 = T_N_B(psi_0,theta_0,gamma_0);
%Tnb_0 = T_N_B_Quaternion(Q_0);
Cne_0 = C_N_E(lambda_0, L_0);
Wnie_0 = Cne_0*Weie;
Wnen_0 = [0;0;0];%静基座数据

%读入初始Wbib_0，计算Wbnb_0
fidIn = fopen([dataPath,'IMU.txt'],'r');
imu_0 = fscanf(fidIn,'%e',[7,1]);
Wbib_0 = [imu_0(5);imu_0(6);imu_0(7)];

Wnin_0 = Wnie_0 + Wnen_0;
Wbnb_0 = Wbib_0 - Tnb_0'*Wnin_0;

Q = Q_0;
Tnb = Tnb_0;
anttitudeAngle = [];
while(fseek(fidIn,2,0) == 0)
    %解算未结束
    fseek(fidIn,-2,0);
    
    %读IMU，计算Wbnb
    imu = fscanf(fidIn,'%e',[7,1]);
    Wbib = [imu(5);imu(6);imu(7)];
    
    Wnin = Wnie_0 + Wnen_0;
    Wbnb = Wbib - Tnb'*Wnin;
    num = num + 1;
    
    if(num ~= 1)
        %已读取数据0，数据1，数据2
        %计算Q
        Q = R_K_4_Quaternion(deltaT,Q,Wbnb_0,Wbnb_1,Wbnb);
        
        %归一化
        Q = Normalization_Q(Q);
        
        %计算Tnb
        Tnb = T_N_B_Quaternion(Q);
        
        %计算psi, theta, gamma
        [psi,theta,gamma] = AnttitudeAngle_Tnb(Tnb);
        anttitudeAngle = [anttitudeAngle;[psi theta gamma]]; %#ok<*AGROW>
        
        Wbnb_0 = Wbnb;
        num = 0;
    else
        %已读取数据0，数据1
        Wbnb_1 = Wbnb;
    end
end

%画图
figure;
ins = fscanf(fidOut,'%e');
%数据对齐
ins = reshape(ins,13,size(ins,1)/13)';
anttitudeAngle_s = degree2radian(ins(:,11:13));
anttitudeAngle = anttitudeAngle(50:50:29950,:);
index = (1:599).*deltaT;

subplot(2,2,1)
plot(index,anttitudeAngle_s(:,1),'m-',index,anttitudeAngle(:,1),'b-')
title('偏航角解算结果')
xlabel('t/s')
ylabel('ψ/rad')
legend('标准','结果','Location','northwest')
hold on;
subplot(2,2,2)
plot(index,anttitudeAngle_s(:,2),'m-',index,anttitudeAngle(:,2),'b-')
title('俯仰角解算结果')
xlabel('t/s')
ylabel('θ/rad')
legend('标准','结果','Location','southeast')
hold on;
subplot(2,2,3)
plot(index,anttitudeAngle_s(:,3),'m-',index,anttitudeAngle(:,3),'b-')
title('横滚角解算结果')
xlabel('t/s')
ylabel('γ/rad')
legend('标准','结果','Location','southwest')
hold on;
subplot(2,2,4)
plot(index,anttitudeAngle_s(:,1) - anttitudeAngle(:,1),'m-',index,anttitudeAngle_s(:,2) - anttitudeAngle(:,2),'b-',index,anttitudeAngle_s(:,3) - anttitudeAngle(:,3),'c-')
title('解算结果误差')
xlabel('t/s')
ylabel('ε/rad')
legend('ψ','θ','γ','Location','northwest');

fclose(fidIn);
fclose(fidOut);
disp('基于四元数的姿态解算结束！');
end
