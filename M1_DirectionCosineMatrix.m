%基于方向余弦矩阵的载体姿态解算程序（ENU坐标系）

function [ ] = M1_DirectionCosineMatrix( )
clc;clear;
close all;
addpath(genpath('Utils'));

%不变量
Weie = [0;0;7.292115e-5];%地球自转角速度
deltaT = 10e-3*2;%惯组数据采样周期

%初始化
dataPath = '../Data/';%数据文件存放地址
fidOut = fopen([dataPath,'INS.txt'],'r');%惯导解算结果（1.提供载体的初始状态；2.解算结果绘图时作为参考值）
initOut = fscanf(fidOut,'%e',[13,1]);
psi_0 = degree2radian(initOut(11));%偏航角ψ
theta_0 = degree2radian(initOut(12));%俯仰角θ
gamma_0 = degree2radian(initOut(13));%横滚角γ
lambda_0 = degree2radian(initOut(2));%经度λ
L_0 = degree2radian(initOut(3));%纬度L
num = 0;

%初始计算
Tnb_0 = T_N_B(psi_0,theta_0,gamma_0);%方向余弦矩阵
Cne_0 = C_N_E(lambda_0, L_0);%位置矩阵
Wnie_0 = Cne_0*Weie;
Wnen_0 = [0;0;0];%静基座数据

%读入初始Wbib_0，计算Wbnb_0
fidIn = fopen([dataPath,'IMU.txt'],'r');%惯组输出的待解算数据
imu_0 = fscanf(fidIn,'%e',[7,1]);
Wbib_0 = [imu_0(5);imu_0(6);imu_0(7)];%陀螺仪数据

Wnin_0 = Wnie_0 + Wnen_0;
Wbnb_0 = Wbib_0 - Tnb_0'*Wnin_0;

%姿态解算
Tnb = Tnb_0;
plotIndex = 1;
figure;
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
        %计算Tnb
        Tnb = R_K_4_DirectionCosineMatrix(deltaT,Tnb,Wbnb_0,Wbnb_1,Wbnb);
        
        %正交化
        Tnb = Orthogonalization_Tnb(Tnb);
        
        %计算psi, theta, gamma
        [psi,theta,gamma] = AnttitudeAngle_Tnb(Tnb);
        
        Wbnb_0 = Wbnb;
        num = 0;
        %画图
        if(mod(plotIndex*2,100) == 0)%每50个解算结果对应1个标准解算结果
            ins = fscanf(fidOut,'%e',[13,1]);
            psi_s = degree2radian(ins(11));
            theta_s = degree2radian(ins(12));
            gamma_s = degree2radian(ins(13));
            
            subplot(2,2,1)
            plot(plotIndex*deltaT,psi_s,'m*',plotIndex*deltaT,psi,'b.')
            title('偏航角解算结果')
            xlabel('t/s')
            ylabel('ψ/rad')
            legend('标准','结果','Location','northwest')
            hold on;
            subplot(2,2,2)
            plot(plotIndex*deltaT,theta_s,'m*',plotIndex*deltaT,theta,'b.')
            title('俯仰角解算结果')
            xlabel('t/s')
            ylabel('θ/rad')
            legend('标准','结果','Location','southeast')
            hold on;
            subplot(2,2,3)
            plot(plotIndex*deltaT,gamma_s,'m*',plotIndex*deltaT,gamma,'b.')
            title('横滚角解算结果')
            xlabel('t/s')
            ylabel('γ/rad')
            legend('标准','结果','Location','southwest')
            hold on;
            subplot(2,2,4)
            plot(plotIndex*deltaT,psi_s - psi,'m.',plotIndex*deltaT,theta_s - theta,'b.',plotIndex*deltaT,gamma_s - gamma,'c.')
            title('解算结果误差')
            xlabel('t/s')
            ylabel('ε/rad')
            legend('ψ','θ','γ','Location','northwest')
            hold on;
            drawnow;
        end
        plotIndex = plotIndex + 1;
    else
        %已读取数据0，数据1
        Wbnb_1 = Wbnb;
    end
end

fclose(fidIn);
fclose(fidOut);
disp('基于方向余弦矩阵的姿态解算结束！');
end
