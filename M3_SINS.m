%捷联式惯性导航系统解算程序（ENU坐标系）

function [ ] = M3_SINS( )
clc;clear;
close all;
addpath(genpath('Utils'));

%常量赋值
Weie = [0;0;7.292115e-5];
%Re = 6378137;f = 1/298.257;g0 = 9.7803;
deltaT = 10e-3;
k1 = 3.828;%高度通道二阶阻尼系统参数（Ref：惯性导航基础P123.）
k2 = 3.2804;

%INS读初值
dataPath = '../Data/';
fidOut = fopen([dataPath,'INS.txt'],'r');
initOut = fscanf(fidOut,'%e',[13,1]);
lambda_0 = degree2radian(initOut(2));
L_0 = degree2radian(initOut(3));
H_0 = initOut(4);
Vn_0 = [initOut(5);initOut(6);initOut(7)];
psi_0 = degree2radian(initOut(11));
theta_0 = degree2radian(initOut(12));
gamma_0 = degree2radian(initOut(13));
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
Wbib_0 = [imu_0(5);imu_0(6);imu_0(7)];
Wbnb_0 = Wbib_0 - Tnb_0'*Wnin_0;
Fb_0 = [imu_0(2);imu_0(3);imu_0(4)];
d_Vn_0 = d_V_N(Tnb_0,Fb_0,Wnen_0,Cne_0,Vn_0,g_0);
%HIGH读Hc_0
fidHc = fopen([dataPath,'HIGH.txt'],'r');
high_0 = fscanf(fidHc,'%e',[3,1]);
Hc_0 = high_0(3);

%H，VnE，VnN，VnU解算（d_Vn）：二级R-K方法
%lambda，L，psi，theta，gamma解算（d_Q,d_Cne）：四级R-K方法
L = L_0;
H = H_0;
H_d = H_0;
Tnb = Tnb_0;
Wnen = Wnen_0;
Cne = Cne_0;
Vn = Vn_0;
VnU_d = Vn_0(3);
Q = Q_0;
d_VnU_d = d_Vn_0(3);
prog = 1;
anttitude_res = [];
location_res = [];
high_res = [];
speed_res = [];
while(fseek(fidIn,2,0) == 0)
    %解算未结束
    fseek(fidIn,-2,0);
    
    g = G_H(H);
    
    %IMU读Wbib，Fb
    imu = fscanf(fidIn,'%e',[7,1]);
    Wbib = [imu(5);imu(6);imu(7)];
    Fb = [imu(2);imu(3);imu(4)];
    num = num + 1;
    
    d_Vn = d_V_N(Tnb,Fb,Wnen,Cne,Vn,g);
    Vn(1) = R_K_2(deltaT,Vn(1),d_Vn_0(1),d_Vn(1));
    Vn(2) = R_K_2(deltaT,Vn(2),d_Vn_0(2),d_Vn(2));
    Vn(3) = R_K_2(deltaT,Vn(3),d_Vn_0(3),d_Vn(3));
    H = R_K_2(deltaT,H,Vn_0(3),Vn(3));
    if(num == 4)
        %HIGH读Hc
        high = fscanf(fidHc,'%e',[3,1]);
        Hc = high(3);
        
        Vn(3) = R_K_2(deltaT*4,VnU_d,d_VnU_d - k2*(H_d - Hc_0),d_Vn(3) - k2*(H - Hc));
        H = R_K_2(deltaT*4,H_d,VnU_d - k1*(H_d - Hc_0),Vn(3) - k1*(H - Hc));
        num = 0;
        
        Hc_0 = Hc;
        VnU_d = Vn(3);
        H_d = H;
        d_VnU_d = d_Vn(3);
    end
    speed_res = [speed_res;Vn']; %#ok<*AGROW>
    high_res = [high_res;H];
    
    d_Vn_0 = d_Vn;
    Vn_0 = Vn;
    
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
        anttitude_res = [anttitude_res;[psi theta gamma]];
        
        %计算lambda，L
        [lambda,L] = Location_Cne(Cne);
        location_res = [location_res;[lambda L]];
        
        disp(prog);
        prog = prog + 1;
    else
        %已读取数据0，数据1
        Wnen_1 = Wnen;
        Wbnb_1 = Wbnb;
    end
end

%画图
figure;
ins = fscanf(fidOut,'%e');
%数据对齐
ins = reshape(ins,13,size(ins,1)/13)';
anttitude_base = degree2radian(ins(:,11:13));
location_base = [degree2radian(ins(:,2:3)) ins(:,4)];
speed_base = ins(:,5:7);
anttitude_res = anttitude_res(50:50:29950,:);
location_res = location_res(50:50:29950,:);
high_res = high_res(100:100:59900,:);
speed_res = speed_res(100:100:59900,:);
index = (1:599).*(deltaT*100);

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
title('姿态解算结果误差')
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
title('速度解算结果误差')
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
title('位置解算结果误差','FontSize',10)
xlabel('t/s')
legend({'λ','L','H'},'FontSize',8,'Location','northeast');

fclose(fidIn);
fclose(fidOut);
fclose(fidHc);
disp('捷联惯导数据解算结束！');
end
