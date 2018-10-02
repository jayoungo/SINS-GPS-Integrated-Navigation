# Inertial Navigation System (INS) and GPS Integrated Navigation
# 惯性导航和GPS组合导航


## 文件结构说明

#### 主程序功能
　　*M1_DirectionCosineMatrix.m*  
　　**基于方向余弦矩阵的载体姿态解算程序**  
　　*M2_Quaternion.m*  
　　**基于四元数的载体姿态解算程序**  
　　*M3_SINS.m*  
　　**捷联式惯性导航系统解算程序**  
　　*M4_InitAlign.m*  
　　**惯性导航系统的初始对准**  
　　*M5_1_SINS_GPS.m*  
　　**SINS/GPS组合导航**  
　　*M5_2_SINS_GPS.m*  
　　**SINS/GPS组合导航（效果更佳）**  

#### Utils
　　被主函数调用的工具类函数集合  

#### Example results
　　解算结果示例（仅供参考）  


## 测试软件版本

　　**MATLAB R2017b**


## 代码标记说明

#### 坐标系
　　*b*  
　　**载体**  
　　*e*  
　　**地球**  
　　*i*  
　　**惯性系**  
　　*n*  
　　**导航系**  

#### 变量通式
　　*Axyzw*  
　　**A^xy_zw**

#### 惯性导航基础
　　王新龙. 惯性导航基础[M]. 西北工业大学出版社, 2013.
