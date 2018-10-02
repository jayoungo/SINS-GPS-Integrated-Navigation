function [ lambda, L ] = Location_Cne( Cne )
%由位置矩阵Cne计算载体经度lambda和纬度L

C_31 = Cne(3,1);
C_32 = Cne(3,2);
C_33 = Cne(3,3);

lambda = atan(C_32/C_31);
if(C_31 < 0)
    if(lambda < 0)
        lambda = lambda + pi;
    else
        lambda = lambda - pi;
    end
end
L = atan(C_33/sqrt(C_31^2 + C_32^2));
end
