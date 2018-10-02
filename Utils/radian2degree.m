function [ d ] = radian2degree( r )
%弧度值（矩阵）r转换为角度值d

d = r.*(180/pi);
end
