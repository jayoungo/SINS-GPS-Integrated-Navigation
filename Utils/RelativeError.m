function [ err ] = RelativeError( S_0, S )
%由旧值S_0和新值S计算相对误差（变化率）的绝对值err

err = abs(S - S_0)/S_0;
end
