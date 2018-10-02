function [ M_new ] = Normalization( M, major )
%对矩阵M按行（major = 1）或列（major = 2）执行向量规范化

if(major == 1)
    M_new = M./repmat(sqrt(sum(M.^2,2)),[1,3]);
elseif(major == 2)
    M_new = M./repmat(sqrt(sum(M.^2)),[3,1]);
else
    error('Input Error! In function Normalization.');
end
end
