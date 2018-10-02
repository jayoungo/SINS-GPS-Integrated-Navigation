function [ M_new ] = Orthogonalization_Schmidt( M )
%使用施密特正交化（Schmidt orthogonalization）方法将3阶方阵M变换为正交矩阵

if(sum(size(M) == [3 3],2) ~= 2)
    error('Input Error! In function Orthogonalization_Schmidt.');
end

M(:,2) = M(:,2) - ((M(:,2)'*M(:,1))/(M(:,1)'*M(:,1))).*M(:,1);
M(:,3) = M(:,3) - ((M(:,3)'*M(:,1))/(M(:,1)'*M(:,1))).*M(:,1) - ((M(:,3)'*M(:,2))/(M(:,2)'*M(:,2))).*M(:,2);

M_new = Normalization(M,2);
end
