function [E] = solve_l1l2(W,lambda)
n = size(W,2); %返回W的第二个维度 列
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);  %返回矩阵w的二范数，如果是向量的话返回的是向量的欧式范数
if nw>lambda
    x = (nw-lambda)*w/nw;  
else
    x = zeros(length(w),1);
end
end