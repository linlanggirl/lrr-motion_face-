function [Z,E] = solve_lrr(X,lambda)
Q = orth(X');   % orth函数是求矩阵正交基
A = X*Q;     % 
[Z,E] = lrra(X,A,lambda);
Z = Q*Z;