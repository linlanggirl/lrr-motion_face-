function [Z,E] = solve_lrr(X,lambda)
Q = orth(X');   % orth�����������������
A = X*Q;     % 
[Z,E] = lrra(X,A,lambda);
Z = Q*Z;