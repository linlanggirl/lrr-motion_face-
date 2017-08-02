% function [] = lrr_face_seg()
clear all;
clc;
data = loadmatfile('yaleb10.mat');   % 1~10对象的照片
X = data.X;     % 2016*640    前1~64列属于1号对象，65~128列属于2号对象...
gnd = data.cids;   % 1*640   前1~64列等于1，65~128列等于2号...
K = max(gnd);     % k = 10
tic;   % tic和toc用来记录matlab命令执行的时间。  tic用来保存当前时间

lambda=0.18;
%run lrr
% Z = solve_lrr(X,0.18);
Q = orth(X');   % orth函数是求矩阵正交基
A = X*Q;     % 
[Z,E] = lrra(X,A,lambda);
Z = Q*Z;

%post processing
[U,S,V] = svd(Z,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;

% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V;

n = size(V,1);
M = zeros(K,K,20);
rand('state',123456789);
for i=1:size(M,3)
    inds = false(n,1);
    while sum(inds)<K
        j = ceil(rand()*n);
        inds(j) = true;
    end
    M(:,:,i) = V(inds,:);
end

idx = kmeans(V,K,'emptyaction','singleton','start',M,'display','off');
toc;     % 使用toc来记录程序完成时间
acc =  1 - missclassGroups(idx,gnd,K)/length(idx);
disp(['seg acc=' num2str(acc)]);