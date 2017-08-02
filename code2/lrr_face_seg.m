% function [] = lrr_face_seg()
clear all;
clc;
data = loadmatfile('yaleb10.mat');   % 1~10�������Ƭ
X = data.X;     % 2016*640    ǰ1~64������1�Ŷ���65~128������2�Ŷ���...
gnd = data.cids;   % 1*640   ǰ1~64�е���1��65~128�е���2��...
K = max(gnd);     % k = 10
tic;   % tic��toc������¼matlab����ִ�е�ʱ�䡣  tic�������浱ǰʱ��

lambda=0.18;
%run lrr
% Z = solve_lrr(X,0.18);
Q = orth(X');   % orth�����������������
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
toc;     % ʹ��toc����¼�������ʱ��
acc =  1 - missclassGroups(idx,gnd,K)/length(idx);
disp(['seg acc=' num2str(acc)]);