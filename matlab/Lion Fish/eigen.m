function [lambda, v, w] = eigen(M)
%[lambda, v, w] = eigen(M)
% M: Popualtion matrix
% lambda: population growth rate
% v: reproductive values
% w: stable stage distribution
[W,D]=eig(M);
D=real(diag(D));
lambda=max(D);
imax=find(D==lambda);
w=W(:,imax);
w=w/sum(w);
[V,D]=eig(M');
D=real(diag(D));
imax=find(D==max(D));
v=abs(V(:,imax))';