function [A B C,curve]=CP_AJD(T,r)
%% TENSORBOX, 2018

% fast CP decomposition using nonsymmetric joint diagonalization;
% The tensor T is assumed to be of order 3 and dimension m x m x p
% and rank r <= min(m,p)
%
[A0 B0 S, curve]=scot(T);
%
ss=diag(sum(abs(S),3));
[is iord]=sort(-ss);
[m m1 p]=size(T);
A=inv(A0);
A=A(:,iord(1:r));
B=inv(conj(B0));
B=B(:,iord(1:r));
C=zeros(p,r);
for k=1:r
    C(:,k)=squeeze(S(iord(k),iord(k),:));
end

function [A B S,curve]=scot(T)
%
% Approximate nonsymmetric joint diagonalization of a tensor
% of order 3 and size m x m x p.
%
[m m1 p]=size(T);
T2=T;
for n=1:p
    T2(:,:,n)=T(:,:,n)';
end
[Ai L]=eig(T(:,:,2)/T(:,:,1));
[Bi L]=eig((T(:,:,2)'/T(:,:,1)'));
A=inv(Ai); B=inv(Bi); %%% initial values
A=A./(sqrt(sum(A.*conj(A),2))*ones(1,m));
B=B./(sqrt(sum(B.*conj(B),2))*ones(1,m));
incr=10;
value=kolik(T,A,B);
curve=value;
iter=0;
while value>1e-8*curve(1) && iter<1000 && abs(incr)>1e-8*value
    iter=iter+1;
    valold=value;
    A=update(T,B);
    % kolik(T,A,B)
    B=update(T2,A);
    value=kolik(T,A,B);
    incr=-value+valold;
    curve=[curve value];
end
% curve
S=T;
for n=1:p
    S(:,:,n)=A*T(:,:,n)*B';
end
%
function A=update(T,B)
[m m1 p]=size(T);
S=zeros(m,m,m);
for n=1:m
    U=zeros(m,m);
    for i=1:p
        u=T(:,:,i)*B(n,:)';
        U=U+u*u';
    end
    S(:,:,n)=U;
end
V=sum(S,3);
for n=1:m
    [W Lam]=eig(S(:,:,n),V-S(:,:,n));
    [Lmin,imin]=max(diag(Lam));
    %    A(j,:)=A(j,:)/norm(A(j,:))
    A(n,:)=W(:,imin)';
    A(n,:)=A(n,:)/norm(A(n,:));
end

function val=kolik(T,A,B)
[m m1 p]=size(T);
mask=ones(m,m)-eye(m);
val=0; val2=0;
for n=1:p
    aux=abs(A*T(:,:,n)*B').^2;
    val=val+sum(sum(mask.*aux));
    val2=val2+sum(diag(aux));
end
val=val/val2;