function sens=numsensA1
% Numerical calculation of sensitivity
%% Parameters given in the paper
s0=0.54;
s1=0.38;
s2=0.78;
s3=0.73;
s4=0.83;
fn=2.37;
d2n=7;
d3n=7;
P=zeros(100,5,5);
P(:,1,1)=linspace(s1,s1+0.0001,100)';
P(:,2,1)=ones(100,1)*s2;
P(:,3,1)=ones(100,1)*s3;
P(:,4,1)=ones(100,1)*s4;
P(:,5,1)=ones(100,1)*fn;
%
P(:,1,2)=ones(100,1)*s1;
P(:,2,2)=linspace(s2,s2+0.0001,100)';
P(:,3,2)=ones(100,1)*s3;
P(:,4,2)=ones(100,1)*s4;
P(:,5,2)=ones(100,1)*fn;
%
P(:,1,3)=ones(100,1)*s1;
P(:,2,3)=ones(100,1)*s2;
P(:,3,3)=linspace(s3,s3+0.0001,100)';
P(:,4,3)=ones(100,1)*s4;
P(:,5,3)=ones(100,1)*fn;
%
P(:,1,4)=ones(100,1)*s1;;
P(:,2,4)=ones(100,1)*s2;
P(:,3,4)=ones(100,1)*s3;
P(:,4,4)=linspace(s4,s4+0.0001,100)';
P(:,5,4)=ones(100,1)*fn;
%
P(:,1,5)=ones(100,1)*s1;;
P(:,2,5)=ones(100,1)*s2;
P(:,3,5)=ones(100,1)*s3;
P(:,4,5)=ones(100,1)*s4;
P(:,5,5)=linspace(fn,fn+0.01,100)';
for PAR = 1:5
    for j=1:100
        s1=P(j,1,PAR);
        s2=P(j,2,PAR);
        s3=P(j,3,PAR);
        s4=P(j,4,PAR);
        fn=P(j,5,PAR);
        G2=s2^d2n*(1-s2)/(1-s2^d2n);
        P2=(1-s2^(d2n-1))/(1-s2^d2n)*s2;
        P3=(1-s3^(d3n-1))/(1-s3^d3n)*s3;
        G3=s3^d3n*(1-s3)/(1-s3^d3n);  
        A=zeros(5,5);
        A(2,1)=s0;
        A(3,2)=s1;
        A(4,3)=G2;
        A(3,3)=P2;
        A(5,4)=G3;
        A(4,4)=P3;
        A(5,5)=s4;
        A(1,5)=fn;
        lambda(j)=max(real(eig(A)));
    end
    DP=P(:,PAR,PAR);
    lambda=lambda(:);
    D=(lambda(2:end)-lambda(1))./(DP(2:end)-DP(1)); % Differences calcualted with different dx values
    STATS=regstats(D,DP(2:end)-DP(1));
    sens(PAR)=STATS.beta(1); % Intercept gives the derivative
end

