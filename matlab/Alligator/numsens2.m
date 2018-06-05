function sens=numsens2
%% Numerical Differentiation to Calculate Sensitivity
%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda)
% Paramters in the paper
Ml=0.350;
Ma=0.052;
Mj=0.165;
rho=0.46;
Dl=30;
Me=0.310;
f=194577;
De=3;
% The finite Survival Rates.
Se=exp(-Me*De);
Sl=exp(-Ml*Dl);
Sj=exp(-Mj);
Sa=exp(-Ma);
% Parameters at which lambda is evaluated
P=zeros(100,4,4);
P(:,1,1)=linspace(Sl,Sl+0.0000001,100)';
P(:,2,1)=ones(100,1)*Sj;
P(:,3,1)=ones(100,1)*Sa;
P(:,4,1)=ones(100,1)*f;
P(:,1,2)=ones(100,1)*Sl;
P(:,2,2)=linspace(Sj,Sj+0.0001,100)';
P(:,3,2)=ones(100,1)*Sa;
P(:,4,2)=ones(100,1)*f;
P(:,1,3)=ones(100,1)*Sl;
P(:,2,3)=ones(100,1)*Sj;
P(:,3,3)=linspace(Sa,Sa+0.0001,100)';
P(:,4,3)=ones(100,1,1)*f;
P(:,1,4)=ones(100,1,1)*Sl;
P(:,2,4)=ones(100,1,1)*Sj;
P(:,3,4)=ones(100,1,1)*Sa;
P(:,4,4)=linspace(f,f+0.01,100)';
for PAR=1:4
    for j=1:100
        Sl=P(j,1,PAR);
        Sj=P(j,2,PAR);
        Sa=P(j,3,PAR);
        f=P(j,4,PAR);
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            A=[0,rho*f*Se*Sa^(27/30)*(Sj/L).^10/sum((Sj/L).^[0:10])*Sj, rho*f*Se*Sa^(27/30);Sl,sum((Sj/L).^[0:9])/sum((Sj/L).^[0:10])*Sj,0;0,(Sj/L).^10/sum((Sj/L).^[0:10])*Sj,Sa];
            L=max(real(eig(A)));
        end
        lambda(j)=L;
    end
    DP=P(:,PAR,PAR);
    lambda=lambda(:);
    D=(lambda(2:end)-lambda(1))./(DP(2:end)-DP(1)); % Differences calcualted with different dx values 
    STATS=regstats(D,DP(2:end)-DP(1));
    sens(PAR)=STATS.beta(1); % Intercept gives the derivative
end