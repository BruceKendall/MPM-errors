function [tbl]= lionfish9Compare
%% Parameters in the paper
Ml=0.350; % Instantaneous mortality of larvae
Ma=0.052; % Instantaneous mortality of adults
Mj=0.165; % Instantaneous mortality of juveniles
rho=0.46; % Proportion female
Dl=30; % Duration of larval stage
Me=0.310; % Instantaneous mortality of eggs
f=194577; % Strictly speaking there is a mistake in this calculation, but we probably keep it as is.
De=3; % Duration of egg stage
%
Se=exp(-Me*De); % Survival over 3 days
Sl=exp(-Ml*Dl); % Monthly survival
Sj=exp(-Mj); % Monthly survival
Sa=exp(-Ma); % Monthly survival

%%  Original matrix with adult survival incorporated into fertility and correction of transition rates with inclusion of fertility for juveniles FAS 1
PA=Sa;
PJ=(10/11)*Sj;% This has been modified (duration in juvenile should be 11 months)
GL=Sl;
GJ=(1/11)*Sj;% This has been modified (duration in juvenile should be 11months)
RA=rho*f*Se*Sa^(27/30);
RJ=(1/11)*rho*f*Se*Sj^(27/30);
A1=[0,RJ,RA;GL,PJ,0;0,GJ,PA];
[lambda1, v1, w1] = eigen(A1);
S1=v1(:)*w1(:)'/(v1(:)'*w1(:)); %Sensitivity Matrix
sens1=[S1(2,1), S1(1,2)*A1(1,2)*27/30/Sj+10/11*S1(2,2)+1/11*S1(3,2), S1(3,3)+S1(1,3)*A1(1,3)*27/30/Sa, S1(1,2)*A1(1,2)/f+S1(1,3)*A1(1,3)/f];

%% Matching proportion (Fixed Stage Duration assuming Lambda = 1) SAS 1
PA=Sa;
PJ=sum(Sj.^[0:9])/sum(Sj.^[0:10])*Sj;% This has been modified to match the proportion transitioning
GL=Sl;
GJ=Sj.^10/sum(Sj.^[0:10])*Sj;% This has been modified to match the proportion transitioning
RA=rho*f*Se*Sa^(27/30);% This has been modified to include adult survival
RJ=Sj.^10/sum(Sj.^[0:10])*rho*f*Se*Sj^(27/30);
A2=[0,RJ,RA;GL,PJ,0;0,GJ,PA];
[lambda2, v2, w2] = eigen(A2);
S2=v2(:)*w2(:)'/(v2(:)'*w2(:)); %Sensitivity Matrix
sens2=[S2(2,1),
    S2(1,2)*((109*Se*Sj^(99/10)*f)/(10*(Sj^10 + Sj^9 + Sj^8 + Sj^7 + Sj^6 + Sj^5 + Sj^4 + Sj^3 + Sj^2 + Sj + 1)) - (Se*Sj^(109/10)*f*(10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1))/(Sj^10 + Sj^9 + Sj^8 + Sj^7 + Sj^6 + Sj^5 + Sj^4 + Sj^3 + Sj^2 + Sj + 1)^2)*rho+...
    S2(2,2)*(10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1)/(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10 + 10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1)+...
    S2(3,2)*(Sj^10*(Sj^10 + 2*Sj^9 + 3*Sj^8 + 4*Sj^7 + 5*Sj^6 + 6*Sj^5 + 7*Sj^4 + 8*Sj^3 + 9*Sj^2 + 10*Sj + 11))/(Sj^10 + Sj^9 + Sj^8 + Sj^7 + Sj^6 + Sj^5 + Sj^4 + Sj^3 + Sj^2 + Sj + 1)^2,
    S2(1,3)*((9*Se*f)/(10*Sa^(1/10)))*rho+ S2(3,3),
    S2(1,2)*A2(1,2)/f+S2(1,3)*A2(1,3)/f]';
%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda) AAS 1
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    A3=[0,                    rho*f*Se*Sj^(27/30)*(Sj/L).^10/sum((Sj/L).^[0:10]),             rho*f*Se*Sa^(27/30);
        Sl,        sum((Sj/L).^[0:9])/sum((Sj/L).^[0:10])*Sj,                                                0;
        0,                       (Sj/L)^10/sum((Sj/L).^[0:10])*Sj ,                                           Sa];
    L=max(real(eig(A3)));
end
[lambda3, v3, w3] = eigen(A3);
%% Grid over which lambda is estimated for numerical differentiation
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
%%
for PAR=1:4 % for the four parameters to which numerial differentiation is to be done
    for j=1:100
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            Sl2=P(j,1,PAR);
            Sj2=P(j,2,PAR);
            Sa2=P(j,3,PAR);
            f2=P(j,4,PAR);
            A=[0,                    rho*f2*Se*Sj2^(27/30)*(Sj2/L).^10/sum((Sj2/L).^[0:10]),    rho*f2*Se*Sa2^(27/30);
                Sl2,        sum((Sj2/L).^[0:9])/sum((Sj2/L).^[0:10])*Sj2,                                                 0;
                0,                       (Sj2/L)^10/sum((Sj2/L).^[0:10])*Sj2 ,                                              Sa2];
            L=max(real(eig(A)));
        end
        lambda(j)=L;
    end
    DP=P(:,PAR,PAR);
    lambda=lambda(:);
    D=(lambda(2:end)-lambda(1))./(DP(2:end)-DP(1)); % Differences calcualted with different dx values
    STATS=regstats(D,DP(2:end)-DP(1));
    sens3(PAR)=STATS.beta(1); % Intercept gives the derivative
end

%%  Original matrix with adult survival incorporated into fertility and correction of transition rates with inclusion of fertility for juveniles FAS 2
PA=Sa;
PJ=(10/11)*Sj;% This has been modified (duration in juvenile should be 11 months)
GL=Sl;
GJ=(1/11)*Sj;% This has been modified (duration in juvenile should be 11months)
RA=rho*f*Se*Sa^(27/30);
A4=[0,0,RA;GL,PJ,0;0,GJ,PA];
[lambda4, v4, w4] = eigen(A4);
S4=v4(:)*w4(:)'/(v4(:)'*w4(:)); %Sensitivity Matrix
sens4=[S4(2,1),
    S4(2,2)*10/11+S4(3,2)/11 ,
    S4(1,3)*((9*Se*f)/(10*Sa^(1/10)))*rho + S4(3,3),
    S4(1,3)*(Sa^(9/10)*Se)*rho ]';

%% Matching proportion (Fixed Stage Duration assuming Lambda = 1) SAS 2
PA=Sa;
PJ=sum(Sj.^[0:9])/sum(Sj.^[0:10])*Sj;% This has been modified to match the proportion transitioning
GL=Sl;
GJ=Sj.^10/sum(Sj.^[0:10])*Sj;% This has been modified to match the proportion transitioning
RA=rho*f*Se*Sa^(27/30);% This has been modified to include adult survival
A5=[0,0,RA;GL,PJ,0;0,GJ,PA];
[lambda5, v5, w5] = eigen(A5);
S5=v5(:)*w5(:)'/(v5(:)'*w5(:)); %Sensitivity Matrix
sens5=[S5(2,1),
    S5(2,2)*(10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1)/(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10 + 10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1)+...
    S5(3,2)*(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10)/(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10 + 10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1),
    S5(1,3)*((9*Se*f)/(10*Sa^(1/10)))*rho+ S5(3,3),
    S5(1,3)*(Sa^(9/10)*Se)*rho]';

%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda) AAS 2
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    A6=[0,                    0,                        rho*f*Se*Sa^(27/30);
        Sl,        sum((Sj/L).^[0:9])/sum((Sj/L).^[0:10])*Sj,             0;
        0,                       (Sj/L)^10/sum((Sj/L).^[0:10])*Sj ,        Sa];
    L=max(real(eig(A6)));
end
[lambda6, v6, w6] = eigen(A6);
for PAR=1:4 % for the four parameters to which numerial differentiation is to be done
    for j=1:100
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            Sl2=P(j,1,PAR);
            Sj2=P(j,2,PAR);
            Sa2=P(j,3,PAR);
            f=P(j,4,PAR);
            A=[0,                    0,                        rho*f*Se*Sa2^(27/30);
                Sl2,        sum((Sj2/L).^[0:9])/sum((Sj2/L).^[0:10])*Sj2,       0;
                0,                       (Sj2/L)^10/sum((Sj2/L).^[0:10])*Sj2 ,  Sa2];
            L=max(real(eig(A)));
        end
        lambda(j)=L;
    end
    DP=P(:,PAR,PAR);
    lambda=lambda(:);
    D=(lambda(2:end)-lambda(1))./(DP(2:end)-DP(1)); % Differences calcualted with different dx values
    STATS=regstats(D,DP(2:end)-DP(1));
    sens6(PAR)=STATS.beta(1); % Intercept gives the derivative
end

%%  Original matrix with adult survival incorporated into fertility and correction of transition rates with inclusion of fertility for juveniles FAS3
PA=Sa;
PJ=(10/11)*Sj;% This has been modified (duration in juvenile should be 11 months)
GL=Sl;
GJ=(1/11)*Sj;% This has been modified (duration in juvenile should be 11months)
RA=rho*f*Se;
A7=[0,0,RA;GL,PJ,0;0,GJ,PA];
[lambda7, v7, w7] = eigen(A7);
S7=v7(:)*w7(:)'/(v7(:)'*w7(:)); %Sensitivity Matrix
sens7=[S7(2,1),
    S7(2,2)*10/11 +    S7(3,2)/11   ,
    S7(3,3),
    S7(1,3)*Se*rho]';
%% Matching proportion (Fixed Stage Duration assuming Lambda = 1) SAS 3
PA=Sa;
PJ=sum(Sj.^[0:9])/sum(Sj.^[0:10])*Sj;% This has been modified to match the proportion transitioning
GL=Sl;
GJ=Sj.^10/sum(Sj.^[0:10])*Sj;% This has been modified to match the proportion transitioning
RA=rho*f*Se;% This has been modified to include adult survival
A8=[0,0,RA;GL,PJ,0;0,GJ,PA];
[lambda8, v8, w8] = eigen(A8);
S8=v8(:)*w8(:)'/(v8(:)'*w8(:)); %Sensitivity Matrix
sens8=[S8(2,1),
    S8(2,2)*(10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1)/(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10 + 10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1)+...
    S8(3,2)*(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10)/(Sj^20 + 2*Sj^19 + 3*Sj^18 + 4*Sj^17 + 5*Sj^16 + 6*Sj^15 + 7*Sj^14 + 8*Sj^13 + 9*Sj^12 + 10*Sj^11 + 11*Sj^10 + 10*Sj^9 + 9*Sj^8 + 8*Sj^7 + 7*Sj^6 + 6*Sj^5 + 5*Sj^4 + 4*Sj^3 + 3*Sj^2 + 2*Sj + 1) ,
    S8(3,3),
    S8(1,3)*Se*rho]';
%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda) AAS 3
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    A9=[0,                    0,                                          rho*f*Se;
        Sl,        sum((Sj/L).^[0:9])/sum((Sj/L).^[0:10])*Sj,          0;
        0,                       (Sj/L)^10/sum((Sj/L).^[0:10])*Sj ,      Sa];
    L=max(real(eig(A9)));
end
[lambda9, v9, w9] = eigen(A9);
for PAR=1:4 % for the four parameters to which numerial differentiation is to be done
    for j=1:100
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            Sl2=P(j,1,PAR);
            Sj2=P(j,2,PAR);
            Sa2=P(j,3,PAR);
            f=P(j,4,PAR);
            A=[0,                    0,                                         rho*f*Se;
                Sl2,        sum((Sj2/L).^[0:9])/sum((Sj2/L).^[0:10])*Sj2,        0;
                0,                       (Sj2/L)^10/sum((Sj2/L).^[0:10])*Sj2 ,   Sa2];
            L=max(real(eig(A)));
        end
        lambda(j)=L;
    end
    DP=P(:,PAR,PAR);
    lambda=lambda(:);
    D=(lambda(2:end)-lambda(1))./(DP(2:end)-DP(1)); % Differences calcualted with different dx values
    STATS=regstats(D,DP(2:end)-DP(1));
    sens9(PAR)=STATS.beta(1); % Intercept gives the derivative
end

PARA=[exp(-Ml*Dl),exp(-Mj),exp(-Ma),f];
elas1=PARA.*sens1/lambda1;
elas2=PARA.*sens2/lambda2;
elas3=PARA.*sens3/lambda3;
elas4=PARA.*sens4/lambda4;
elas5=PARA.*sens5/lambda5;
elas6=PARA.*sens6/lambda6;
elas7=PARA.*sens7/lambda7;
elas8=PARA.*sens8/lambda8;
elas9=PARA.*sens9/lambda9;

%% Damping Ratio
D1=sort(eig(A1),'descend');
D2=sort(eig(A2),'descend');
D3=sort(eig(A3),'descend');
D4=sort(eig(A4),'descend');
D5=sort(eig(A5),'descend');
D6=sort(eig(A6),'descend');
D7=sort(eig(A7),'descend');
D8=sort(eig(A8),'descend');
D9=sort(eig(A9),'descend');
DR(1)=D1(1)/abs(D1(2));
DR(2)=D2(1)/abs(D2(2));
DR(3)=D3(1)/abs(D3(2));
DR(4)=D4(1)/abs(D4(2));
DR(5)=D5(1)/abs(D5(2));
DR(6)=D6(1)/abs(D6(2));
DR(7)=D7(1)/abs(D7(2));
DR(8)=D8(1)/abs(D8(2));
DR(9)=D9(1)/abs(D9(2));

%% Generation Time
GT(1)=lambda1*(v1*w1)/(v1(1)*A1(1,:)*w1);
GT(2)=lambda2*(v2*w2)/(v2(1)*A2(1,:)*w2);
GT(3)=lambda3*(v3*w3)/(v3(1)*A3(1,:)*w3);
GT(4)=lambda4*(v4*w4)/(v4(1)*A4(1,:)*w4);
GT(5)=lambda5*(v5*w5)/(v5(1)*A5(1,:)*w5);
GT(6)=lambda6*(v6*w6)/(v6(1)*A8(1,:)*w6);
GT(7)=lambda7*(v7*w7)/(v7(1)*A7(1,:)*w7);
GT(8)=lambda8*(v8*w8)/(v8(1)*A8(1,:)*w8);
GT(9)=lambda9*(v9*w9)/(v9(1)*A9(1,:)*w9);

tbl=table([lambda1;lambda2;lambda3;lambda4;lambda5;lambda6;lambda7;lambda8;lambda9],[w1';w2';w3';w4';w5';w6';w7';w8';w9'],[v1;v2;v3;v4;v5;v6;v7;v8;v9],[sens1;sens2;sens3;sens4;sens5;sens6;sens7;sens8;sens9],[elas1;elas2;elas3;elas4;elas5;elas6;elas7;elas8;elas9],GT',DR','VariableNames',{'Lambda','Stable_Stage_Distribution','Reproductive_Values','Sensitivity','Elasticity','Generation_Time','Damping_Ratio'})
%writetable(tbl,'lionfish9models.csv')

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
