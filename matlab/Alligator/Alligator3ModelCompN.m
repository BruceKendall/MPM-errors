function [tbl]=Alligator3ModelCompN
%% Parameters in the paper for the North population
se=0.54;
sl=0.38;
sj=0.78;
ss=0.73;
sa=0.83;
f=2.37;
d1=1;
d2=7;
d3=7;

%% Correction of fertility and reduction of stage to four stage for north population FAS1
G1=sl;
P2=(d2-1)/d2*sj;
G2=1/d2*sj;
P3=(d3-1)/d3*ss;
G3=1/d3*ss;
P4=sa;
F3=ss^(9/12)*se*f*1/d3;
F4=sa^(9/12)*se*f;
A1=[0,0,F3,F4;
    G1,P2,0,0;
    0, G2,P3,0;
    0,0,G3,P4];
[lambda1, v1, w1] = eigen(A1);
S1=v1(:)*w1(:)'/(v1(:)'*w1(:)); %Sensitivity Matrix
sens1=[S1(2,1), S1(2,2)*(d2 - 1)/d2+S1(3,2)*1/d2, S1(3,3)*(d3 - 1)/d3+S1(4,3)*1/d3+S1(1,3)*((3*f)/(4*d3*ss^(1/4)))*se, S1(4,4)+S1(1,4)*(3*f*se)/(4*sa^(1/4)),   S1(1,3)*(ss^(3/4)/d3)*se   + S1(1,4)*sa^(3/4)*se];

%% Correction of fertility and reduction of stage to four stage for north population SAS 1
    G1=sl;    
    P2=sum((sj).^[0:(d2-2)])/sum((sj).^[0:(d2-1)])*sj;
    G2=(sj).^(d2-1)/sum((sj).^[0:(d2-1)])*sj;%
    P3=sum((ss).^[0:(d3-2)])/sum((ss).^[0:(d3-1)])*ss;
    G3=(ss).^(d3-1)/sum((ss).^[0:(d3-1)])*ss;%
    P4=sa;
    F3=ss^(9/12)*se*f*(ss).^(d3-1)/sum((ss).^[0:(d3-1)]);
    F4=sa^(9/12)*se*f;
A2=[0,0,F3,F4;
    G1,P2,0,0;
    0, G2,P3,0;
    0,0,G3,P4];
[lambda2, v2, w2] = eigen(A2);
S2=v2(:)*w2(:)'/(v2(:)'*w2(:)); %Sensitivity Matrix
sens2 (1)=S2(2,1) ;
sens2 (2)=   S2(2,2)*(6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1)/(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6 + 6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1) + ...
    S2(3,2)*(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6)/(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6 + 6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1) ;
sens2 (3)=   S2(3,3)*(6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1)/(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6 + 6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1) + ...
    S2(4,3)*(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6)/(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6 + 6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1) + ...
    S2(1,3)*((27*f*ss^(23/4))/(4*(ss^6 + ss^5 + ss^4 + ss^3 + ss^2 + ss + 1)) - (f*ss^(27/4)*(6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1))/(ss^6 + ss^5 + ss^4 + ss^3 + ss^2 + ss + 1)^2)*se ;
sens2 (4)=    S2(4,4) +...
    S2(1,4)*(3*f*se)/(4*sa^(1/4));
sens2 (5)=    S2(1,3)*(ss^(27/4)/(ss^6 + ss^5 + ss^4 + ss^3 + ss^2 + ss + 1))*se  +...
    S2(1,4)*sa^(3/4)*se;

%% North population fixed duration with labmda AAS 1
A3=zeros(4,4);
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    G1=sl;
    P2=sum((sj/L).^[0:(d2-2)])/sum((sj/L).^[0:(d2-1)])*sj;
    G2=(sj/L).^(d2-1)/sum((sj/L).^[0:(d2-1)])*sj;%
    P3=sum((ss/L).^[0:(d3-2)])/sum((ss/L).^[0:(d3-1)])*ss;
    G3=(ss/L).^(d3-1)/sum((ss/L).^[0:(d3-1)])*ss;%
    P4=sa;
    F3=ss^(9/12)*se*f*(ss/L).^(d3-1)/sum((ss/L).^[0:(d3-1)]);
    F4=sa^(9/12)*se*f;
    A3=[0,0,F3,F4;
        G1,P2,0,0;
        0, G2,P3,0;
        0,0,G3,P4];
    L=max(real(eig(A3)));
end
[lambda3, v3, w3] = eigen(A3);
%% Grid over which lambda is estimated for numerical differentiation
P=zeros(100,5,5);
P(:,1,1)=linspace(sl,sl+0.0001,100)';
P(:,2,1)=ones(100,1)*sj;
P(:,3,1)=ones(100,1)*ss;
P(:,4,1)=ones(100,1)*sa;
P(:,5,1)=ones(100,1)*f;
%
P(:,1,2)=ones(100,1)*sl;
P(:,2,2)=linspace(sj,sj+0.0001,100)';
P(:,3,2)=ones(100,1)*ss;
P(:,4,2)=ones(100,1)*sa;
P(:,5,2)=ones(100,1)*f;
%
P(:,1,3)=ones(100,1)*sl;
P(:,2,3)=ones(100,1)*sj;
P(:,3,3)=linspace(ss,ss+0.0001,100)';
P(:,4,3)=ones(100,1)*sa;
P(:,5,3)=ones(100,1)*f;
%
P(:,1,4)=ones(100,1)*sl;
P(:,2,4)=ones(100,1)*sj;
P(:,3,4)=ones(100,1)*ss;
P(:,4,4)=linspace(sa,sa+0.0001,100)';
P(:,5,4)=ones(100,1)*f;
%
P(:,1,5)=ones(100,1)*sl;
P(:,2,5)=ones(100,1)*sj;
P(:,3,5)=ones(100,1)*ss;
P(:,4,5)=ones(100,1)*sa;
P(:,5,5)=linspace(f,f+0.01,100)';
%
for PAR = 1:5
    for j=1:100
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            sl2=P(j,1,PAR);
            sj2=P(j,2,PAR);
            ss2=P(j,3,PAR);
            sa2=P(j,4,PAR);
            f2=P(j,5,PAR);
            G1=sl2;
            P2=sum((sj2/L).^[0:(d2-2)])/sum((sj2/L).^[0:(d2-1)])*sj2;
            G2=(sj2/L).^(d2-1)/sum((sj2/L).^[0:(d2-1)])*sj2;%
            P3=sum((ss2/L).^[0:(d3-2)])/sum((ss2/L).^[0:(d3-1)])*ss2;
            G3=(ss2/L).^(d3-1)/sum((ss2/L).^[0:(d3-1)])*ss2;%
            P4=sa2;
            F3=ss2^(9/12)*se*f2*(ss2/L).^(d3-1)/sum((ss2/L).^[0:(d3-1)]);
            F4=sa2^(9/12)*se*f2;
            A=[0,0,F3,F4;
                G1,P2,0,0;
                0, G2,P3,0;
                0,0,G3,P4];
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

%% Correction of fertility and reduction of stage to four stage for north population FAS2
G1=sl;
P2=(d2-1)/d2*sj;
G2=1/d2*sj;
P3=(d3-1)/d3*ss;
G3=1/d3*ss;
P4=sa;
F4=sa^(9/12)*se*f;
A4=[0,0,0,F4;
    G1,P2,0,0;
    0, G2,P3,0;
    0,0,G3,P4];
[lambda4, v4, w4] = eigen(A4);
S4=v4(:)*w4(:)'/(v4(:)'*w4(:)); %Sensitivity Matrix
sens4=[S4(2,1), S4(2,2)*(d2 - 1)/d2+S4(3,2)*1/d2, S4(3,3)*(d3 - 1)/d3+S4(4,3)*1/d3, S4(4,4)+S4(1,4)*(3*f*se)/(4*sa^(1/4)),  S4(1,4)*sa^(3/4)*se];

%% Correction of fertility and reduction of stage to four stage for north population SAS 2
G1=sl;    
    P2=sum((sj).^[0:(d2-2)])/sum((sj).^[0:(d2-1)])*sj;
    G2=(sj).^(d2-1)/sum((sj).^[0:(d2-1)])*sj;%
    P3=sum((ss).^[0:(d3-2)])/sum((ss).^[0:(d3-1)])*ss;
    G3=(ss).^(d3-1)/sum((ss).^[0:(d3-1)])*ss;%
    P4=sa;
    F4=sa^(9/12)*se*f;
    A5=[0,0,0,F4;
        G1,P2,0,0;
        0, G2,P3,0;
        0,0,G3,P4];
[lambda5, v5, w5] = eigen(A5);
S5=v5(:)*w5(:)'/(v5(:)'*w5(:)); %Sensitivity Matrix
sens5(1)=S5(2,1);
sens5(2)=    S5(2,2)*(6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1)/(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6 + 6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1) +...
    S5(3,2)*(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6)/(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6 + 6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1);
sens5(3)=    S5(3,3)*(6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1)/(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6 + 6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1)+...
    S5(4,3)*(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6)/(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6 + 6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1);
sens5(4)=    S5(4,4)+...
    S5(1,4)*(3*f*se)/(4*sa^(1/4));
sens5(5)=   S5(1,4)*sa^(3/4)*se;
%% North population fixed duration with labmda AAS 2
A6=zeros(4,4);
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    G1=sl;
    P2=sum((sj/L).^[0:(d2-2)])/sum((sj/L).^[0:(d2-1)])*sj;
    G2=(sj/L).^(d2-1)/sum((sj/L).^[0:(d2-1)])*sj;%
    P3=sum((ss/L).^[0:(d3-2)])/sum((ss/L).^[0:(d3-1)])*ss;
    G3=(ss/L).^(d3-1)/sum((ss/L).^[0:(d3-1)])*ss;%
    P4=sa;
    F4=sa^(9/12)*se*f;
    A6=[0,0,0,F4;
        G1,P2,0,0;
        0, G2,P3,0;
        0,0,G3,P4];
    L=max(real(eig(A6)));
end
[lambda6, v6, w6] = eigen(A6);
for PAR = 1:5
    for j=1:100
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            sl2=P(j,1,PAR);
            sj2=P(j,2,PAR);
            ss2=P(j,3,PAR);
            sa2=P(j,4,PAR);
            f2=P(j,5,PAR);
            G1=sl2;
            P2=sum((sj2/L).^[0:(d2-2)])/sum((sj2/L).^[0:(d2-1)])*sj2;
            G2=(sj2/L).^(d2-1)/sum((sj2/L).^[0:(d2-1)])*sj2;%
            P3=sum((ss2/L).^[0:(d3-2)])/sum((ss2/L).^[0:(d3-1)])*ss2;
            G3=(ss2/L).^(d3-1)/sum((ss2/L).^[0:(d3-1)])*ss2;%
            P4=sa2;
            F4=sa2^(9/12)*se*f2;
            A=[0,0,0,F4;
                G1,P2,0,0;
                0, G2,P3,0;
                0,0,G3,P4];
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

%% Correction of fertility and reduction of stage to four stage for north population FAS 3
G1=sl;
P2=(d2-1)/d2*sj;
G2=1/d2*sj;
P3=(d3-1)/d3*ss;
G3=1/d3*ss;
P4=sa;
F4=se*f;
A7=[0,0,0,F4;
    G1,P2,0,0;
    0, G2,P3,0;
    0,0,G3,P4];
[lambda7, v7, w7] = eigen(A7);
S7=v7(:)*w7(:)'/(v7(:)'*w7(:)); %Sensitivity Matrix
sens7=[S7(2,1), S7(2,2)*(d2 - 1)/d2+S7(3,2)*1/d2, S7(3,3)*(d3 - 1)/d3+S7(4,3)*1/d3, S7(4,4),   S7(1,4)*se];

%% Correction of fertility and reduction of stage to four stage for north population SAS 3
G1=sl;    
    P2=sum((sj).^[0:(d2-2)])/sum((sj).^[0:(d2-1)])*sj;
    G2=(sj).^(d2-1)/sum((sj).^[0:(d2-1)])*sj;%
    P3=sum((ss).^[0:(d3-2)])/sum((ss).^[0:(d3-1)])*ss;
    G3=(ss).^(d3-1)/sum((ss).^[0:(d3-1)])*ss;%
    P4=sa;
    F4=se*f;
A8=[0,0,0,F4;
    G1,P2,0,0;
    0, G2,P3,0;
    0,0,G3,P4];
[lambda8, v8, w8] = eigen(A8);
S8=v8(:)*w8(:)'/(v8(:)'*w8(:)); %Sensitivity Matrix
sens8(1)=S8(2,1);
sens8(2)=    S8(2,2)*(6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1)/(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6 + 6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1) +...
    S8(3,2)*(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6)/(sj^12 + 2*sj^11 + 3*sj^10 + 4*sj^9 + 5*sj^8 + 6*sj^7 + 7*sj^6 + 6*sj^5 + 5*sj^4 + 4*sj^3 + 3*sj^2 + 2*sj + 1);
sens8(3)=   S8(3,3)*(6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1)/(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6 + 6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1)+...
    S8(4,3)*(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6)/(ss^12 + 2*ss^11 + 3*ss^10 + 4*ss^9 + 5*ss^8 + 6*ss^7 + 7*ss^6 + 6*ss^5 + 5*ss^4 + 4*ss^3 + 3*ss^2 + 2*ss + 1);
sens8(4)=    S8(4,4);
sens8(5)=    S8(1,4)*se;

%% North population fixed duration with labmda AAS 3
A9=zeros(4,4);
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    G1=sl;
    P2=sum((sj/L).^[0:(d2-2)])/sum((sj/L).^[0:(d2-1)])*sj;
    G2=(sj/L).^(d2-1)/sum((sj/L).^[0:(d2-1)])*sj;%
    P3=sum((ss/L).^[0:(d3-2)])/sum((ss/L).^[0:(d3-1)])*ss;
    G3=(ss/L).^(d3-1)/sum((ss/L).^[0:(d3-1)])*ss;%
    P4=sa;
    F4=se*f;
    A9=[0,0,0,F4;
        G1,P2,0,0;
        0, G2,P3,0;
        0,0,G3,P4];
    L=max(real(eig(A9)));
end
[lambda9, v9, w9] = eigen(A9);
for PAR = 1:5
    for j=1:100
        L=1;
        for k=1:100 % Iterative method to seek convergence of lambda
            sl2=P(j,1,PAR);
            sj2=P(j,2,PAR);
            ss2=P(j,3,PAR);
            sa2=P(j,4,PAR);
            f2=P(j,5,PAR);
            G1=sl2;
            P2=sum((sj2/L).^[0:(d2-2)])/sum((sj2/L).^[0:(d2-1)])*sj2;
            G2=(sj2/L).^(d2-1)/sum((sj2/L).^[0:(d2-1)])*sj2;%
            P3=sum((ss2/L).^[0:(d3-2)])/sum((ss2/L).^[0:(d3-1)])*ss2;
            G3=(ss2/L).^(d3-1)/sum((ss2/L).^[0:(d3-1)])*ss2;%
            P4=sa2;
            F4=se*f2;
            A=[0,0,0,F4;
                G1,P2,0,0;
                0, G2,P3,0;
                0,0,G3,P4];
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


%% Elasticity
PARA=[sl,sj,ss,sa,f];
elas1=PARA.*sens1/lambda1;
elas2=PARA.*sens2/lambda2;
elas3=PARA.*sens3/lambda3;
elas4=PARA.*sens4/lambda4;
elas5=PARA.*sens5/lambda5;
elas6=PARA.*sens6/lambda6;
elas7=PARA.*sens7/lambda7;
elas8=PARA.*sens8/lambda8;
elas9=PARA.*sens9/lambda9;
% 
% 
D1=sort(eig(A1),'descend');
D2=sort(eig(A2),'descend');
D3=sort(eig(A3),'descend');
D4=sort(eig(A4),'descend');
D5=sort(eig(A5),'descend');
D6=sort(eig(A6),'descend');
D7=sort(eig(A7),'descend');
D8=sort(eig(A8),'descend');
D9=sort(eig(A9),'descend');
% 
DR(1)=D1(1)/abs(D1(2));
DR(2)=D2(1)/abs(D2(2));
DR(3)=D3(1)/abs(D3(2));
DR(4)=D4(1)/abs(D4(2));
DR(5)=D5(1)/abs(D5(2));
DR(6)=D6(1)/abs(D6(2));
DR(7)=D7(1)/abs(D7(2));
DR(8)=D8(1)/abs(D8(2));
DR(9)=D9(1)/abs(D9(2));
% 
% %% Generation Time
% 
GT(1)=lambda1*(v1*w1)/(v1(1)*A1(1,:)*w1);
GT(2)=lambda2*(v2*w2)/(v2(1)*A2(1,:)*w2);
GT(3)=lambda3*(v3*w3)/(v3(1)*A3(1,:)*w3);
GT(4)=lambda4*(v4*w4)/(v4(1)*A4(1,:)*w4);
GT(5)=lambda5*(v5*w5)/(v5(1)*A5(1,:)*w5);
GT(6)=lambda6*(v6*w6)/(v6(1)*A6(1,:)*w6);
GT(7)=lambda7*(v7*w7)/(v7(1)*A7(1,:)*w7);
GT(8)=lambda8*(v8*w8)/(v8(1)*A8(1,:)*w8);
GT(9)=lambda9*(v9*w9)/(v9(1)*A9(1,:)*w9);
% 
% 
tbl=table([lambda1;lambda2;lambda3;lambda4;lambda5;lambda6;lambda7;lambda8;lambda9],[w1';w2';w3';w4';w5';w6';w7';w8';w9'],[v1;v2;v3;v4;v5;v6;v7;v8;v9],[sens1;sens2;sens3;sens4;sens5;sens6;sens7;sens8;sens9],[elas1;elas2;elas3;elas4;elas5;elas6;elas7;elas8;elas9],GT',DR','VariableNames',{'Lambda','Stable_Stage_Distribution','Reproductive_Values','Sensitivity','Elasticity','Generation_Time','Damping_Ratio'})
writetable(tbl,'alligatorNorth.csv')
% 
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
end