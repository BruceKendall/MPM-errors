clear
%% Parameters in the paper for the North population
s0=0.54;
s1=0.38;
s2=0.78;
s3=0.73;
s4=0.83;
fn=2.37;
d0=0.25;
d1n=1;
d2n=7;
d3n=7;

%% North Population Origninal Matrix
G2=s2^d2n*(1-s2)/(1-s2^d2n);
P2=(1-s2^(d2n-1))/(1-s2^d2n)*s2;
P3=(1-s3^(d3n-1))/(1-s3^d3n)*s3;
G3=s3^d3n*(1-s3)/(1-s3^d3n);
A1=zeros(5,5);
A1(2,1)=s0;
A1(3,2)=s1;
A1(4,3)=G2;
A1(3,3)=P2;
A1(5,4)=G3;
A1(4,4)=P3;
A1(5,5)=s4;
A1(1,5)=fn;
%% Correction of fertility and reduction of stage to four stage for north population
P2=(1-s2^(d2n-1))/(1-s2^d2n)*s2;
G2=s2^d2n*(1-s2)/(1-s2^d2n);
P3=(1-s3^(d3n-1))/(1-s3^d3n)*s3;
G3=s3^d3n*(1-s3)/(1-s3^d3n);
A2=zeros(4,4);
A2(2,1)=s1;
A2(2,2)=P2;
A2(3,2)=G2;
A2(3,3)=P3;
A2(4,3)=G3;
A2(4,4)=s4;
A2(1,3)=s3^(9/12)*s0*fn*s3^(d3n-1)*(1-s3)/(1-s3^d3n);
A2(1,4)=s4^(9/12)*s0*fn;

%% North population fixed duration with labmda
A3=zeros(4,4);
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    P2=sum((s2/L).^[0:(d2n-2)])/sum((s2/L).^[0:(d2n-1)])*s2;
    G2=(s2/L).^(d2n-1)/sum((s2/L).^[0:(d2n-1)])*s2;%
    P3=sum((s3/L).^[0:(d3n-2)])/sum((s3/L).^[0:(d3n-1)])*s3;
    G3=(s3/L).^(d3n-1)/sum((s3/L).^[0:(d3n-1)])*s3;%
    A3(2,1)=s1;
    A3(2,2)=P2;
    A3(3,2)=G2;
    A3(3,3)=P3;
    A3(4,3)=G3;
    A3(4,4)=s4;
    A3(1,3)=s3^(9/12)*s0*fn*(s3/L).^(d3n-1)/sum((s3/L).^[0:(d3n-1)]);
    A3(1,4)=s4^(9/12)*s0*fn;
    L=max(real(eig(A3)));
end
%% Leslie Matrx for North Population
A4=zeros(16,16);
A4(2,1)=s1;
A4(3,2)=s2;
A4(4,3)=s2;
A4(5,4)=s2;
A4(6,5)=s2;
A4(7,6)=s2;
A4(8,7)=s2;
A4(9,8)=s2;
A4(10,9)=s3;
A4(11,10)=s3;
A4(12,11)=s3;
A4(13,12)=s3;
A4(14,13)=s3;
A4(15,14)=s3;
A4(16,15)=s3;
A4(16,16)=s4;
A4(1,15)=s3^(9/12)*fn*s0;
A4(1,16)=s4^(9/12)*fn*s0;

%% South Population Parameters
fs=5.98;
d2s=3;%South stage duration
d3s=3;%This was modified from the table to match the results in the paper

%% South Population "original matrix"%
P2=(1-s2^(d2s-1))/(1-s2^d2s)*s2;
G2=s2^d2s*(1-s2)/(1-s2^d2s);
P3=(1-s3^(d3s-1))/(1-s3^d3s)*s3;
G3=s3^d3s*(1-s3)/(1-s3^d3s);
A5=zeros(5,5);
A5(2,1)=s0;
A5(3,2)=s1;
A5(4,3)=G2;
A5(3,3)=P2;
A5(5,4)=G3;
A5(4,4)=P3;
A5(5,5)=s4;
A5(1,5)=fs;

%% Correction of fertility and reduction of stage to four stage for south population
P2=(1-s2^(d2s-1))/(1-s2^d2s)*s2;
G2=s2^d2s*(1-s2)/(1-s2^d2s);
P3=(1-s3^(d3s-1))/(1-s3^d3s)*s3;
G3=s3^d3s*(1-s3)/(1-s3^d3s);
fs2=s3^(9/12)*fs*s0*s3^(d3s-1)*(1-s3)/(1-s3^d3s); %survival of adults/9months+ eggsurvial/3months x fecundity x sex ratio%
fs3=s4^(9/12)*fs*s0; %survival of adults/9months+ eggsurvial/3months x fecundity x sex ratio%
A6=zeros(4,4);
A6(2,1)=s1;
A6(2,2)=P2;
A6(3,2)=G2;
A6(3,3)=P3;
A6(4,3)=G3;
A6(4,4)=s4;
A6(1,3)=fs2;
A6(1,4)=fs3;

%% South population fixed duration with lambda
A7=zeros(4,4);
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    P2=sum((s2/L).^[0:(d2s-2)])/sum((s2/L).^[0:(d2s-1)])*s2;
    G2=(s2/L).^(d2s-1)/sum((s2/L).^[0:(d2s-1)])*s2;%
    P3=sum((s3/L).^[0:(d3s-2)])/sum((s3/L).^[0:(d3s-1)])*s3;
    G3=(s3/L).^(d3s-1)/sum((s3/L).^[0:(d3s-1)])*s3;%
    A7(2,1)=s1;
    A7(2,2)=P2;
    A7(3,2)=G2;
    A7(3,3)=P3;
    A7(4,3)=G3;
    A7(4,4)=s4;
    A7(1,3)=s3^(9/12)*fs*s0*(s3/L).^(d3s-1)/sum((s3/L).^[0:(d3s-1)]);
    A7(1,4)=s4^(9/12)*fs*s0;
    L=max(real(eig(A7)));
end

%% Leslie Matrix for South Population
A8=zeros(8,8);
A8(2,1)=s1;
A8(3,2)=s2;
A8(4,3)=s2;
A8(5,4)=s2;
A8(6,5)=s3;
A8(7,6)=s3;
A8(8,7)=s3;
A8(8,8)=s4;
A8(1,7)=s3^(9/12)*fs*s0;
A8(1,8)=s4^(9/12)*fs*s0;

%% Eigenvalue and Eigenvectors
[lambda1, v1, w1] = eigen(A1);
[lambda2, v2, w2] = eigen(A2);
[lambda3, v3, w3] = eigen(A3);
[lambda4, v4, w4] = eigen(A4);
[lambda5, v5, w5] = eigen(A5);
[lambda6, v6, w6] = eigen(A6);
[lambda7, v7, w7] = eigen(A7);
[lambda8, v8, w8] = eigen(A8);

%% Lambda%
LAM=[lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda8];
bar(LAM)
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6';'M7';'M8'])
xlabel('Model','fontsize',16)
ylabel('\lambda','fontsize',16)
AX=axis;
axis([0.5,8.5,0.6,1.1])

%% Stable stage distribution
figure(2)
subplot1(4,1,'Gap',[0.02 0.1])
subplot1(1)
title('Hatclings','fontsize',12)
bar([w1(2),w2(1),w3(1),w4(1),w5(2),w6(1),w7(1),w8(1)])
AX=axis;
axis([0.5,8.5,AX(3),0.45])
set(gca,'fontsize',12)
subplot1(2)
title('Juveniles','fontsize',12)
bar([w1(3),w2(2),w3(2),sum(w4(2:8)),w5(3),w6(2),w7(2),sum(w8(2:4))])
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gca,'fontsize',12)
subplot1(3)
title('Subadults','fontsize',12)
bar([w1(4),w2(3),w3(3),sum(w4(9:15)),w5(4),w6(3),w7(3),sum(w8(5:7))])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6,7,8])
ylabel('Proportion of Individuals in a Population','fontsize',12)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
subplot1(4)
title('Adults','fontsize',12)
bar([w1(5),w2(4),w3(4),w4(16),w5(5),w6(4),w7(4),w8(8)])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6';'M7';'M8'])
xlabel('Model','fontsize',12)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
w4v2=[w4(1),sum(w4(2:8)),sum(w4(9:15)),w4(16)];
w8v2=[w8(1),sum(w8(2:4)),sum(w8(5:7)),w8(8)];

%% Reproductive values
figure(3)
v1=v1/v1(2);
v2=v2/v2(1);
v3=v3/v3(1);
v4=v4/v4(1);
v5=v5/v5(2);
v6=v6/v6(1);
v7=v7/v7(1);
v8=v8/v8(1);
subplot1(3,1,'Gap',[0.02 0.1])
subplot1(1)
title('Juveniles','fontsize',12)
bar([v1(3),v2(2),v3(2),(v4(2:8)*w4(2:8)/(sum(w4(2:8)))),v5(3),v6(2),v7(2),(v8(2:4)*w8(2:4)/(sum(w8(2:4))))])
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
subplot1(2)
title('Subadults','fontsize',12)
bar([v1(4),v2(3),v3(3),v4(9:15)*w4(9:15)/sum(w4(9:15)),v5(4),v6(3),v7(3),(v8(5:7)*w8(5:7)/(sum(w8(5:7))))])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6,7,8])
ylabel('Reproductive Value','fontsize',12)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
subplot1(3)
title('Adults','fontsize',12)
bar([v1(5),v2(4),v3(4),v4(16),v5(5),v6(4),v7(4),v8(8)])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6';'M7';'M8'])
xlabel('Model','fontsize',12)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
v4v2=[v4(1),(v4(2:8)*w4(2:8)/(sum(w4(2:8)))),v4(9:15)*w4(9:15)/sum(w4(9:15)),v4(16)];
v8v2=[v8(1),(v8(2:4)*w8(2:4)/(sum(w8(2:4)))),(v8(5:7)*w8(5:7)/(sum(w8(5:7)))),v8(8)];
 
%% Sensitivity
figure(4)
sens1=numsensA1;
sens2=numsensA2;
sens3=numsensA3;
S4=v4(:)*w4(:)'/(v4(:)'*w4(:)); %Sensitivity Matrix
sens4=[S4(2,1),S4(3,2)+S4(4,3)+S4(5,4)+S4(6,5)+S4(7,6)+S4(8,7)+S4(9,8),S4(10,9)+S4(11,10)+S4(12,11)+S4(13,12)+S4(14,13)+S4(15,14)+S4(16,15)+S4(1,15)*9/12*s3^(-3/12)*s0*fn,S4(16,16)+S4(1,16)*9/12*s4^(-3/12)*s0*fn,    S4(1,15)*s3^(9/12)*s0+S4(1,16)*s4^(9/12)*s0];
sens5=numsensA5;
sens6=numsensA6;
sens7=numsensA7;
S8=v8(:)*w8(:)'/(v8(:)'*w8(:)); %Sensitivity Matrix
sens8=[S8(2,1), S8(3,2)+S8(4,3)+S8(5,4),S8(6,5)+S8(7,6)+S8(8,7)+S8(1,7)*9/12*s3^(-3/12)*s0*fs,   S8(8,8)+S8(1,8)*9/12*s4^(-3/12)*s0*fs,       S8(1,7)*s3^(9/12)*s0+S8(1,8)*s4^(9/12)*s0];

subplot1(8,1,'Gap',[0.02 0.05])
subplot1(1)
title('Model 1','fontsize',12)
bar(sens1)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(2)
title('Model 2','fontsize',12)
bar(sens2)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(3)
title('Model 3','fontsize',12)
bar(sens3)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(4)
title('Model 4','fontsize',12)
bar(sens4)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])
ylabel('Sensitivity','fontsize',12)

subplot1(5)
title('Model 5','fontsize',12)
bar(sens5)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])
      
subplot1(6)
title('Model 6','fontsize',12)
bar(sens6)
AX=axis;
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(7)
title('Model 7','fontsize',12)
bar(sens7)
AX=axis;
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(8)
title('Model 8','fontsize',12)
bar(sens8)
AX=axis;
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6],'xticklabel',['S_h';'S_j';'S_s';'S_a';'F_ '])
AX=axis;
axis([0.5,5.5,0,1])

set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')   

 %% Elasticity
figure(5)
PARA=[0.38,0.78,0.73,0.83,2.37];
elas1=PARA.*sens1/lambda1;
elas2=PARA.*sens2/lambda2;
elas3=PARA.*sens3/lambda3;
elas4=PARA.*sens4/lambda4;
PARA=[0.38,0.78,0.73,0.83,5.98];
elas5=PARA.*sens5/lambda5;
elas6=PARA.*sens6/lambda6;
elas7=PARA.*sens7/lambda7;
elas8=PARA.*sens8/lambda8;
subplot1(8,1,'Gap',[0.02 0.05])
subplot1(1)
title('Model 1','fontsize',12)
bar(elas1)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(2)
title('Model 2','fontsize',12)
bar(elas2)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(3)
title('Model 3','fontsize',12)
bar(elas3)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(4)
title('Model 4','fontsize',12)
bar(elas4)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])
ylabel('Elasticity','fontsize',12)

subplot1(5)
title('Model 5','fontsize',12)
bar(elas5)
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])
      
subplot1(6)
title('Model 6','fontsize',12)
bar(elas6)
AX=axis;
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(7)
title('Model 7','fontsize',12)
bar(elas7)
AX=axis;
set(gca,'fontsize',12,'xtick',[1,2,3,4,5])
AX=axis;
axis([0.5,5.5,0,1])

subplot1(8)
title('Model 8','fontsize',12)
bar(elas8)
AX=axis;
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6],'xticklabel',['S_h';'S_j';'S_s';'S_a';'F_ '])
AX=axis;
axis([0.5,5.5,0,1])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')   

%% Damping Ratio
figure (6)
D1=sort(eig(A1),'descend');
D2=sort(eig(A2),'descend');
D3=sort(eig(A3),'descend');
D4=sort(eig(A4),'descend');
D5=sort(eig(A5),'descend');
D6=sort(eig(A6),'descend');
D7=sort(eig(A7),'descend');
D8=sort(eig(A8),'descend');
DR(1)=D1(1)/abs(D1(2));
DR(2)=D2(1)/abs(D2(2));
DR(3)=D3(1)/abs(D3(2));
DR(4)=D4(1)/abs(D4(2));
DR(5)=D5(1)/abs(D5(2));
DR(6)=D6(1)/abs(D6(2));
DR(7)=D7(1)/abs(D7(2));
DR(8)=D8(1)/abs(D8(2));
bar(DR)
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6';'M7';'M8'])
xlabel('Model','fontsize',16)
ylabel('Damping Ratio','fontsize',16)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])

%% Generation Time
figure (7)
GT(1)=lambda1*(v1*w1)/(v1(1)*A1(1,:)*w1);
GT(2)=lambda2*(v2*w2)/(v2(1)*A2(1,:)*w2);
GT(3)=lambda3*(v3*w3)/(v3(1)*A3(1,:)*w3);
GT(4)=lambda4*(v4*w4)/(v4(1)*A4(1,:)*w4);
GT(5)=lambda5*(v5*w5)/(v5(1)*A5(1,:)*w5);
GT(6)=lambda6*(v6*w6)/(v6(1)*A6(1,:)*w6);
GT(7)=lambda7*(v7*w7)/(v7(1)*A7(1,:)*w7);
GT(8)=lambda8*(v8*w8)/(v8(1)*A8(1,:)*w8);
bar(GT)
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6';'M7';'M8'])
xlabel('Model','fontsize',16)
ylabel('Generation Time','fontsize',16)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])

figure (1)
print -dtiff -r300 alligatorlambda.tif
figure(2)
print -dtiff -r300 alligatorSSD.tif
figure(3)
print -dtiff -r300 alligatorRV.tif
figure(4) 
print -dtiff -r300 alligatorSensitivity.tif
figure(5) 
print -dtiff -r300 alligatorElasticity.tif
figure(6) 
print -dtiff -r300 alligatorDR.tif
figure(7) 
print -dtiff -r300 alligatorGT.tif

tbl=table([lambda1;lambda2;lambda3;lambda4;lambda5;lambda6;lambda7;lambda8],[w1(1:end)';0,w2';0,w3';0,w4v2;w5(1:end)';0, w6';0, w7';0, w8v2],[v1(1:end);0,v2;0, v3;0, v4v2;v5(1:end);0,v6;0,v7;0,v8v2],[sens1;sens2;sens3;sens4;sens5;sens6;sens7;sens8],[elas1;elas2;elas3;elas4;elas5;elas6;elas7;elas8],GT',DR','VariableNames',{'Lambda','Stable_Stage_Distribution','Reproductive_Values','Sensitivity','Elasticity','Generation_Time','Damping_Ratio'})
writetable(tbl,'alligator.csv')
 
