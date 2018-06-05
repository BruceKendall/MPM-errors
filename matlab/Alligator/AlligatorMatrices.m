clear
%% North Population Origninal Matrix
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
f3=s4^(9/12)*s0*fn; %survival of adults/9months+ eggsurvial/3months x fecundity %
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
A2(1,4)=f3;
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
    A3(1,4)=f3;
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
A4(1,16)=f3;
%% South Population "original matrix"%
fs=5.98;
d2s=3;%South stage duration
d3s=3;%This was modified from the table to match the results in the paper
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
f4=s4^(9/12)*fs*s0; %survival of adults/9months+ eggsurvial/3months x fecundity x sex ratio%
P2=(1-s2^(d2s-1))/(1-s2^d2s)*s2;
G2=s2^d2s*(1-s2)/(1-s2^d2s);
P3=(1-s3^(d3s-1))/(1-s3^d3s)*s3;
G3=s3^d3s*(1-s3)/(1-s3^d3s);
A6=zeros(4,4);
A6(2,1)=s1;
A6(2,2)=P2;
A6(3,2)=G2;
A6(3,3)=P3;
A6(4,3)=G3;
A6(4,4)=s4;
A6(1,4)=f4;
%% South population fixed duration with lambda
A7=zeros(4,4);
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    P2=sum((s2/L).^[0:d2s-2])/sum((s2/L).^[0:d2s-1])*s2;
    G2=(s2/L).^(d2s-1)/sum((s2/L).^[0:d2s-1])*s2;%
    P3=sum((s3/L).^[0:d3s-2])/sum((s3/L).^[0:d3s-1])*s3;
    G3=(s3/L).^(d3s-1)/sum((s3/L).^[0:d3s-1])*s3;%
    A7(2,1)=s1;
    A7(2,2)=P2;
    A7(3,2)=G2;
    A7(3,3)=P3;
    A7(4,3)=G3;
    A7(4,4)=s4;
    A7(1,4)=f4;
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
A8(1,8)=f4;
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
axis([0.5,8.5,AX(3),AX(4)])
set(gca,'fontsize',12)
subplot1(2)
title('Juveniles','fontsize',12)
bar([w1(3),w2(2),w3(2),sum(w4(2:8)),w5(3),w6(2),w7(2),sum(w8(2:4))])
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gca,'fontsize',12)
subplot1(3)
title('SubAdults','fontsize',12)
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
% title('Hatclings','fontsize',12)
% bar([v1(2),v2(1),v3(1),v4(1),v5(2),v6(1),v7(1),v8(1)])
% AX=axis;
% axis([0.5,8.5,AX(3),AX(4)])
% subplot1(2)
title('Juveniles','fontsize',12)
bar([v1(3),v2(2),v3(2),mean(v4(2:8)),v5(3),v6(2),v7(2),mean(v8(2:4))])
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
subplot1(2)
title('SubAdults','fontsize',12)
bar([v1(4),v2(3),v3(3),mean(v4(9:15)),v5(4),v6(3),v7(3),mean(v8(5:7))])
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
 %% Sensitivity
figure(4)
S1=v1(:)*w1(:)'/(v1(:)'*w1(:)); %Sensitivity Matrix
sens1=[S1(2,2)*A1(2,2)/(A1(2,2)+A1(3,2))+S1(3,2)*A1(3,2)/(A1(2,2)+A1(3,2)),S1(3,3)*A1(3,3)/(A1(3,3)+A1(4,3))+S1(4,3)*A1(4,3)/(A1(3,3)+A1(4,3)),S1(4,4)*A1(4,4)/(A1(4,4)+A1(5,4))+S1(5,4)*A1(5,4)/(A1(4,4)+A1(5,4)),S1(5,5),S1(1,5)];
S2=v2(:)*w2(:)'/(v2(:)'*w2(:)); %Sensitivity Matrix
sens2=[S2(2,1),S2(2,2)*A2(2,2)/(A2(2,2)+A2(3,2))+S2(3,2)*A2(3,2)/(A2(2,2)+A2(3,2)),S2(3,3)*A2(3,3)/(A2(3,3)+A2(4,3))+S2(4,3)*A2(4,3)/(A2(3,3)+A2(4,3)),S2(4,4)+S2(1,4)*(9/12)*s4^(-3/12)*s0*fn,S2(1,4)*s4^(9/12)*s0];
% S3=v3(:)*w3(:)'/(v3(:)'*w3(:)); %Sensitivity Matrix
% sens3=[S3(2,1),S3(2,2)*A3(2,2)/(A3(2,2)+A3(3,2))+S3(3,2)*A3(3,2)/(A3(2,2)+A3(3,2)),S3(3,3)*A3(3,3)/(A3(3,3)+A3(4,3))+S3(4,3)*A3(4,3)/(A3(3,3)+A3(4,3)),S3(4,4),S3(1,4)];
sens3=numsens2;
S4=v4(:)*w4(:)'/(v4(:)'*w4(:)); %Sensitivity Matrix
sens4=[S4(2,1),S4(3,2)+S4(4,3)+S4(5,4)+S4(6,5)+S4(7,6)+S4(8,7)+S4(9,8),S4(10,9)+S4(11,10)+S4(12,11)+S4(13,12)+S4(14,13)+S4(15,14)+S4(16,15),S4(16,16),S4(1,16)];
S5=v5(:)*w5(:)'/(v5(:)'*w5(:)); %Sensitivity Matrix
sens5=[S5(2,2)*A5(2,2)/(A5(2,2)+A5(3,2))+S5(3,2)*A5(3,2)/(A5(2,2)+A5(3,2)),S5(3,3)*A5(3,3)/(A5(3,3)+A5(4,3))+S5(4,3)*A5(4,3)/(A5(3,3)+A5(4,3)),S5(4,4)*A5(4,4)/(A5(4,4)+A5(5,4))+S5(5,4)*A5(5,4)/(A5(4,4)+A5(5,4)),S5(5,5),S5(1,5)];
S6=v6(:)*w6(:)'/(v6(:)'*w6(:)); %Sensitivity Matrix
sens6=[S6(2,1),S6(2,2)*A6(2,2)/(A6(2,2)+A6(3,2))+S6(3,2)*A6(3,2)/(A6(2,2)+A6(3,2)),S6(3,3)*A6(3,3)/(A6(3,3)+A6(4,3))+S6(4,3)*A6(4,3)/(A6(3,3)+A6(4,3)),S6(4,4)+S6(1,4)*(9/12)*s4^(-3/12)*s0*fs,S6(1,4)*s4^(9/12)*s0];
% S7=v7(:)*w7(:)'/(v7(:)'*w7(:)); %Sensitivity Matrix
% sens7=[S7(2,1),S7(2,2)*A7(2,2)/(A7(2,2)+A7(3,2))+S7(3,2)*A7(3,2)/(A7(2,2)+A7(3,2)),S7(3,3)*A7(3,3)/(A7(3,3)+A7(4,3))+S7(4,3)*A7(4,3)/(A7(3,3)+A7(4,3)),S7(4,4),S7(1,4)];
sens7=numsens3;
S8=v8(:)*w8(:)'/(v8(:)'*w8(:)); %Sensitivity Matrix
sens8=[S8(2,1),S8(3,2)+S8(4,3)+S8(5,4),S8(6,5)+S8(7,6)+S8(8,7),S8(8,8),S8(1,8)];
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
