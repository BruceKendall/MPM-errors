%% Parameters in the paper
Ml=0.350;
Ma=0.052;
Mj=0.165;
rho=0.46;
Dl=30;
Me=0.310;
f=194577; % Strictly speaking there is a mistake in this calculation, but we probably keep it as is. 
De=3;
%% Original matrix (Table 1)
PA=exp(-Ma);
PJ=(11/12)*exp(-Mj);
GL=exp(-Ml*Dl);
GJ=(1/12)*exp(-Mj);
RA=rho*f*exp(-Me*De);
A1=[0,0,RA;GL,PJ,0;0,GJ,PA];
%% Original matrix with adult survival incorporated into fertility
PA=exp(-Ma);
PJ=(11/12)*exp(-Mj);
GL=exp(-Ml*Dl);
GJ=(1/12)*exp(-Mj);
RA=rho*f*exp(-Me*De)*exp(-Ma*27/30); % This has been modified to include adult survival
A2=[0,0,RA;GL,PJ,0;0,GJ,PA];
%% Original matrix with correction of transition rates (Geometric Distribution)
PA=exp(-Ma);
PJ=(9/10)*exp(-Mj);% This has been modified (duration in juvenile should be 10 months)
GL=exp(-Ml*Dl);
GJ=(1/10)*exp(-Mj);% This has been modified (duration in juvenile should be 10 months)
RA=rho*f*exp(-Me*De); 
A3=[0,0,RA;GL,PJ,0;0,GJ,PA];
%% Matching proportion (Fixed Stage Duration assuming Lambda = 1)
PA=exp(-Ma);
PJ=sum(exp(-Mj).^[0:8])/sum(exp(-Mj).^[0:9])*exp(-Mj);% This has been modified to match the proportion transitioning
GL=exp(-Ml*Dl);
GJ=exp(-Mj).^9/sum(exp(-Mj).^[0:9])*exp(-Mj);% This has been modified to match the proportion transitioning
RA=rho*f*exp(-Me*De)*exp(-Ma*27/30);% This has been modified to include adult survival
A4=[0,0,RA;GL,PJ,0;0,GJ,PA];
%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda)
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    A5=[0,0, rho*f*exp(-Me*De)*exp(-Ma*27/30);exp(-Ml*Dl),sum((exp(-Mj)/L).^[0:8])/sum((exp(-Mj)/L).^[0:9])*exp(-Mj),0;0,(exp(-Mj)/L).^9/sum((exp(-Mj)/L).^[0:9])*exp(-Mj),exp(-Ma)];
    L=max(real(eig(A5)));
end
%% Leslie Matrix (This assumes it takes 12 months to start reproducing)
A6=zeros(12,12);
A6(2,1)=exp(-Ml*Dl);
A6(3,2)=exp(-Mj);
A6(4,3)=exp(-Mj);
A6(5,4)=exp(-Mj);
A6(6,5)=exp(-Mj);
A6(7,6)=exp(-Mj);
A6(8,7)=exp(-Mj);
A6(9,8)=exp(-Mj);
A6(10,9)=exp(-Mj);
A6(11,10)=exp(-Mj);
A6(12,11)=exp(-Mj);
A6(12,12)=exp(-Ma);
A6(1,12)=rho*f*exp(-Me*De)*exp(-Ma*27/30);
%% Original matrix with correction of transition rates (Geometric Distribution)
PA=exp(-Ma);
PJ=(10/11)*exp(-Mj);% This has been modified (duration in juvenile should be 11 months)
GL=exp(-Ml*Dl);
GJ=(1/11)*exp(-Mj);% This has been modified (duration in juvenile should be 11months)
RJ=GJ*RA;
RA=rho*f*exp(-Me*De); 
A7=[0,RJ,RA;GL,PJ,0;0,GJ,PA];
%% Matching proportion (Fixed Stage Duration assuming Lambda = 1)
PA=exp(-Ma);
PJ=sum(exp(-Mj).^[0:9])/sum(exp(-Mj).^[0:10])*exp(-Mj);% This has been modified to match the proportion transitioning
GL=exp(-Ml*Dl);
GJ=exp(-Mj).^10/sum(exp(-Mj).^[0:10])*exp(-Mj);% This has been modified to match the proportion transitioning
RA=rho*f*exp(-Me*De)*exp(-Ma*27/30);% This has been modified to include adult survival
RJ=RA*GJ;
A8=[0,RJ,RA;GL,PJ,0;0,GJ,PA];
%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda)
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    A9=[0,                    rho*f*exp(-Me*De)*exp(-Ma*27/30)*(exp(-Mj)/L).^10/sum((exp(-Mj)/L).^[0:10])*exp(-Mj),                        rho*f*exp(-Me*De)*exp(-Ma*27/30);
        exp(-Ml*Dl),        sum((exp(-Mj)/L).^[0:9])/sum((exp(-Mj)/L).^[0:10])*exp(-Mj),                                                                  0;
        0,                       (exp(-Mj)/L).^10/sum((exp(-Mj)/L).^[0:10])*exp(-Mj),                                                                             exp(-Ma)];
    L=max(real(eig(A9)));
end
%% Eigenvalue and Eigenvectors
[lambda1, v1, w1] = eigen(A1);
[lambda2, v2, w2] = eigen(A2);
%  [lambda3, v3, w3] = eigen(A3);
%  [lambda4, v4, w4] = eigen(A4);
%  [lambda5, v5, w5] = eigen(A5);
[lambda6, v6, w6] = eigen(A6);
[lambda3, v3, w3] = eigen(A7);
[lambda4, v4, w4] = eigen(A8);
[lambda5, v5, w5] = eigen(A9);
%% Lambda
LAM=[lambda1,lambda2,lambda3,lambda4,lambda5,lambda6];
bar(LAM)
axis([0.5,6.5,0.9,1.2])
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6'])
xlabel('Model','fontsize',16)
ylabel('\lambda','fontsize',16)
%% Stable stage distribution
figure(2)
subplot1(2,1,'Gap',[0.02 0.1])
subplot1(1)
% title('Larvae','fontsize',12)
% bar([w1(1),w2(1),w3(1),w4(1),w5(1)])
% axis([0.5,5.5,0.999,1.001])
% set(gca,'fontsize',12)
% subplot1(2)
title('Juveniles','fontsize',12)
bar([w1(2),w2(2),w3(2),w4(2),w5(2),sum(w6(2:11))])
AX=axis;
axis([0.5,6.5,AX(3),AX(4)])
set(gca,'fontsize',12)
ylabel('Proportion of Individuals in a Popualtion','fontsize',12)
subplot1(2)
title('Adults','fontsize',12)
bar([w1(3),w2(3),w3(3),w4(3),w5(3),sum(w6(12))])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6'])
xlabel('Model','fontsize',12)
ylabel('Proportion of Individuals in a Popualtion','fontsize',12)
AX=axis;
axis([0.5,6.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
 %% Reproductive values
figure(3)
v1=v1/v1(1);
v2=v2/v2(1);
v3=v3/v3(1);
v4=v4/v4(1);
v5=v5/v5(1);
v6=v6/v6(1);
subplot1(2,1,'Gap',[0.02 0.1])
subplot1(1)
% title('Larvae','fontsize',12)
% bar([w1(1),w2(1),w3(1),w4(1),w5(1)])
% axis([0.5,5.5,0.999,1.001])
% set(gca,'fontsize',12)
% subplot1(2)
title('Juveniles','fontsize',12)
bar([v1(2),v2(2),v3(2),v4(2),v5(2),mean(v6(2:11))])
AX=axis;
axis([0.5,6.5,AX(3),AX(4)])
set(gca,'fontsize',12)
ylabel('Reproductive value','fontsize',12)
subplot1(2)
title('Adults','fontsize',12)
bar([v1(3),v2(3),v3(3),v4(3),v5(3),v6(12)])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5,6],'xticklabel',['M1';'M2';'M3';'M4';'M5';'M6'])
xlabel('Model','fontsize',12)
ylabel('Reproductive value','fontsize',12)
AX=axis;
axis([0.5,6.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
 %% Sensitivity
figure(4)
S1=v1(:)*w1(:)'/(v1(:)'*w1(:)); %Sensitivity Matrix
sens1=[S1(2,1),S1(2,2)*A1(2,2)/(A1(2,2)+A1(3,2))+S1(3,2)*A1(3,2)/(A1(2,2)+A1(3,2)),S1(3,3),S1(1,3)*rho*exp(-Me*De)];
S2=v2(:)*w2(:)'/(v2(:)'*w2(:)); %Sensitivity Matrix
sens2=[S2(2,1),S2(2,2)*A2(2,2)/(A2(2,2)+A2(3,2))+S2(3,2)*A2(3,2)/(A2(2,2)+A2(3,2)),S2(3,3)+S2(1,3)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S2(1,3)*rho*exp(-Me*De)*exp(-Ma*27/30)];
S3=v3(:)*w3(:)'/(v3(:)'*w3(:)); %Sensitivity Matrix
sens3=[S3(2,1),S3(2,2)*A3(2,2)/(A3(2,2)+A3(3,2))+S3(3,2)*A3(3,2)/(A3(2,2)+A3(3,2)),S3(3,3),S3(1,3)*rho*exp(-Me*De)];
S4=v4(:)*w4(:)'/(v4(:)'*w4(:)); %Sensitivity Matrix
sens4=[S4(2,1),S4(2,2)*A4(2,2)/(A4(2,2)+A4(3,2))+S4(3,2)*A4(3,2)/(A4(2,2)+A4(3,2)),S4(3,3)+S4(1,3)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S4(1,3)*rho*exp(-Me*De)*exp(-Ma*27/30)];
% S5=v5(:)*w5(:)'/(v5(:)'*w5(:)); %Sensitivity Matrix
% sens5=[S5(2,1),S5(2,2)*A5(2,2)/(A5(2,2)+A5(3,2))+S5(3,2)*A5(3,2)/(A5(2,2)+A5(3,2)),S5(3,3)+S5(1,3)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S5(1,3)*rho*exp(-Me*De)*exp(-Ma*27/30)];
sens5=numsens2;
S6=v6(:)*w6(:)'/(v6(:)'*w6(:)); %Sensitivity Matrix
sens6=[S6(2,1),S6(3,2)+S6(4,3)+S6(5,4)+S6(6,5)+S6(7,6)+S6(8,7)+S6(9,8)+S6(10,9)+S6(11,10)+S6(12,11),S6(12,12)+S6(1,12)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S6(1,12)*rho*exp(-Me*De)*exp(-Ma*27/30)];% I am not sure if we could just sum??
subplot1(6,1,'Gap',[0.02 0.05])
subplot1(1)
title('Model 1','fontsize',12)
bar(sens1)
set(gca,'fontsize',12,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])

subplot1(2)
title('Model 2','fontsize',12)
bar(sens2)
set(gca,'fontsize',12,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])

subplot1(3)
title('Model 3','fontsize',12)
bar(sens3)
set(gca,'fontsize',12,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])
ylabel('Sensitivity','fontsize',12)

subplot1(4)
title('Model 4','fontsize',12)
bar(sens4)
AX=axis;
axis([0.5,4.5,10^-10,10^10])
set(gca,'fontsize',12,'YScale','log')

subplot1(5)
title('Model 5','fontsize',12)
bar(sens5)
set(gca,'fontsize',12,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')       
      
subplot1(6)
title('Model 6','fontsize',12)
bar(sens6)
AX=axis;
axis([0.5,4.5,10^-10,10^10])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5],'xticklabel',['S_L';'S_J';'S_A';'f  '],'YScale','log')
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')       
      
 %% Elasticity
figure(5)
PARA=[exp(-Ml*Dl),exp(-Mj),exp(-Ma),f];
elas1=PARA.*sens1/lambda1;
elas2=PARA.*sens2/lambda2;
elas3=PARA.*sens3/lambda3;
elas4=PARA.*sens4/lambda4;
elas5=PARA.*sens5/lambda5;
elas6=PARA.*sens6/lambda6;
subplot1(6,1,'Gap',[0.02 0.05])
subplot1(1)
title('Model 1','fontsize',12)
bar(elas1)
axis([0.5,4.5,0,1])
set(gca,'fontsize',12)

subplot1(2)
title('Model 2','fontsize',12)
bar(elas2)
axis([0.5,4.5,0,1])
set(gca,'fontsize',12)

subplot1(3)
title('Model 3','fontsize',12)
bar(elas3)
axis([0.5,4.5,0,1])
set(gca,'fontsize',12)
ylabel('Elasticity','fontsize',12)

subplot1(4)
title('Model 4','fontsize',12)
bar(elas4)
axis([0.5,4.5,0,1])
set(gca,'fontsize',12)

subplot1(5)
title('Model 5','fontsize',12)
bar(elas5)
axis([0.5,4.5,0,1])
set(gca,'fontsize',12)
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')       
      
subplot1(6)
title('Model 6','fontsize',12)
bar(elas6)
axis([0.5,4.5,0,1])
set(gca,'fontsize',12,'xtick',[1,2,3,4,5],'xticklabel',['S_L';'S_J';'S_A';'f  '])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')       

figure (1)
print -dtiff -r300 lionfishlambda.tif
figure(2)
print -dtiff -r300 lionfishSSD.tif
figure(3)
print -dtiff -r300 lionfishRV.tif
figure(4) 
print -dtiff -r300 lionfishSensitivity.tif
figure(5) 
print -dtiff -r300 lionfishElasticity.tif

