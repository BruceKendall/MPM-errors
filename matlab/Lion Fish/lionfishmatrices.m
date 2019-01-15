clear
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

%% Original matrix with adult survival incorporated into fertility and correction of transition rates
PA=exp(-Ma);
PJ=(10/11)*exp(-Mj);
GL=exp(-Ml*Dl);
GJ=(1/11)*exp(-Mj);
RA=rho*f*exp(-Me*De)*exp(-Ma*27/30); % This has been modified to include adult survival
A2=[0,0,RA;GL,PJ,0;0,GJ,PA];

%%  Original matrix with adult survival incorporated into fertility and correction of transition rates with inclusion of fertility for juveniles
PA=exp(-Ma);
PJ=(10/11)*exp(-Mj);% This has been modified (duration in juvenile should be 11 months)
GL=exp(-Ml*Dl);
GJ=(1/11)*exp(-Mj);% This has been modified (duration in juvenile should be 11months)
RA=rho*f*exp(-Me*De)*exp(-Ma*27/30);
RJ=(1/11)*rho*f*exp(-Me*De)*exp(-Mj*27/30);
A3=[0,RJ,RA;GL,PJ,0;0,GJ,PA];

%% Matching proportion (Fixed Stage Duration assuming Lambda = 1)
PA=exp(-Ma);
PJ=sum(exp(-Mj).^[0:9])/sum(exp(-Mj).^[0:10])*exp(-Mj);% This has been modified to match the proportion transitioning
GL=exp(-Ml*Dl);
GJ=exp(-Mj).^10/sum(exp(-Mj).^[0:10])*exp(-Mj);% This has been modified to match the proportion transitioning
RA=rho*f*exp(-Me*De)*exp(-Ma*27/30);% This has been modified to include adult survival
RJ=exp(-Mj).^10/sum(exp(-Mj).^[0:10])*rho*f*exp(-Me*De)*exp(-Mj*27/30);
A4=[0,RJ,RA;GL,PJ,0;0,GJ,PA];

%% Matching Proportion with lambda (Fixed Stage Duration with discounting with Lambda)
L=1;
for k=1:1000 % Iterative method to seek convergence of lambda
    A5=[0,                    rho*f*exp(-Me*De)*exp(-Mj*27/30)*(exp(-Mj)/L).^10/sum((exp(-Mj)/L).^[0:10]),                        rho*f*exp(-Me*De)*exp(-Ma*27/30);
        exp(-Ml*Dl),        sum((exp(-Mj)/L).^[0:9])/sum((exp(-Mj)/L).^[0:10])*exp(-Mj),                                                                  0;
        0,                       (exp(-Mj)/L).^10/sum((exp(-Mj)/L).^[0:10])*exp(-Mj),                                                                             exp(-Ma)];
    L=max(real(eig(A5)));
end
%% Leslie Matrix (This assumes it takes 12 months to start reproducing) without survival on fertility
A6=zeros(13,13);
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
A6(13,12)=exp(-Mj);
A6(13,13)=exp(-Ma);
A6(1,12)=rho*f*exp(-Me*De);
A6(1,13)=rho*f*exp(-Me*De);

%% Leslie Matrix (This assumes it takes 12 months to start reproducing) w/o reproduction on 13th class
A7=zeros(13,13);
A7(2,1)=exp(-Ml*Dl);
A7(3,2)=exp(-Mj);
A7(4,3)=exp(-Mj);
A7(5,4)=exp(-Mj);
A7(6,5)=exp(-Mj);
A7(7,6)=exp(-Mj);
A7(8,7)=exp(-Mj);
A7(9,8)=exp(-Mj);
A7(10,9)=exp(-Mj);
A7(11,10)=exp(-Mj);
A7(12,11)=exp(-Mj);
A7(13,12)=exp(-Mj);
A7(13,13)=exp(-Ma);
A7(1,13)=rho*f*exp(-Me*De)*exp(-Ma*27/30);

%% Leslie Matrix (This assumes it takes 12 months to start reproducing)
A8=zeros(13,13);
A8(2,1)=exp(-Ml*Dl);
A8(3,2)=exp(-Mj);
A8(4,3)=exp(-Mj);
A8(5,4)=exp(-Mj);
A8(6,5)=exp(-Mj);
A8(7,6)=exp(-Mj);
A8(8,7)=exp(-Mj);
A8(9,8)=exp(-Mj);
A8(10,9)=exp(-Mj);
A8(11,10)=exp(-Mj);
A8(12,11)=exp(-Mj);
A8(13,12)=exp(-Mj);
A8(13,13)=exp(-Ma);
A8(1,12)=rho*f*exp(-Me*De)*exp(-Mj*27/30);
A8(1,13)=rho*f*exp(-Me*De)*exp(-Ma*27/30);

%% Eigenvalue and Eigenvectors
[lambda1, v1, w1] = eigen(A1);
[lambda2, v2, w2] = eigen(A2);
[lambda3, v3, w3] = eigen(A3);
[lambda4, v4, w4] = eigen(A4);
[lambda5, v5, w5] = eigen(A5);
[lambda6, v6, w6] = eigen(A6);
[lambda7, v7, w7] = eigen(A7);
[lambda8, v8, w8] = eigen(A8);

%% Lambda
LAM=[lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda8];
bar(LAM)
axis([0.5,8.5,0.9,1.2])
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['L1';'L2';'L3';'L4';'L5';'L6';'L7';'L8'])
xlabel('Model','fontsize',16)
ylabel('\lambda','fontsize',16)
set(gcf,'PaperPosition',[0.25,0.25,10.5,8],'PaperSize',[8.5,11],'PaperOrientation','landscape')   

%% Stable stage distribution
figure(2)
subplot1(2,1,'Gap',[0.02 0.1])
subplot1(1)
title('Juveniles','fontsize',16)
bar([w1(2),w2(2),w3(2),w4(2),w5(2),sum(w6(2:12)),sum(w7(2:12)),sum(w8(2:12))])
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gca,'fontsize',16)
ylabel('Proportion of Individuals','fontsize',16)
subplot1(2)
title('Adults','fontsize',16)
bar([w1(3),w2(3),w3(3),w4(3),w5(3),w6(13),w7(13),w8(13)])
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['L1';'L2';'L3';'L4';'L5';'L6';'L7';'L8'])
xlabel('Model','fontsize',16)
ylabel('Proportion of Individuals','fontsize',16)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
% set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
% set(gcf,'PaperPosition',[0.25,0.25,10.5,8],'PaperSize',[8.5,11],'PaperOrientation','portrait')   
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
w6v2=[w6(1),sum(w6(2:12)),w6(13)];
w7v2=[w7(1),sum(w7(2:12)),w7(13)];
w8v2=[w8(1),sum(w8(2:12)),w8(13)];
%% Reproductive values
figure(3)
v1=v1/v1(1);
v2=v2/v2(1);
v3=v3/v3(1);
v4=v4/v4(1);
v5=v5/v5(1);
v6=v6/v6(1);
v7=v7/v7(1);
v8=v8/v8(1);
subplot1(2,1,'Gap',[0.02 0.1])
subplot1(1)
title('Juveniles','fontsize',16)
bar([v1(2),v2(2),v3(2),v4(2),v5(2),v6(2:12)*w6(2:12)/sum(w6(2:12)),v7(2:12)*w7(2:12)/sum(w7(2:12)),v8(2:12)*w8(2:12)/sum(w8(2:12))])
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gca,'fontsize',16)
ylabel('Reproductive Value','fontsize',16)
subplot1(2)
title('Adults','fontsize',16)
bar([v1(3),v2(3),v3(3),v4(3),v5(3),v6(13),v7(13),v8(13)])
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['L1';'L2';'L3';'L4';'L5';'L6';'L7';'L8'])
xlabel('Model','fontsize',16)
ylabel('Reproductive Value','fontsize',16)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,8,10.5],'PaperSize',[8.5,11],'PaperOrientation','portrait')
 v6v2=[v6(1),v6(2:12)*w6(2:12)/sum(w6(2:12)),v6(13)];
  v7v2=[v7(1),v7(2:12)*w7(2:12)/sum(w7(2:12)),v7(13)];
   v8v2=[v8(1),v8(2:12)*w8(2:12)/sum(w8(2:12)),v8(13)];
%% Sensitivity
figure(4)
S1=v1(:)*w1(:)'/(v1(:)'*w1(:)); %Sensitivity Matrix
sens1=[S1(2,1),S1(2,2)*A1(2,2)/(A1(2,2)+A1(3,2))+S1(3,2)*A1(3,2)/(A1(2,2)+A1(3,2)),S1(3,3),S1(1,3)*rho*exp(-Me*De)];
S2=v2(:)*w2(:)'/(v2(:)'*w2(:)); %Sensitivity Matrix
sens2=[S2(2,1),S2(2,2)*A2(2,2)/(A2(2,2)+A2(3,2))+S2(3,2)*A2(3,2)/(A2(2,2)+A2(3,2)),  S2(3,3)+S2(1,3)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),   S2(1,3)*rho*exp(-Me*De)*exp(-Ma*27/30)];
S3=v3(:)*w3(:)'/(v3(:)'*w3(:)); %Sensitivity Matrix
sens3=[S3(2,1),S3(2,2)*A3(2,2)/(A3(2,2)+A3(3,2))+S3(3,2)*A3(3,2)/(A3(2,2)+A3(3,2))+S3(1,2)/11*rho*f*exp(-Me*De)*27/30*exp(Mj*3/30),S3(3,3)+S3(1,3)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S3(1,2)*rho*exp(-Me*De)*exp(-Mj*27/30)/11+S3(1,3)*rho*exp(-Me*De)*exp(-Ma*27/30)];
S4=v4(:)*w4(:)'/(v4(:)'*w4(:)); %Sensitivity Matrix
sens4=numsens4; % Numberical calculation of sensitivity
sens5=numsens5; % Numberical calculation of sensitivity
S6=v6(:)*w6(:)'/(v6(:)'*w6(:)); %Sensitivity Matrix
sens6=[S6(2,1),S6(3,2)+S6(4,3)+S6(5,4)+S6(6,5)+S6(7,6)+S6(8,7)+S6(9,8)+S6(10,9)+S6(11,10)+S6(12,11)+S6(13,12),S6(13,13),S6(1,12)*rho*exp(-Me*De)+S6(1,13)*rho*exp(-Me*De)];
S7=v7(:)*w7(:)'/(v7(:)'*w7(:)); %Sensitivity Matrix
sens7=[S7(2,1),S7(3,2)+S7(4,3)+S7(5,4)+S7(6,5)+S7(7,6)+S7(8,7)+S7(9,8)+S7(10,9)+S7(11,10)+S7(12,11)+S7(13,12),S7(13,13)+S7(1,13)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S7(1,13)*rho*exp(-Me*De)*exp(-Ma*27/30)];
S8=v8(:)*w8(:)'/(v8(:)'*w8(:)); %Sensitivity Matrix
sens8=[S8(2,1),S8(3,2)+S8(4,3)+S8(5,4)+S8(6,5)+S8(7,6)+S8(8,7)+S8(9,8)+S8(10,9)+S8(11,10)+S8(12,11)+S8(13,12)+S8(1,12)*rho*f*exp(-Me*De)*27/30*exp(Mj*3/30),S8(13,13)+S8(1,13)*rho*f*exp(-Me*De)*27/30*exp(Ma*3/30),S8(1,12)*rho*exp(-Me*De)*exp(-Mj*27/30)+S8(1,13)*rho*exp(-Me*De)*exp(-Ma*27/30)];

subplot1(8,1,'Gap',[0.02 0.05])

subplot1(1)
title('L1','fontsize',16)
bar(sens1)
set(gca,'fontsize',16,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])

subplot1(2)
title('L2','fontsize',16)
bar(sens2)
set(gca,'fontsize',16,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])

subplot1(3)
title('L3','fontsize',16)
bar(sens3)
set(gca,'fontsize',16,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])


subplot1(4)
title('L4','fontsize',16)
bar(sens4)
AX=axis;
axis([0.5,4.5,10^-10,10^10])
set(gca,'fontsize',16,'YScale','log')

subplot1(5)
title('L5','fontsize',16)
bar(sens5)
set(gca,'fontsize',16,'YScale','log')
AX=axis;
ylabel('Sensitivity','fontsize',16)
axis([0.5,4.5,10^-10,10^10])
    
      
subplot1(6)
title('L6','fontsize',16)
bar(sens6)
set(gca,'fontsize',16,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])

subplot1(7)
title('L7','fontsize',16)
bar(sens7)
set(gca,'fontsize',16,'YScale','log')
AX=axis;
axis([0.5,4.5,10^-10,10^10])

subplot1(8)
title('L8','fontsize',16)
bar(sens8)
AX=axis;
axis([0.5,4.5,10^-10,10^10])
xlabel('Parameters','fontsize',16)
set(gca,'fontsize',16,'xtick',[1,2,3,4,5],'xticklabel',['S_L';'S_J';'S_A';'f  '],'YScale','log')
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
elas7=PARA.*sens7/lambda7;
elas8=PARA.*sens8/lambda8;

subplot1(8,1,'Gap',[0.02 0.05])
subplot1(1)
title('L1','fontsize',16)
bar(elas1)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)

subplot1(2)
title('L2','fontsize',16)
bar(elas2)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)

subplot1(3)
title('L3','fontsize',16)
bar(elas3)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)


subplot1(4)
title('L4','fontsize',16)
bar(elas4)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)

subplot1(5)
title('L5','fontsize',16)
bar(elas5)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)
ylabel('Elasticity','fontsize',16)
     
subplot1(6)
title('L6','fontsize',16)
bar(elas6)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)

subplot1(7)
title('L7','fontsize',16)
bar(elas7)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16)

subplot1(8)
title('L8','fontsize',16)
bar(elas8)
axis([0.5,4.5,0,1])
set(gca,'fontsize',16,'xtick',[1,2,3,4,5],'xticklabel',['S_L';'S_J';'S_A';'f  '])
xlabel('Parameters','fontsize',16)

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
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['L1';'L2';'L3';'L4';'L5';'L6';'L7';'L8'])
xlabel('Model','fontsize',16)
ylabel('Damping Ratio','fontsize',16)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,10.5,8],'PaperSize',[8.5,11],'PaperOrientation','landscape')   

%% Generation Time
figure (7)
GT(1)=lambda1*(v1*w1)/(v1(1)*A1(1,:)*w1);
GT(2)=lambda2*(v2*w2)/(v2(1)*A2(1,:)*w2);
GT(3)=lambda3*(v3*w3)/(v3(1)*A3(1,:)*w3);
GT(4)=lambda4*(v4*w4)/(v4(1)*A4(1,:)*w4);
GT(5)=lambda5*(v5*w5)/(v5(1)*A5(1,:)*w5);
GT(6)=lambda6*(v6*w6)/(v6(1)*A8(1,:)*w6);
GT(7)=lambda7*(v7*w7)/(v7(1)*A7(1,:)*w7);
GT(8)=lambda8*(v8*w8)/(v8(1)*A8(1,:)*w8);

bar(GT)
set(gca,'fontsize',16,'xtick',[1,2,3,4,5,6,7,8],'xticklabel',['L1';'L2';'L3';'L4';'L5';'L6';'L7';'L8'])
xlabel('Model','fontsize',16)
ylabel('Generation Time','fontsize',16)
AX=axis;
axis([0.5,8.5,AX(3),AX(4)])
set(gcf,'PaperPosition',[0.25,0.25,10.5,8],'PaperSize',[8.5,11],'PaperOrientation','landscape')   

% figure (1)
% print -dtiff -r300 lionfishlambda.tif
% figure(2)
% print -dtiff -r300 lionfishSSD.tif
% figure(3)
% print -dtiff -r300 lionfishRV.tif
% figure(4) 
% print -dtiff -r300 lionfishSensitivity.tif
% figure(5) 
% print -dtiff -r300 lionfishElasticity.tif
% figure(6) 
% print -dtiff -r300 lionfishDR.tif
% figure(7) 
% print -dtiff -r300 lionfishGT.tif
% tbl=table([lambda1;lambda2;lambda3;lambda4;lambda5;lambda6;lambda7;lambda8],[w1';w2';w3';w4';w5';w6v2;w7v2;w8v2],[v1;v2;v3;v4;v5;v6v2;v7v2;v8v2],[sens1;sens2;sens3;sens4;sens5;sens6;sens7;sens8],[elas1;elas2;elas3;elas4;elas5;elas6;elas7;elas8],GT',DR','VariableNames',{'Lambda','Stable_Stage_Distribution','Reproductive_Values','Sensitivity','Elasticity','Generation_Time','Damping_Ratio'})
% writetable(tbl,'lionfish.csv')
