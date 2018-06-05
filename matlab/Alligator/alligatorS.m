s0=0.54;
s1=0.38;
s2=0.78;
s3=0.73;
s4=0.83;
f=5.98;
d0=0.25;
d1=1;
d2=3;
d3=3;
P2=(1-s2^(d2-1))/(1-s2^d2)*s2;
G2=s2^d2*(1-s2)/(1-s2^d2);
P3=(1-s3^(d3-1))/(1-s3^d3)*s3;
G3=s3^d3*(1-s3)/(1-s3^d3);
A=zeros(5,5);
A(2,1)=s0;
A(3,2)=s1;
A(4,3)=G2;
A(3,3)=P2;
A(5,4)=G3;
A(4,4)=P3;
A(5,5)=s4;
A(1,5)=f;
[lambda, v, w] = eigen(A);