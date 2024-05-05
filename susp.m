function [y2_mean,th2_mean,y2_sp,th2_sp] = susp(k1,k2,k3,k4,c1,c2,c3,c4)
%%
m1=1.93;%kg
m2=3.95;
m3=17.5;
p=450;%mm
q=450;
v=278;%10km/h=278mm/s로 뒀음.
h=55;%mm/휠베이스의 두께
J=m3*((p+q)^2+h^2)/12;%m3*((p+q)^2+h^2)/12,1693802.083

M=[m1 0 0 0;0 m2 0 0; 0 0 J 0;0 0 0 m3];
C=[c1+c2 0 c1*p -c1; 
   0 c3+c4 -c3*q -c3;
   c1*p -c3*q c1*p^2+c3*q^2 c3*q-c1*p;
   -c1 -c3 c3*q-p*c1 c1+c3];
K=[k1+k2 0 k1*p -k1; 
   0 k3+k4 -k3*q -k3;
   k1*p -k3*q k1*p^2+k3*q^2 k3*q-k1*p;
   -k1 -k3 k3*q-p*k1 k1+k3];

%modal analysis
K_modal=M^(-1/2)*K*M^(-1/2);
[Evector, Evalue]=eig(K_modal);
%Evalue의 양의 제곱근이 wi이다.
%modal force f(자갈길의 높이 20mm, 주기는 1초)
w=2*pi;
A=20;
%초기화
P=zeros(4,1);
F=zeros(4,1);
r=zeros(4,1);
ze=zeros(4,1);
X=zeros(4,1);    
wn=sqrt(Evalue);
y=zeros(2001,1);
th=zeros(2001,1);


for t=0:0.01:20%임의의 시간 20초 동안 측정
    for i=1:4
        wi=wn(i,1);
        ze(i,1)=c1/(2*k1)*wi;
        P(i,1)=atan(2*ze(i,1)*wi/(wi^2-w^2));
    end
    for j=1:4
        F=transpose(Evector)*(M^(-1/2))*[(c2)*A*w*cos(w*(t-(p+q)/v)-P(j,1))+(k2)*A*sin(w*(t-(p+q)/v)-P(j,1));
                                              (c4)*A*w*(cos(w*t-P(j,1)))+(k4)*A*(sin(w*t-P(j,1)));
                                               0;
                                               0];        

        r(j,1)=F(j,1)/(sqrt((wi^2-w^2)^2+(2*ze(j,1)*wi*w)^2));
        X=M^(-1/2)*Evector*r;
    end   
    y(round(100*t+1),1)=X(3,1);%y
    th(round(100*t+1),1)=X(4,1);%theta
end  

y2=diff(diff(y)/0.01)/0.01;%개수가 t 집합 개수보다 적다 (1999개)    
th2=diff(diff(th/0.01))/0.01;%개수가 t 집합 개수보다 적다 (1999개)
y2_mean=mean(y2);
th2_mean=mean(th2);
y2_sp=var(y2,1);
th2_sp=var(th2,1);