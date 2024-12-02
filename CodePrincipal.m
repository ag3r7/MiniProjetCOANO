clear
clc

M1=400;
M2=40;
m=10;
e=0.1;
k1=400;
k2=1300;
w=4000*pi/30; %omega
g=9.81;

Z1_0=0;
Z2_0=0;
Zp1_0=0;
Zp2_0=0;
X0=[Z1_0;Z2_0;Zp1_0;Zp2_0];

t0=0;
tf=10;
dt=0.05*pi/w; %faire étude stab + dépendance de omega
t=t0:dt:tf;
N=length(t);

X=zeros(4,N);
X(:,1)=X0;

b=m*e*w^2*sin(w*t)/(M1+m);


A=[    0          0          1          0;
       0          0          0          1;
   -k1/(M1+m) k1/(M1+m)      0          0;
     k1/M2   -(k1+k2)/M2     0          0];

B0=[0;0;-g;-g];

Q_Eul_exp=(eye(4)+dt*A);
Q_Eul_imp=(eye(4)-dt*A)^-1;

Q_crank_1=((eye(4)-0.5*dt*A)^-1);
Q_crank_2=(eye(4)+0.5*dt*A);

Q_crank=Q_crank_1*Q_crank_2;

schema="EulerExp";

if (schema=="EulerExp")
    for i=1:N-1
        B=B0;
        B(3)=B(3)+b(i);
        X(:,i+1)=Q_Eul_exp*X(:,i)+dt*(B);
    end
end

hold on
plot(t,X(1,:),'r')
plot(t,X(2,:),'g')