clc
close all
clear all

N=250;
h=2*pi/N-1;
kz=3;
epi=1;
mu=1;
gamma=.9;
check_frequency=2;

Dif=zeros(N);
d=eye(N); % for the Dif
I=eye(N);
O_3n=zeros(3*N);
I_3n=eye(3*N);
O=zeros(N);
eig_val=zeros(6*N,1);

tic; 

for s=1:N-1
   Dif(s*N+s)=1;
end
D_E=(1/h)*(Dif-d);
D_E(N)=(1/h)*1;
D_H=(-1)*transpose(D_E);

A_1=[O -i*kz*I O; i*kz*I O -D_E; O D_E O];
A_2=O_3n;
A_3=O_3n;
A_4=[O -i*kz*I O; i*kz*I O -D_H; O D_H O];
A=[A_1 A_2; A_3 A_4];

B_1=gamma*I_3n;
B_2=i*mu*I_3n;
B_3=-i*epi*I_3n;
B_4=gamma*I_3n;
B=[B_1 B_2; B_3 B_4];

%calculate eig
[vec,lambda]=eig(A,B);
for s=1:6*N
    eig_val(s)=lambda(s,s);
end
M1=mat2cell(vec,[N N N N N N]);
vec_Ex=M1{1,1};
vec_Ey=M1{2,1};
vec_Ez=M1{3,1};
vec_Hx=M1{4,1};
vec_Hy=M1{5,1};
vec_Hz=M1{6,1};

val_want=find(eig_val<check_frequency & eig_val>-check_frequency & eig_val~=0);
val_show=eig_val(eig_val<check_frequency & eig_val>-check_frequency & eig_val~=0);
toc;