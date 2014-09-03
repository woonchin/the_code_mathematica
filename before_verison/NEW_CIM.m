clc
close all
clear all

N=2500;
M=200;%medium length
epi_minus=1;
epi_plus=1+.002;
mu_minus=1;
mu_plus=1+.002;
alpha=1/2; %put interface at middle
beta=1/2;

h=2*pi/(N-1);
Dx=zeros(N);
O=zeros(N/2);
I=eye(N);
I_half=eye(N/2);% interface on the half
value_Ez=zeros(N,1);

tic;

rho_plus_1=epi_plus/(alpha*epi_plus+beta*epi_minus);
rho_minus_1=epi_minus/(alpha*epi_plus+beta*epi_minus);
rho_plus_2=epi_minus/(alpha*epi_minus+beta*epi_plus);
rho_minus_2=epi_plus/(alpha*epi_minus+beta*epi_plus);
for i=1:N-1
   Dx(i*N+i)=1;
end
Dx=Dx+transpose(Dx);
Dx=Dx-2*I;
Dx(N*(N-1)+1)=1;
Dx(N)=1;
Dx= (-1/epi_minus)*(1/h^2)*Dx;
%---assume interface_1 on the N/2~N/2+1
Dx(N/2,:)=0;
Dx(N/2+1,:)=0;
% on the grid_N/2
Dx(N/2,N/2-1)= -epi_minus/h^2;
Dx(N/2,N/2)= epi_minus*(1+rho_plus_1)/h^2;
Dx(N/2,N/2+1)= -epi_minus*rho_plus_1/h^2;
%on the grid_N/2+1
Dx(N/2+1,N/2)= -epi_plus*rho_minus_1/h^2;
Dx(N/2+1,N/2+1)= epi_plus*(1+rho_minus_1)/h^2;
Dx(N/2+1,N/2+2)= -epi_plus/h^2;
%---assume interface_2 on the N/2+M~N/2+M+1
for v=N/+2:N/2+M-1 % medium permitivity
    Dx(v,:)= (epi_minus/epi_plus)*Dx(v,:);
end
 Dx(N/2+M,:)=0;
 Dx(N/2+M+1,:)=0;
 %on the grid N/2+M
 Dx(N/2+M,N/2+M-1)= -epi_plus/h^2;
 Dx(N/2+M,N/2+M)= epi_plus*(1+rho_plus_2)/h^2;
 Dx(N/2+M,N/2+M+1)= -epi_plus*rho_plus_2/h^2;
 %on the grid N/2+M+1
 Dx(N/2+M+1,N/2+M)= -epi_minus*rho_minus_2*h^2;
 Dx(N/2+M+1,N/2+M+1)= epi_minus*(1+rho_minus_2)/h^2;
 Dx(N/2+M+1,N/2+M+2)= -epi_minus/h^2;
 I_f=[mu_minus*I_half O; O mu_minus*I_half];
for  t=N/2+1:N/2+M
    I_f(t,t)=(mu_plus/mu_minus)*I_f(t,t);
end
%--eigenvalue
[vec_Ez,lambda_Ez]=eig(Dx,I_f);

for s=1:N
    value_Ez(s)=sqrt(lambda_Ez(s,s));
end
val_want_Ez=find(value_Ez <2 & value_Ez > -2);
value_Ez_show=value_Ez(value_Ez <2 & value_Ez > -2);
toc;