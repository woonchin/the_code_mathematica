clc
close all
clear all

%assume the interface_1 at N/2~N/2+1 and interface_2 at N/2+M ~N/2+M+1
N=300;
M=100;
h=2*pi/N-1;
kz=2;
alpha=1/3;%(x_head - x_i) /h assume put on the middle
beta=2/3; %(x_i+1-x_head)/h

epi_minus=1;
mu_minus=1;
gamma_minus=.9;

epi_plus=50;
mu_plus=50;
gamma_plus=0.2;

check_frequency=1;

Dif=zeros(N);
d=eye(N); % for the Dif
I=eye(N);
O_3n=zeros(3*N);
O=zeros(N);
eig_val=zeros(6*N,1);

tic; 

for s=1:N-1
   Dif(s*N+s)=1;
end
D_E=(1/h)*(Dif-d);
D_E(N)=(1/h)*1;
D_H=(-1)*transpose(D_E);

% %interface_1&2 
D_E(N/2,:)=0;
D_E(N/2+M,:)=0;
%for at N/2
D_E(N/2,N/2)=(1/h)*(-1/alpha);
D_E(N/2,N/2+1)=(1/h)*((1+beta)/alpha);
D_E(N/2,N/2+2)=(1/h)*(-beta/alpha);
%for at N/2+M
D_E(N/2+M,N/2+M)=(1/h)*(-1/alpha);
D_E(N/2+M,N/2+1+M)=(1/h)*((1+beta)/alpha);
D_E(N/2+M,N/2+2+M)=(1/h)*(-beta/alpha);
% %interface_1&2 
D_H(N/2+1,:)=0;
D_H(N/2+M+1,:)=0;
%for at N/2+1
D_H(N/2+1,N/2-1)=(1/h)*(alpha/beta); 
D_H(N/2+1,N/2)=(1/h)*((-alpha-1)/beta);
D_H(N/2+1,N/2+1)=(1/h)*(1/beta);
%for atN/2+M+1
D_H(N/2+1+M,N/2-1+M)=(1/h)*(alpha/beta); 
D_H(N/2+1+M,N/2+M)=(1/h)*((-alpha-1)/beta);
D_H(N/2+1+M,N/2+1+M)=(1/h)*(1/beta);

A_1=[O -i*kz*I O; i*kz*I O -D_E; O D_E O];
A_2=O_3n;
A_3=O_3n;
A_4=[O -i*kz*I O; i*kz*I O -D_H; O D_H O];
A=[A_1 A_2; A_3 A_4];

B_1=gamma_minus*I;
B_2=i*mu_minus*I;
B_3=-i*epi_minus*I;
B_4=gamma_minus*I;
for s=N/2+1:N/2+M
    B_1(s,s)=(gamma_plus/gamma_minus)*B_1(s,s);
    B_2(s,s)=(mu_plus/mu_minus)*B_2(s,s);
    B_3(s,s)=(epi_plus/epi_minus)*B_3(s,s);
    B_4(s,s)=(gamma_plus/gamma_minus)*B_4(s,s);
end
B_1=[B_1 O O;O B_1 O;O O B_1 ];
B_2=[B_2 O O;O B_2 O;O O B_2 ];
B_3=[B_3 O O;O B_3 O;O O B_3 ];
B_4=[B_4 O O;O B_4 O;O O B_4 ];
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
% 
% check_B1=zeros(3*N,1);
% check_B2=zeros(3*N,1);
% check_B3=zeros(3*N,1);
% check_B4=zeros(3*N,1);
% for s=1:3*N
%     check_B1(s)=B_1(s,s);
%     check_B2(s)=B_2(s,s);
%     check_B3(s)=B_3(s,s);
%     check_B4(s)=B_4(s,s);
% end
toc;