clc
close all
clear all

%assume the interface_1 at N/2~N/2+1 and interface_2 at N/2+M ~N/2+M+1
N=500;
M=N/10;
h=2*pi/N-1;
kz=2;

epi_minus=4;
mu_minus=4;
gamma_minus=2;

epi_plus=4;
mu_plus=4;
gamma_plus=2;

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

% %interface_1&2 
D_E(N/2,:)=0;
D_E(N/2+M,:)=0;
%for at N/2
D_E(N/2,N/2)=(1/h)*1;
D_E(N/2,N/2-1)=(1/h)*-1;
%for at N/2+M
D_E(N/2+M,N/2+M)=(1/h)*1;
D_E(N/2+M,N/2+M-1)=(1/h)*-1;
% %interface_1&2 
D_H(N/2+1,:)=0;
D_H(N/2+M+1,:)=0;
%for at N/2+1
D_H(N/2+1,N/2+1)=-1;
D_H(N/2+1,N/2+2)=1;
%for atN/2+M+1
D_H(N/2+M+1,N/2+M+1)=-1;
D_H(N/2+M+1,N/2+M+2)=1;

A_1=[O -i*kz*I O; i*kz*I O -D_E; O D_E O];
A_2=O_3n;
A_3=O_3n;
A_4=[O -i*kz*I O; i*kz*I O -D_H; O D_H O];
A=[A_1 A_2; A_3 A_4];

B_1=gamma_minus*I_3n;
B_2=i*mu_minus*I_3n;
B_3=-i*epi_minus*I_3n;
B_4=gamma_minus*I_3n;
for s=N/2+1:N/2+M
    B_1(s,s)=(gamma_plus/gamma_minus)*B_1(s,s);
    B_1(s+N,s+N)=(gamma_plus/gamma_minus)*B_1(s+N,s+N);
    B_1(s+2*N,s+2*N)=(gamma_plus/gamma_minus)*B_1(s+2*N,s+2*N);
    
    B_2(s,s)=(mu_plus/mu_minus)*B_2(s,s);
    B_2(s+N,s+N)=(mu_plus/mu_minus)*B_2(s+N,s+N);
    B_2(s+2*N,s+2*N)=(mu_plus/mu_minus)*B_2(s+2*N,s+2*N);
    
    B_3(s,s)=(epi_plus/epi_minus)*B_3(s,s);
    B_3(s+N,s+N)=(epi_plus/epi_minus)*B_3(s+N,s+N);
    B_3(s+2*N,s+2*N)=(epi_plus/epi_minus)*B_3(s+2*N,s+2*N);
    
    B_4(s,s)=(gamma_plus/gamma_minus)*B_4(s,s);
    B_4(s+N,s+N)=(gamma_plus/gamma_minus)*B_4(s+N,s+N);
    B_4(s+2*N,s+2*N)=(gamma_plus/gamma_minus)*B_4(s+2*N,s+2*N);
end
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