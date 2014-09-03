clc
close all 
clear all

N=1000;
h=2*pi/(N-1);

Dif=zeros(N);
i=sqrt(-1);
O=zeros(N);
I=eye(N);
value_EzHy=zeros(N,1);
value_EyHz=zeros(N,1);

tic;
%------create matrix-----
%--------d/dx E(x)-----
for i=1:N-1
        Dif(i*N+i)=1;
end
Dif_E=(1/h)*(Dif-I);
Dif_E(2*N)=(1/h)*1;
%--------d/dx H(x)-----
Dif_n=transpose(Dif);
Dif_H=(1/(h))*(Dif+(-1)*Dif_n);
Dif_H((N-2)*N+1)=-1/(h);
Dif_H(2*N)=1/(h);

D=[-Dif_E O;O -Dif_H];
B=i*[O I;I O];  
[vec_EzHy,lambda_EzHy]=eig(D,B);
[vec_EyHz,lambda_EyHz]=eig(-D,B);
%--------------eigenvalue ------------------
for s=1:2*N
    value_EzHy(s)=lambda_EzHy(s,s);
    value_EyHz(s)=lambda_EyHz(s,s);
end
M1=mat2cell(vec_EzHy,[N N]);
M2=mat2cell(vec_EyHz,[N N]);
vec_Ez=M1{1,1};
vec_Hy=M1{2,1};
vec_Ey=M2{1,1};
vec_Hz=M2{2,1};

val_want_EyHz=find(value_EyHz <1 & value_EyHz > -1);
val_want_EzHy=find(value_EzHy <1 & value_EzHy > -1);
toc;