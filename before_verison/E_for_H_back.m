clc
 clear all
 close all

N=100;
h=2*pi/(N-1);

Dif=zeros(N);
d=eye(N); % for the Dif
I=eye(N);
O=zeros(N);
value_EzHy=zeros(N,1);
value_EyHz=zeros(N,1);
tic;

%----------create difference matrix  backard-------
for i=1:N-1
   Dif(i*N+i)=1;
end
Dif_E=(1/h)*(Dif-d);
Dif_E(N)=(1/h)*1;

Dif_H=(-1)*transpose(Dif_E);


D=[-Dif_E O;O -Dif_H];
B=j*[O I;I O];  

%--------------eigenvalue ------------------
[vec_EzHy,lambda_EzHy]=eig(D,B);
[vec_EyHz,lambda_EyHz]=eig(-D,B);
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

val_want_EyHz=find(value_EyHz <2 & value_EyHz > -2);
val_want_EyHz_1=value_EyHz(value_EyHz <2 & value_EyHz > -2);
val_want_EzHy=find(value_EzHy <2 & value_EzHy > -2);
val_want_EzHy_1=value_EzHy(value_EzHy <2 & value_EzHy > -2);
toc;