clear all
close all

N=15;
h=100;
i=1;
k_z=10;
D=zeros(N);
I=eye(N);
O=zeros(N);
O_B=zeros(3*N);
for i=1:N-1
    D(16+(i-1)*(N+1))=1/(2*h);
    D(2+(i-1)*(N+1))=-1/(2*h);  
end
a=[O -i*k_z*I O;-i*k_z*I O D; O D O];
A=[a O_B;O_B a];