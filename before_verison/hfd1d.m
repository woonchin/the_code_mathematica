function  k=hfd1d(length,N)

h=length/N;
A=spalloc(N-1,...
          N-1,...
          3*(N-1));
d=-2/h^2;
s=1/h^2;

for i=1:N-1
    A(i,i)=d;
end
for i=1:N-2
    A(i,i+1)=s;
    A(i+1,i)=s;
end
lamda=eig(A);
k=sqrt(sort(lamda));