clc
close all
N=1000;
Y=[991 1000];
t=1;
% Ez Hy
while t<3
    
    for i=1:N
        A1(i)=real(vec_Ez(:,Y(t));
        A2(i)=real(vec_Ez(:,Y(t)+1));
        A3(i)=real(vec_Ez(:,Y(t)+2));
        A4(i)=real(vec_Ez(:,Y(t)+3));
        A5(i)=real(vec_Ez(:,Y(t)+4));
        A6(i)=real(vec_Ez(:,Y(t)+5));
        A7(i)=real(vec_Ez(:,Y(t)+6));
        A8(i)=real(vec_Ez(:,Y(t)+7));
        A9(i)=real(vec_Ez(:,Y(t)+8));
        A10(i)=real(vec_Ez(:,Y(t)+9));
    end
    subplot(10,1,1);
    plot(A1);
    title('Ey');
    subplot(10,1,2);
    plot(A2);
    subplot(10,1,3);
    plot(A3);
    subplot(10,1,4);
    plot(A4);
    subplot(10,1,5);
    plot(A5);
    subplot(10,1,6);
    plot(A6);
    subplot(10,1,7);
    plot(A7);
    subplot(10,1,8);
    plot(A8);
    subplot(10,1,9);
    plot(A9);
    subplot(10,1,10);
    plot(A10);
    figure;

    for i=1:N
        B1(i)=imag(vec_Hy(:,Y(t)));
        B2(i)=imag(vec_Hy(:,Y(t)+1));
        B3(i)=imag(vec_Hy(:,Y(t)+2));
        B4(i)=imag(vec_Hy(:,Y(t)+3));
        B5(i)=imag(vec_Hy(:,Y(t)+4));
        B6(i)=imag(vec_Hy(:,Y(t)+5));
        B7(i)=imag(vec_Hy(:,Y(t)+6));
        B8(i)=imag(vec_Hy(:,Y(t)+7));
        B9(i)=imag(vec_Hy(:,Y(t)+8));
        B10(i)=imag(vec_Hy(:,Y(t)+9));
        
    end
    subplot(10,1,1);
    plot(B1);
    title('Hz');
    subplot(10,1,2);
    plot(B2);
    subplot(10,1,3);
    plot(B3);
    subplot(10,1,4);
    plot(B4);
    subplot(10,1,5);
    plot(B5);
     subplot(10,1,6);
    plot(B6);
     subplot(10,1,7);
    plot(B7);
    subplot(10,1,8);
    plot(B8);
     subplot(10,1,9);
    plot(B9);
     subplot(10,1,10);
    plot(B10);
    figure;
    t=t+1;
end 

q=1;
%Ey Hz
while q<3
    for i=1:N
        A1(i)=real(vec_Ey(:,Y(q)));
        A2(i)=real(vec_Ey(:,Y(q)+1));
        A3(i)=real(vec_Ey(:,Y(q)+2));
        A4(i)=real(vec_Ey(:,Y(q)+3));
        A5(i)=real(vec_Ey(:,Y(q)+4));
        A6(i)=real(vec_Ey(:,Y(q)+5));
        A7(i)=real(vec_Ey(:,Y(q)+6));
        A8(i)=real(vec_Ey(:,Y(q)+7));
        A9(i)=real(vec_Ey(:,Y(q)+8));
        A10(i)=real(vec_Ey(:,Y(q)+9));
    end
    
    subplot(10,1,1);
    plot(A1);
    title('Ey');
    subplot(10,1,2);
    plot(A2);
    subplot(10,1,3);
    plot(A3);
    subplot(10,1,4);
    plot(A4);
    subplot(10,1,5);
    plot(A5);
    subplot(10,1,6);
    plot(A6);
    subplot(10,1,7);
    plot(A7);
    subplot(10,1,8);
    plot(A8);
    subplot(10,1,9);
    plot(A9);
    subplot(10,1,10);
    plot(A10);
    figure;
    
    for i=1:N
        B1(i)=imag(vec_Hy(:,Y(q)));
        B2(i)=imag(vec_Hy(:,Y(q)+1));
        B3(i)=imag(vec_Hy(:,Y(q)+2));
        B4(i)=imag(vec_Hy(:,Y(q)+3));
        B5(i)=imag(vec_Hy(:,Y(q)+4));
        B6(i)=imag(vec_Hy(:,Y(q)+5));
        B7(i)=imag(vec_Hy(:,Y(qt)+6));
        B8(i)=imag(vec_Hy(:,Y(q)+7));
        B9(i)=imag(vec_Hy(:,Y(q)+8));
        B10(i)=imag(vec_Hy(:,Y(q)+9));
        
    end
    subplot(10,1,1);
    plot(B1);
    title('Hz');
    subplot(10,1,2);
    plot(B2);
    subplot(10,1,3);
    plot(B3);
    subplot(10,1,4);
    plot(B4);
    subplot(10,1,5);
    plot(B5);
     subplot(10,1,6);
    plot(B6);
     subplot(10,1,7);
    plot(B7);
    subplot(10,1,8);
    plot(B8);
     subplot(10,1,9);
    plot(B9);
     subplot(10,1,10);
    plot(B10);
    figure;
    q=q+1;
end 