clc
close all
clear all
Num=[100 200 300 400 500 600 700 800 900 1000];
cnt=1;
criterion_time=zeros(length(Num),1);
kz=2;
kx=2;
alpha=1/2;%(x_head - x_i) /h assume put on the middle
beta=1/2; %(x_i+1-x_head)/h
 epi_minus=1;
 mu_minus=1;
 gamma_minus=0;

 epi_plus=1;
 mu_plus=10;
 gamma_plus=.5;
 check_frequency=2;
 %%
while cnt<length(Num)+1
    
    %assume the interface_1 at N/2~N/2+1 and interface_2 at N/2+M ~N/2+M+1
    tic;
    N=Num(cnt);
    M=N/5;
    h=2*pi/(N-1);
    Dif=zeros(N);
    d=eye(N); % for the Dif
    I=eye(N);
    O_3n=zeros(3*N);
    O=zeros(N);
    eig_val=zeros(6*N,1);
    %%  difference matrix A with Ghost state
        for s=1:N-1
           Dif(s*N+s)=1;
        end
    D_E=(1/h)*(Dif-d);
    D_E(N)=(1/h)*1;
    D_H=(-1)*transpose(D_E);

    % interface_1&2 
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
     %interface_1&2 
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
    A_4=[O -1i*kz*I O; i*kz*I O -D_H; O D_H O];
    A=[A_1 A_2; A_3 A_4];
    %% B
    B_1=gamma_minus*I;
    B_2=i*mu_minus*I;
    B_3=-i*epi_minus*I;
    B_4=gamma_minus*I;
        for s=N/2+1:N/2+M
            B_1(s,s)=gamma_plus*I(s,s);%fix for gamma=0
            B_2(s,s)=i*mu_plus*I(s,s);
            B_3(s,s)=-i*epi_plus*I(s,s);
            B_4(s,s)=gamma_plus*I(s,s);
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
     time=toc;
    % %% Null Basis
    % ccnt=1;
    % det_num=zeros(length(val_want),1);
    %     for q=1:length(val_want) 
    %         [M,Det]=NUA(gamma_plus,epi_plus,mu_plus,gamma_minus,epi_minus,mu_minus,kx,kz,val_show(q));
    %          det_num(q)=Det;
    %         if  abs(Det)< 1e-10
    %             basis{q}=null(M);            
    %             ccnt=ccnt+1;
    %         end
    %     end
    %% check with jump condition
    cri_left_1=zeros(length(val_want),1);
    cri_left_2=zeros(length(val_want),1);
    cri_right_1=zeros(length(val_want),1);
    cri_right_2=zeros(length(val_want),1);
        for t=1:length(val_want)
        Hx1_plus=vec_Hx(N/2+1,val_want(t));
        Hx1_minus=vec_Hx(N/2,val_want(t));
        Ex1_plus=vec_Ex(N/2+1,val_want(t));
        Ex1_minus=vec_Ex(N/2,val_want(t));
        
        Hx2_plus=vec_Hx(N/2+M,val_want(t));
        Hx2_minus=vec_Hx(N/2+M+1,val_want(t));
        Ex2_plus=vec_Ex(N/2+M,val_want(t));
        Ex2_minus=vec_Ex(N/2+M+1,val_want(t));
        
        %[mu Hx]-i*[r Ex]=0
        %[epi Ex]+i*[r Hx]=0
        cri_left_1(t)=i*val_show(t)*(( mu_plus*Hx1_plus-mu_minus*Hx1_minus)-i*( gamma_plus*Ex1_plus-gamma_minus*Ex1_minus ));
        cri_left_2(t)=i*val_show(t)*(( epi_plus*Ex1_plus-epi_minus*Ex1_minus )+i*( gamma_plus*Hx1_plus-gamma_minus*Hx1_minus ));        
        cri_right_1(t)=i*val_show(t)*(( mu_minus*Hx2_minus-mu_plus*Hx2_plus)-i*( gamma_minus*Ex2_minus-gamma_plus*Ex2_plus ));
        cri_right_2(t)=i*val_show(t)*(( epi_minus*Ex2_minus-epi_plus*Ex2_plus )+i*( gamma_minus*Hx2_minus-gamma_plus*Hx2_plus ));    
        end
        criterion_left_err1{1,cnt}=cri_left_1;
        criterion_left_err2{1,cnt}=cri_left_2;
        criterion_right_err1{1,cnt}=cri_right_1;
        criterion_right_err2{1,cnt}=cri_right_2;
       
        criterion_time(cnt,1)=time;
        cnt=cnt+1;
end