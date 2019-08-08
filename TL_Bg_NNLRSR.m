% ||L||*+||Z||*+a*||H1||+b*||S1||+c*||H2||+d*||S2||+e*||S3||
% s.t. F=L+S1,L=H1; L=A*Z+S2;Z=H2,H2>=0
%%加入边界惩罚项
%s.t. S3=F(:,bg)-A; namely: S3=S1(:,bg);
function [L,S1,A,Z2,S2]=TL_Bg_NNLRSR(featMat,bg,a)
    addpath(genpath('./inexact_alm_rpca'));
    [LL,SS1,~] = inexact_alm_rpca(featMat, a(2));
    bgfeat=LL(bg,:);
    feat=LL;
    L=feat';
    S1=SS1';
%     A=bgfeat';
    A=(featMat(bg,:))';
    F=featMat';
    
    a1=a(1);
    a2=a(2);
    a3=a(3);
    a4=a(4);
    

    [m,n]=size(L);
    [mm,nn]=size(A);

% H1=zeros(nn,n);
    H1=L;

    Y_1=0.*L;
    Y_2=0.*L;


%     Z2=zeros(nn,n);
%     H2=zeros(nn,n);
%     S2=0.*L;
    [Z2,H2,S2]=LADMAP(L,A,a3,a4);
    
    Y_3=0.*L;
    Y_4=zeros(nn,n);
    %%
    u_1=0.1;
    u_1max=10^10;
    p1_0=1.1;
    chy11=1e-5;
%     chy12=1e-2;
    chy12=0.015;

    %%
    u_2=0.1;
    u_2max=10^10;
    p2_0=1.1;
    chy21=1e-5;
%     chy22=1e-2;
    chy22=0.015;

    eta1=(norm(A,2)).^2;
    %%

    term11=1;
    term12=0;
    term21=1;
    term22=0;
    t=0;
    tt1=1e6;

while ( (term11>=chy11) | (term12>=chy12) |(term21>=chy21) | (term22>=chy22) ) & ((term11>=chy11) |t<200)

T1=F-S1+Y_1./u_1;
T2=A*Z2+S2-Y_3./u_2;
T3=H1-Y_2./u_1;
%% L
L_kn1=L;
thresh_L=1./(2*u_1+u_2);
L=(u_1*(T1+T3)+u_2*T2)./(2*u_1+u_2);

[U_L,sigma_L,V_L]=svd(L);
mn_L=max(size(sigma_L));
thresh_LD=diag(thresh_L.*ones(mn_L,1));
thresh_LD=thresh_LD(1:size(sigma_L,1),1:1:size(sigma_L,2));
T_L=sigma_L-thresh_LD;
T_L1=max(T_L,zeros(size(sigma_L)));
L=U_L*T_L1*V_L';

%% H1
H_kn1=H1;
thresh_H1=a1/u_1;
H1=L+Y_2/u_1;

[U_H1,sigma_H1,V_H1]=svd(H1);
mn_H1=max(size(sigma_H1));
thresh_H1D=diag(thresh_H1.*ones(mn_H1,1));
thresh_H1D=thresh_H1D(1:size(sigma_H1,1),1:1:size(sigma_H1,2));
T_H1=sigma_H1-thresh_H1D;
T_H11=max(T_H1,zeros(size(sigma_H1)));
H1=U_H1*T_H11*V_H1';


%% S1
S_kn1=S1;
thresh_S1=a2/u_1;
S1=F-L+Y_1/u_1;

[U_S1,sigma_S1,V_S1]=svd(S1);
mn_S1=max(size(sigma_S1));
thresh_S1D=diag(thresh_S1.*ones(mn_S1,1));
thresh_S1D=thresh_S1D(1:size(sigma_S1,1),1:1:size(sigma_S1,2));
T_S1=sigma_S1-thresh_S1D;
T_S11=max(T_S1,zeros(size(sigma_S1)));
S1=U_S1*T_S11*V_S1';


%% Z2
Z_kn2=Z2;
% A=L(:,bg);
eta1=norm(A,2).^2;
thresh_Z2=(eta1*u_2).^(-1);
Z2=Z2+((A'*(F-A*Z2-S2+Y_3/u_2))-(Z2-H2+Y_4/u_2))/eta1;

[U_Z2,sigma_Z2,V_Z2]=svd(Z2);
mn_Z2=max(size(sigma_Z2));
thresh_Z2D=diag(thresh_Z2.*ones(mn_Z2,1));
thresh_Z2D=thresh_Z2D(1:size(sigma_Z2,1),1:1:size(sigma_Z2,2));
T_Z2=sigma_Z2-thresh_Z2D;
T_Z21=max(T_Z2,zeros(size(sigma_Z2)));
Z2=U_Z2*T_Z21*V_Z2';


%% H2
H_kn2=H2;
thresh_H2=a3/u_2;
H2=Z2+Y_4/u_2;

[U_H2,sigma_H2,V_H2]=svd(H2);
mn_H2=max(size(sigma_H2));
thresh_H2D=diag(thresh_H2.*ones(mn_H2,1));
thresh_H2D=thresh_H2D(1:size(sigma_H2,1),1:1:size(sigma_H2,2));
T_H2=sigma_H2-thresh_H2D;
T_H21=max(T_H2,zeros(size(sigma_H2)));
H2=U_H2*T_H21*V_H2';
H2=max(H2,0.*H2);

%% S2
S_kn2=S2;
thresh_S2=a4/u_2;
S2=L-A*Z2+Y_3/u_2;

[U_S2,sigma_S2,V_S2]=svd(S2);
mn_S2=max(size(sigma_S2));
thresh_S2D=diag(thresh_S2.*ones(mn_S2,1));
thresh_S2D=thresh_S2D(1:size(sigma_S2,1),1:1:size(sigma_S2,2));
T_S2=sigma_S2-thresh_S2D;
T_S21=max(T_S2,zeros(size(sigma_S2)));
S2=U_S2*T_S21*V_S2';

Y_1=Y_1+u_1*(F-L-S1);
Y_2=Y_2+u_1*(L-H1);
Y_3=Y_3+u_2*(L-A*Z2-S2);
Y_4=Y_4+u_2*(Z2-H2);


if term12<chy12
    p1=p1_0;
else
    p1=1;
end
u_1=min(u_1max,p1*u_1);
term11=norm(F-L-S1,'fro')/norm(F,'fro');
term12=u_1*max(max(norm(L-L_kn1,'fro'),norm(H1-H_kn1,'fro')),norm(S1-S_kn1,'fro'))/norm(F,'fro');

if term22<chy22
    p2=p2_0;
else
    p2=1;
end
u_2=min(u_2max,p2*u_2);
term21=norm(L-A*Z2-S2,'fro')/norm(L,'fro');
term22=u_2*max(max(sqrt(eta1)*norm(Z2-Z_kn2,'fro'),norm(H2-H_kn2,'fro')),norm(S2-S_kn2,'fro'))/norm(L,'fro');
%%

% pause
% disp([term11,term12,term21,term22]);

tt=sum([term12,term21,term22]);
if tt<=tt1
    L1=L;
    S11=S1;
    A1=A;
    Z12=Z2;
    S12=S2;
%     L,S1,A,Z2,S2
    tt1=tt;
end
t=t+1;
end

%      time2=toc;
%     disp(time2);
    L=L1;
    S1=S11;
    A=A1;
    Z2=Z12;
    S2=S12;



