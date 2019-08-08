function [Z,H,E]=LADMAP(X,A,beta,lambda)
% Z0=H0=E0=Y10=Y20=0;
[m,n]=size(X);
[mm,nn]=size(A);
Z_k1=zeros(nn,n);
H_k1=zeros(nn,n);
E_k1=0.*X;
Y_1k=0.*X;
Y_2k=zeros(nn,n);
u_k=0.1;
u_max=10^10;
p0=1.1;
chy1=1e-6;
chy2=1e-2;
eta1=(norm(A,2)).^2;
% k=0;

term1=1;
term2=0;
t=0;

while ( (term1>=chy1) | (term2>=chy2) )
Z_kn1=Z_k1;
Z_k1=Z_k1+((A'*(X-A*Z_k1-E_k1+Y_1k/u_k))-(Z_k1-H_k1+Y_2k/u_k))/eta1;
thresh_Z=(eta1*u_k).^(-1);
% [U_z,sigma_Z,V_z]=svd(Z_k1,'econ');
[U_z,sigma_Z,V_z]=svd(Z_k1);
mn_z=max(size(sigma_Z));
thresh_ZD=diag(thresh_Z.*ones(mn_z,1));
thresh_ZD=thresh_ZD(1:size(sigma_Z,1),1:1:size(sigma_Z,2));
T_Z=sigma_Z-thresh_ZD;
T_Z1=max(T_Z,zeros(size(sigma_Z)));
Z_k1=U_z*T_Z1*V_z';

H_kn1=H_k1;
H_k1=Z_k1+Y_2k/u_k;
thresh_H=beta./u_k;
% [U_h,sigma_H,V_h]=svd(H_k1,'econ');
[U_h,sigma_H,V_h]=svd(H_k1);
mn_h=max(size(sigma_H));
thresh_HD=diag(thresh_H.*ones(mn_h,1));
thresh_HD=thresh_HD(1:size(sigma_H,1),1:1:size(sigma_H,2));
T_H=sigma_H-thresh_HD;
T_H1=max(T_H,zeros(size(sigma_H)));
H_k1=U_h*T_H1*V_h';
H_k1=max(H_k1,0.*H_k1);

E_kn1=E_k1;
E_k1=X-A*Z_k1+Y_1k/u_k;
thresh_E=lambda./u_k;
% [U_e,sigma_E,V_e]=svd(E_k1,'econ');
[U_e,sigma_E,V_e]=svd(E_k1);
mn_e=max(size(sigma_E));
thresh_ED=diag(thresh_E.*ones(mn_e,1));
thresh_ED=thresh_ED(1:size(sigma_E,1),1:1:size(sigma_E,2));
T_E=sigma_E-thresh_ED;
T_E1=max(T_E,zeros(size(sigma_E)));
E_k1=U_e*T_E1*V_e';
% E_k1=u_k*(E_k1)./(2*lambda+u_k);


Y_1k=Y_1k+u_k*(X-A*Z_k1-E_k1);
Y_2k=Y_2k+u_k*(Z_k1-H_k1);

if term2<chy2
    p=p0;
else
    p=1;
end
u_k=min(u_max,p*u_k);
term1=norm(X-A*Z_k1-E_k1,'fro')/norm(X,'fro');
term2=u_k*max(max(sqrt(eta1)*norm(Z_k1-Z_kn1,'fro'),norm(H_k1-H_kn1,'fro')),norm(E_k1-E_kn1,'fro'))/norm(X,'fro');
t=t+1;
% disp(t)
end

Z=Z_k1;
H=H_k1;
E=E_k1;
