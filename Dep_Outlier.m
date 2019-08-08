%% 加载数据，并归一化适合区域
close all;
clear all;
clc;
%% Path settings
inputImgPath = 'INPUT_IMG';                 % input image path
inputImgPath1 = 'INPUT_IMG1';                 % input depth_map path
inputImgPath2 = 'GLRNNR_saliency';                 % input 2D_saliency_map path
resSalPath = 'Depth_saliency';                     % Outlier_result path
if ~exist(resSalPath, 'file')
    mkdir(resSalPath);
end
% addpath
imgFiles = imdir(inputImgPath);
imgFiles1 = imdir(inputImgPath1);
imgFiles2 = imdir(inputImgPath2);
%%
load('camera_parms.mat');
load('name_cols110.mat');
name_T={Imgname11.name};
name_T=name_T';
for indImg = 1:113
    % read image
    tic
    imgPath = fullfile(inputImgPath, imgFiles(indImg).name);
    imgPath1 = fullfile(inputImgPath1, imgFiles1(indImg).name);
    imgPath2 = fullfile(inputImgPath2, imgFiles2(indImg).name);

    img.name = imgPath((strfind(imgPath,'\')+1):end);      
    A1 = imread(imgPath);
    D1 =double(imread(imgPath1));
    S1 =imread(imgPath2);

    
%% 加一步处理
SSS=D1;
STS=SSS(20:end,20:end);
ml=7;
    
SSS(abs(SSS-mean(SSS(:)))>ml*sqrt(var(STS(:))))=nan;
D1=SSS;
D1(isnan(SSS))=0;
%%
    
%% 加载相机内参数
load('camera_2019parms.mat');
load('name_cols110.mat');
name_T={Imgname11.name};
name_T=name_T';
t=find(ismember(name_T,imgFiles(indImg).name));

xx0=cols(t);
%% 图像预处理评估
AA=A1;
DD=(D1.*p1+p2).*~(~D1);
[m,n]=size(D1);
[x,y]=meshgrid(1:n,1:m);

%% 计算空间真实位置
xx=round(x+xx0-x0);
xx_C=xx./f.*DD;


%%

P=[];
yy=0.07*y;

%%%%%%%%%%%%%根据二维图计算初始显著区域位置%%%%%%%%%%%%%%%%%%%%%%

counts =hist(double(S1(:)),255);
[lv,em] = otsuthresh(counts);
Sal_init_DSR=im2bw(S1,0.5*lv);

[P,Ra,Circle_1]=Inscribed_Circle(Sal_init_DSR);
P=P(end:-1:1);

if 2.*Ra/m>0.7 | 2.*Ra/n>0.7 
    Ra=round(0.35*min(m,n));
elseif P(1)-Ra<5 | P(1)-Ra>n-5 | P(2)-Ra<5 | P(2)-Ra>m-5 
    Ra=Ra-2;
end
%%%%%%%%%%%%%%得到最大内接圆区域%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pos=[P(1),P(2)];

%%均匀化
Bw=0.*DD;
Bw((x-P(1)).^2+(y-P(2)).^2<=Ra.^2&(x-P(1)).^2+(y-P(2)).^2>=(Ra-2).^2)=1;

%% y方向运动精度0.07mm
x1=xx_C(Pos(2),Pos(1));
y1=0.07*y(Pos(2),Pos(1));
% yy=0.07*y;

sd=sqrt((xx_C-x1).^2+(yy-y1).^2);
U=Bw.*sd;
U=unique(U(:));
U(U==0)=[];
Rb=min(U);

Brw=0.*DD;
Brw((xx_C-x1).^2+(yy-y1).^2<Rb.^2)=1;

% 插值均匀性采样
es_x=Brw.*xx_C;
es_y=Brw.*yy;
es_z=Brw.*DD;

es = [es_x(:),es_y(:),es_z(:)];
es(es(:,1)==0&es(:,2)==0,:)=[];
%% 多次随机采样均值
num_Dian=size(es,1);
angle=zeros(10,10);

R_NJ=floor(num_Dian);
t=1;
xc=0;
yc=0;
t_num=3;%20;
while t<t_num
    [ex,ey]=DrawPoint(x1,y1,Rb,R_NJ);
    ex=ex';
    ey=ey';
    %%重投影到图像上进行插值
    es = [es_x(:),es_y(:),es_z(:)];
    es(es(:,1)==0&es(:,2)==0,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%
    ez=griddata(es(:,1),es(:,2),es(:,3),ex,ey);%插值
    % edata=[ex,ey,ez];
    %%%%%%%投影平面%%%%%%%%%%%%%%%
    [fitresult, gof] = createFit_exyz(ex, ey, ez);
    %% 找出对应投影点位置： z=f(x,y) = p00 + p10*x + p01*y ; p10*x+p01*y-z+p00=0;
    a=fitresult.p10;
    b=fitresult.p01;
    c=-1;
    d=fitresult.p00;
    ed=abs(a.*ex+b.*ey+c.*ez+d)./sqrt(a.^2+b.^2+c.^2);
    ez=ed;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    ex(isnan(ez))=[];
    ey(isnan(ez))=[];
    ez(isnan(ez))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 先根据圆内数据拟合出方向，es

es=[ex,ey,ez];
% v=[xcc-x1,ycc-y1];
% v=v./norm(v);
tt=1;
sum_Loss=[];
for tn=-3:0.05:3
v=[cos(tn*pi/180),sin(tn*pi/180)];
es2=es(:,1:2)-[x1,y1];
es2=es2.*v;
es2=sum(es2');
es2=es2';

%% 缩减点间距

es3=0.33*n.*(es2-min(es2))./(max(es2)-min(es2));
es3=round(es3);

idx=unique(es3);
num_idx=length(idx);

for i=1:num_idx
    L_xy=find(es3==idx(i));
    temp=es(L_xy,:);
%%%%%%%%%%%%%%使用RANSAC 算法拟合直线%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Loss_N(i)=var(temp(:,3));
end
    sum_Loss(tt)=sum(Loss_N);
    tt=tt+1;
end
%%%%%%%%%%%%%%%%%%%%确定方向%%%%%%%%%%%%%%%%%%%%%%%%%

sum_Loss1=medfilt1(sum_Loss,3);
SV_1=mean(find(sum_Loss1==min(sum_Loss1)));

theta=(SV_1-1)*0.05-3;

%% 所有的值
tn=theta;
%%

%% 更改后
ess1=[xx_C(:),yy(:),DD(:),x(:),y(:)];

es=[xx_C(:),yy(:),double(DD(:)),x(:),y(:)];
essD=[xx_C(:),yy(:),double(DD(:)),x(:),y(:)];
%%
es(ess1(:,3)==0,:)=[];
essD(ess1(:,3)==0,:)=[];

%%
es1=es(:,4:5);
es=es(:,1:3);

ex=es(:,1);
ey=es(:,2);
%%
v=[cos(tn*pi/180),sin(tn*pi/180)];
es2=es(:,1:2)-[mean(ex(:)),mean(ey(:))];
es2=es2.*v;
es2=sum(es2');
es2=es2';

%%
es3=0.33*n.*(es2-min(es2))./(max(es2)-min(es2));
es3=round(es3);

idx=unique(es3);
num_idx=length(idx);

for i=1:num_idx
        L_xy=find(es3==idx(i));
        temp=es(L_xy,:);
        tempD=essD(L_xy,:);
        tempD(tempD(:,3)==0,:)=[];
        if length(tempD)>8
            d_thresh=Ransuc_L(tempD,temp);
            es(L_xy,3)=d_thresh;
        else 
            
%%%%%%%%%%%%%%使用RANSAC 算法拟合直线%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if var(temp(:,3))<0.6       
            Loss_N=median(temp(:,3));
            es(L_xy,4)=abs(es(L_xy,3)-Loss_N);    
%         else 
            tempp=temp(:,3);
            d1=mean(temp(:,3));
            dv1=var(tempp(tempp<d1));
            dv2=var(tempp(tempp>=d1));
            if dv1<=dv2
                d21=tempp(tempp<d1);
            else 
                d21=tempp(tempp>=d1);
            end
                d22=median(d21);
                dv21=var(d21(d21<d22));
                dv22=var(d21(d21>=d22));
            if dv21<=dv22
                d_thresh=median(d21(d21<d22));
            else 
                d_thresh=median(d21(d21>=d22));
            end   
            es(L_xy,3)=abs(es(L_xy,3)-d_thresh);
            
        end
        end
            %% 数据分割后，找到 var小的 中位数        
        end
% end
Sal=0.*DD;
dat=es(:,3);%+0*es(:,4);

%%
for ii=1:length(es(:,3))
Sal(es1(ii,2),es1(ii,1))=dat(ii);
end
Sal=Sal-min(dat(:));
Sal(Sal<0)=0;
Sal_L{t}=Sal;
angle(t)=theta;
    t=t+1;
end

aa=abs(DD-medfilt2(DD,[7,7]));
Mask=abs(DD-medfilt2(DD,[7,7]))>var(aa(:));
se1=strel('disk',10);%这里是创建一个半径为5的平坦型圆盘结构元素
Mask2=imdilate(Mask,se1);
se1=strel('disk',15);%这里是创建一个半径为5的平坦型圆盘结构元素
Mask3=imdilate(~DD,se1);

%% 
Sal_DL=Sal_L{1}+Sal_L{2};
Sal_DL=~Mask3.*~Mask2.*Sal_DL.*im2bw(DD);
%%
DSL=Sal_DL(20:end-20,20:end-20);
Sal_DL=DSL;
Sal_DLx=Sal_DL.*x(20:end-20,20:end-20);
Sal_DLy=Sal_DL.*y(20:end-20,20:end-20);
Sal_DLxy=sum(Sal_DL(:));
Sx=sum(Sal_DLx(:))/Sal_DLxy;
Sy=sum(Sal_DLy(:))/Sal_DLxy;
sigxy=min(var(x(:)),var(y(:)));
Sal_DL=Sal_L{1}+Sal_L{2};
%%
Sal_DL=~Mask3.*~Mask2.*Sal_DL.*im2bw(DD);
%%
Sal_DL1=Sal_DL.*(1+exp(-((x-Sx).^2+(y-Sy).^2)./sigxy.^2));
Sal_DL=Sal_DL1;
 DSL=Sal_DL(20:end-20,20:end-20);
Sal_DL=(Sal_DL-min(DSL(:)))./range(DSL(:));
salPath = fullfile(resSalPath, strcat(img.name(1:end-4), '.png'));  
imwrite(Sal_DL,salPath);
toc
disp(toc)
end