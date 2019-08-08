close all;
clear;
clc;
tic;
% Gt_path='.\MSRA-B\';
Gt_path='.\my_GT1210\';
% Salcut_path='.\CA\';
Salcut_path='.\NLNL_SPR200\';
Salcut_path1='.\NLNL_SPR300\';
Salcut_path2='.\NLNL_SPR350\';
Salcut_path3='.\NLNL_SPR400\';
Salcut_path4='.\NLNL_SPR700\';
Gt_list=dir(strcat(Gt_path,'*.png'));             % 单通道

% load('NLNL_SPR77.mat')
% a=optimresults.x;
% a=ones(5,1);
a=[1,0,0,1,0];

Salcut_list=Gt_list;
for i=1:length(Gt_list)

Salcut_list(i).name=[Salcut_list(i).name(1:end-4),'.png'];
% Salcut_list1(i).name=[Salcut_list(i).name(1:end-4),'.png'];
% Salcut_list2(i).name=[Salcut_list(i).name(1:end-4),'.png'];
end
imgNum=length(Salcut_list);                       % 图片数量

if imgNum~=length(Gt_list)
    error('图片数量不等');
end
%--------------------------------------------------%
%---------------------计算公式—————————-—-%
%  精确率P=(M∩G)/M          召回率R=(M∩G)/G       %
%     P=TP/(TP+FP)               R=TP/(TP+FN)      %
%--------------------------------------------------%
EPS=1e-200;              % Epsilon (zero value)无穷小
Num_Threshold=256;
betaSqr=0.3;
Precision=zeros(imgNum,Num_Threshold);
Recall=zeros(imgNum,Num_Threshold);
TPR=zeros(imgNum,Num_Threshold);
FPR=zeros(imgNum,Num_Threshold);
MAE=zeros(1,imgNum);

for i=1:imgNum
    imgGT_gray=imread(strcat(Gt_path,Gt_list(i).name));     % 读取ground-truth
%     disp(Gt_list(i).name);
    [rows,cols]=size(imgGT_gray);   % 获取单通道图像的大小,ground-truth(单通道)
    imgSal_gray=double(imread(strcat(Salcut_path,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图
    %% 调整分辨率
    imgSal_gray1=double(imread(strcat(Salcut_path1,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图
    imgSal_gray2=double(imread(strcat(Salcut_path2,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图
    imgSal_gray3=double(imread(strcat(Salcut_path3,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图
    imgSal_gray4=double(imread(strcat(Salcut_path4,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图

    
    if size(imgSal_gray1,1)~=rows
        imgGT_gray=imgGT_gray(1:5:end,1:5:end);
    end
        S=a(1).*imgSal_gray+a(2).*imgSal_gray1+a(3)*imgSal_gray2+a(4).*imgSal_gray3+a(5)*imgSal_gray4;
        S(isnan(S))=0;
%         imgSal_gray=(S-min(S(:)))./range(S(:));
%         imgSal_gray=im2uint8(imgSal_gray);
        lll=(S-min(S(:)))./range(S(:));
        
        Tname=Salcut_list(i).name;
        Tname1=split(Tname,'.');
        Tname1=Tname1{1};
        imwrite(lll,['./Fuse2D/',Tname1,'.png']);
        
end