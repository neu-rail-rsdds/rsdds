close all;
clear;
clc;

mkdir('Final_saliency');

Salcut_path='.\GLRNNRD_saliency\';
Salcut_path1='.\Depth_saliency\';
Gt_list=dir(strcat(Salcut_path,'*.png'));             % 单通道

a=ones(3,1);

Salcut_list=Gt_list;
for i=1:length(Gt_list)
Salcut_list(i).name=[Salcut_list(i).name(1:end-4),'.png'];
Salcut_list1(i).name=[Salcut_list(i).name(1:end-4),'.png'];
end
imgNum=length(Salcut_list);                       % 图片数量

if imgNum~=length(Gt_list)
    error('图片数量不等');
end

for i=1:imgNum
    imgSal_gray=double(imread(strcat(Salcut_path,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图
    imgSal_gray1=double(imread(strcat(Salcut_path1,Salcut_list1(i).name)));   % 读取显著图,显著图为单通道灰度图
    S=a(1).*imgSal_gray+a(2).*imgSal_gray1+a(3).*sqrt(imgSal_gray.*imgSal_gray1);
        S(isnan(S))=0;
        lll=(S-min(S(:)))./range(S(:));
        Tname=Salcut_list(i).name;
        Tname1=split(Tname,'.');
        Tname1=Tname1{1};
        imwrite(lll,['./Final_saliency/',Tname1,'.png']);
        
end