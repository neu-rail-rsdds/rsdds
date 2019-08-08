clear all;
clc;
close all;
name1=ls('.\INPUT_IMG\*.bmp');
name2=ls('.\Depth_saliency\*.png');
num_l=length(name1);
for iii=1:num_l
    Img=imread(['.\INPUT_IMG\',name1(iii,:)]);
    Saliency_D=double(imread(['.\Depth_saliency\',name2(iii,:)]));
    %%
    T=[200,300,350,400,700]; %% Multi-scale superpixel segmentation
    [m,n,~]=size(Img);
    S=zeros([m,n]);
    for j=1:2:length(T)
        SL=NLNL_FUSE(Img,Saliency_D,T(j));
        S=S+SL;
    end
    %%
    S(isnan(S))=0;
    lll=(S-min(S(:)))./range(S(:));
    %%
    Tname=name1(iii,:);
    Tname1=split(Tname,'.');
    Tname1=Tname1{1};
    imwrite(lll,['./GLRNNRD_saliency/',Tname1,'.png']);
end