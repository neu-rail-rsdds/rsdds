%% Calculating a two-dimensional saliency map
clear all;
clc;
close all;
name1=ls('.\INPUT_IMG\*.bmp');
num_l=length(name1);
mkdir('GLRNNR_saliency');
for iii=1:num_l
    Img=imread(['.\INPUT_IMG\',name1(iii,:)]);
    T=[200,300,350,400,700]; %% Multi-scale superpixel segmentation
    [m,n,~]=size(Img);
    S=zeros([m,n]);
    for j=1:2:length(T)
        [SL,time]=NLNL_SPR(Img,T(j));
        S=S+SL;
        disp(time);
    end
    S(isnan(S))=0;
    lll=(S-min(S(:)))./range(S(:));
    Tname=name1(iii,:);
    Tname1=split(Tname,'.');
    Tname1=Tname1{1};
    imwrite(lll,['./GLRNNR_saliency/',Tname1,'.png']);
end