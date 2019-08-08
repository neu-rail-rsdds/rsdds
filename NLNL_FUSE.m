%% 获取初始结果

%% the code of "the entropy of gridents of Depth for Saliency"
%% 加载图像
% clear all;
% clc;
% num=1; %%rows: 2375 ; cols:3791;
% Img=imread(['.\INPUT_IMG\',num2str(num),'.bmp']);
% 

function lll=NLNL_FUSE(Img,Saliency_D,spnum)

[height,width]=size(Img);
[mm,nn]=size(Saliency_D);
if height~=mm 
    Img=Img(1:5:end,1:5:end,:);
end

[sup.label, sup.num]=superpixels(im2double(Img),spnum,'Compactness',30);
%% 边缘_背景
L1=sup.label;
L_B=1+0.*L1;
L_B(2:end-1,2:end-1)=0;
LB_Idx=L_B.*L1;
LB_Idx=unique(LB_Idx(:));
LB_Idx(LB_Idx(:)==0)=[];
bsnum=unique(LB_Idx);
supNum = max(sup.label(:));

%% 特征提取

  %% STEP-3. Extract features
    % get superpixel statistics
sup.pixIdx = cell(sup.num, 1);
sup.pixNum = zeros(sup.num,1);
for i = 1:sup.num
     temp = find(sup.label==i);
     sup.pixIdx{i} = temp;
     sup.pixNum(i) = length(temp);
end
    %%
addpath('./extractFeatures'); 

% featImg = ExtractFeature(im2single(img.RGB));
featImg = ExtractFeature(im2single(Img));
for i = 1:3
    featImg(:,:,i) = mat2gray(featImg(:,:,i)).*255;
end
featImg(:,:,end+1)=Saliency_D;
featMat = GetMeanFeat(featImg, sup.pixIdx);  
featMat = featMat./255;
colorFeatures = featMat(:,1:3);
medianR = median(colorFeatures(:,1)); medianG = median(colorFeatures(:,2)); medianB = median(colorFeatures(:,3));
featMat(:,1:3) = (featMat(:,1:3)-1.2*repmat([medianR, medianG, medianB],size(featMat,1),1))*1.5;


%% 全局低秩



%%
    addpath(genpath('./inexact_alm_rpca'));
% [A_hat E_hat iter] = inexact_alm_rpca(featMat, 0.12);
%     [A_hat E_hat iter] = inexact_alm_rpca(featMat, 0.12);
    
    % ||L||*+||Z||*+a*||H1||+b*||S1||+c*||H2||+d*||S2||
    % s.t. F=L+S1,L=H1; L=A*Z+S2;Z=H2,H2>=0
    %%加入边界惩罚项
    %s.t. S3=F(:,bg)-A;
    %%
    a=[0,1.2,0.02,1.2];
    [L,S1,AA,Z2,S2]=TL_Bg_NNLRSR(featMat,bsnum,a);
    %%
    recError = sum((S1+S2).^2);
    recError = (recError- min(recError(:)))/(max(recError(:)) - min(recError(:)));
    lll=0.*L1;
%     ll=sum(abs(E)./sum(E));
%     ll=sum(abs(E));
    for i=1:supNum
    lll(L1==i)=recError(i);
%    lll(L1==i)=ll(i);
    end
%     imshow(lll,[]); 
 