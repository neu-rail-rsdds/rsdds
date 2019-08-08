% function [P,R,Tpr,Fpr,FMeasure]=Evalution(Gt_path,Salcut_path) %带返回值
%本程序的功能是评价显著性算法的优劣(每次执行一种方法的结果)
%评价指标：PR曲线、ROC曲线、MAE、FMeasure
%评价指标：AUC(Area Under ROC Curve)、MeanFMeasure、MaxFMeasure
%获取文件夹下某一格式的所有图像
close all;
clear;
clc;
tic;
% Gt_path='.\MSRA-B\';
Gt_path='.\Ground Truth\';
% Salcut_path='.\CA\';
Salcut_path='.\Final_saliency\';
Gt_list=dir(strcat(Gt_path,'*.png'));             % 单通道

Salcut_list=Gt_list;
for i=1:length(Gt_list)
Salcut_list(i).name=[Salcut_list(i).name(1:end-4),'.png'];
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
    disp(Gt_list(i).name);
    [rows,cols]=size(imgGT_gray);   % 获取单通道图像的大小,ground-truth(单通道)
    imgSal_gray=uint8(imread(strcat(Salcut_path,Salcut_list(i).name)));   % 读取显著图,显著图为单通道灰度图
    %% 调整分辨率
    imgSal_gray=imresize(imgSal_gray,size(imgGT_gray));
%     imgSal=imread(strcat(Salcut_path,Salcut_list(i).name));   % 读取显著图,显著图为单通道灰度图
%     imgSal_gray=rgb2gray(imgSal);                             % 显著图转灰度图
    imgSal_gray_double=im2double(imgSal_gray);                % 显著性像素范围[0,1]
%     imgGT_double=double(imgGT_gray)./255; %归一化处理[0,1]
    thresh = graythresh(imgGT_gray);
%     for j = 1:255
        imgGT_bw=imgGT_gray> thresh ;
%         imgGT_bw=imgGT_gray>128;                                       % ground-truth二值化
        imgGT_not=~imgGT_bw;                    % 取非
        G=length(find(imgGT_bw==1));            % ground-truth目标mask大小
        G_not=length(find(imgGT_not==1)); 
        for thr=1:Num_Threshold
            imgSal_bw=imgSal_gray >= (thr-1);   % 固定阈值二值化
            M=length(find(imgSal_bw==1));
            tpM=imgGT_bw & imgSal_bw;           % 显著图与ground-truth的交集
            fpM=imgGT_not & imgSal_bw;          % 显著图与ground-truth补集的交集
            TP=length(find(tpM==1)); 
            FP=length(find(fpM==1));
    %         Precision(i,thr)=TP/M;
            Precision(i,thr)=TP/(M+EPS);
            Recall(i,thr)=TP/(G+EPS);       
            TPR(i,thr)=TP/(G+EPS);      
            FPR(i,thr)=FP/(G_not+EPS);      
        end
        tmp=abs(imgSal_gray_double-imgGT_bw);
        MAE(i)=sum(tmp(:))/(rows*cols);          % 第i张图片的平均绝对误差
        
%     end
end

% result(P/R/Tpr/Fpr/Mae)
P=sum(Precision)/imgNum;  % 精确率
R=sum(Recall)/imgNum;     % 召回率
Tpr=sum(TPR)/imgNum;      % 真阳性率
Fpr=sum(FPR)/imgNum;      % 假阳性率
Mae=sum(MAE)/imgNum;      % 平均绝对误差
% 其他指标(FMeasure、MeanFMeasure、MaxFMeasure)
FMeasure=((1+betaSqr)*(P.*R))./(betaSqr*P+R);  % F_beta
MeanFMeasure=sum(FMeasure)/numel(FMeasure);    % MeanFMeasure
MaxFMeasure=max(FMeasure);                     % MaxFMeasure
% 计算ROC曲线下的面积AUC
AUC_tmp=(Fpr(1,1:end-1)-Fpr(1,2:end)).*(Tpr(1,1:end-1)+Tpr(1,2:end));
AUC=sum(AUC_tmp)/2;                            % AUC
%----------------法二----------------------%
% AUC=0;
% for t=2:Num_Threshold
%     AUC=AUC+((Fpr(t-1)-Fpr(t))*(Tpr(t-1)+Tpr(t)))/2;
% end
%------------------------------------------%
%%
% 绘制评价曲线
%-----------PR曲线--------------%
figure(1);
% plot(R,P,'k*','LineWidth',2); 
plot(R,P,'linewidth', 0.5); 
xlabel('Recall');
ylabel('Precision');
% legend('SR');       % 以RC为例
grid on;
axis([0 1 0 1]);
title('Precision recall curve');
%-----------ROC曲线--------------%
figure(2);
% plot(Fpr,Tpr,'b*','LineWidth',2); 
plot(Fpr,Tpr,'linewidth', 0.5); 
xlabel('False positive rate');
ylabel('True positive rate');
% legend('SR');       % 以RC为例
grid on;
axis([0 1 0 1]);
title('ROC curve');
%-----------性能评估-------------%
figure(3);
% methodLabels=categorical({'RC'});  %以RC为例
bar([MeanFMeasure,MaxFMeasure,AUC]);
legend('Mean F_\beta', 'Max F_\beta', 'AUC');
grid on;
title('性能指标');
%%
%%
% 保存计算结果
save(Salcut_path(1:end-1),'R','P','Fpr','Tpr','Mae','MeanFMeasure','MaxFMeasure','AUC' );
toc;

