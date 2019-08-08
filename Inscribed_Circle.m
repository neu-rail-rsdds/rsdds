function [P,R_p,Circle_1]=Inscribed_Circle(BW)

%% 随机撒500个点，计算里面最大内接圆点坐标
%% 到每个区间的凸多边形轮廓点位置的距离，X最近点，Y最近点，300X500计算量
    [m,n]=size(BW);
    num=15000;  
    mm=rand(1,num)*m;
    mm(round(mm)==0)=1;
    nn=rand(1,num)*n;
    nn(round(nn)==0)=1;
    pos=[round(mm);round(nn)];
    pos=pos';
    Poss=sub2ind2([m,n],pos);

    %% 删除 BW内点，
    [Bx,By]=meshgrid(1:n,1:m);
    idx=[Bx(1,:),Bx(end,:),Bx(:,1)',Bx(:,end)'];
    idy=[By(1,:),By(end,:),By(:,1)',By(:,end)'];
    P_id=[idy',idx'];
    PP=unique(P_id,'rows');
    PP(PP(:,1)==0,:)=[];
    
    %%%
    Bx=Bx.*BW;
    By=By.*BW;
    %%%%%  提取四周边界坐标   

    %%
    pos_B=[By(:),Bx(:)];
    pos_B(Bx(:)==0,:)=[];
    Poss_B=sub2ind2([m,n],pos_B);
    pos_xy=setdiff(Poss,Poss_B);
    %%
    [pos_Bx,pos_By]=ind2sub([m,n],pos_xy);
    
    %% 计算 最大内接圆
    %%%%% 计算 各区域质心 到 随机点的距离
    pos_Bxy=[pos_Bx,pos_By];
    
    B = bwboundaries(BW);
    B=cell2mat(B);
    Dist_B=pdist2(pos_Bxy,[B;PP],'euclidean');
    
    %%%% 筛选距离中的 最小 最大距离
    MD_B=min(Dist_B');
    pos_MDB=find(MD_B==max(MD_B));
    P=pos_Bxy(pos_MDB(1),:);
    R_p=max(MD_B);
    P_CR=P;
    R_pp=R_p;
    
    BW1=BW;
    Circle_1=zeros(size(BW));
    [Bx,By]=meshgrid(1:n,1:m);
%     BW1((Bx-P(2)).^2+(By-P(1)).^2<max(MD_B).^2)=1;
    BW1((Bx-P(2)).^2+(By-P(1)).^2<=R_p.^2)=1;
    Circle_1((Bx-P(2)).^2+(By-P(1)).^2<=R_p.^2)=1;
    end