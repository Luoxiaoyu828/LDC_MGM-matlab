function [outcat] = extroclump_parameters(NCLUST,xx,clustInd,centInd,data)
% Input:
% NCLUST: 检测算法得到的类别个数
% data: 聚类的数据
% clustInd: 聚类的编号，同一个类用一个相同的数字标记
% xx: data的坐标系
% centInd: centInd = [centIndex, clust_id] 代表聚类中心点在data中的索引以及聚类的编号(ID)
% 
% Output:
% Peak1,Peak2,Peak3,Cen1,Cen2,Cen3,Size1,Size2,Size3,Sum,Peak,Volume


dim = size(xx,2);
if dim==3
    clustSum = zeros(1,NCLUST);
    clustVolume = zeros(1,NCLUST);
    clustPeak = zeros(NCLUST,1);
    clump_Cen = zeros(NCLUST,dim);
    clustSize = zeros(NCLUST,dim);
    [size_x,size_y,size_z] = size(data);
    
    clump_Peak = xx(centInd(:,1),:);
    precent=1;
%     cl_result = zeros(size(data));
    for i = 1 : NCLUST
        cl_i = zeros(size_x,size_y,size_z);
        cl_1_index_= xx((clustInd==centInd(i,2)),:);
       
        clustNum=size(cl_1_index_,1);

        for j =1:clustNum
            cl_i(cl_1_index_(j,1), cl_1_index_(j,2), cl_1_index_(j,3))=1;
        end
        
        CC = bwconncomp(cl_i);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        big_clust_idx = CC.PixelIdxList{idx};
        
        cl_1_index=xx(big_clust_idx,:);  % 尺寸：clustNum * 3
        cl_i=zeros(size_x,size_y,size_z);
        clustNum=size(cl_1_index,1);
        clump_sum_=zeros(clustNum,1); % 尺寸：clustNum * 1
        for j =1:clustNum
            cl_i(cl_1_index(j,1), cl_1_index(j,2), cl_1_index(j,3))=1;
            clump_sum_(j,1)=data(cl_1_index(j,1), cl_1_index(j,2), cl_1_index(j,3));
        end
        
        clustsum = sum(clump_sum_);
        clump_Cen(i,:) = clump_sum_' * cl_1_index / clustsum;
        clustVolume(i) = clustNum;
        clustSum(i) = clustsum;
        %     peak是在滤波之后的数据得到的
        %     clustPeak(i)=data_filter(cluster_points(i,1),cluster_points(i,2),cluster_points(i,3));
        clustPeak(i) = data(clump_Peak(i,1),clump_Peak(i,2),clump_Peak(i,3));
        x_i = cl_1_index - clump_Cen(i,:);
        clustSize(i,:)=sqrt(((clump_sum_'*x_i(:,:).^2) / clustsum)-(clump_sum_'*x_i(:,:)./clustsum).^2);
        
        temp = cl_i .* data;
        temp_sort= sort(temp(:),'descend');
        
        inde_end = round(precent*size(cl_1_index,1));
        temp_mean=mean(temp_sort(inde_end));
        temp(temp<temp_mean)=0;
        temp(temp>=temp_mean)=i;
%         cl_result=cl_result+temp;
    end
else
    clustSum=zeros(1,NCLUST);
    clustVolume=zeros(1,NCLUST);
    clustPeak=zeros(NCLUST,1);
    clump_Cen = zeros(NCLUST,dim);
    clustSize = zeros(NCLUST,dim);
    [size_x,size_y]=size(data);
    
    clump_Peak=xx(centInd(:,1),:);
    precent=1;
    cl_result = zeros(size(data));
    for i =1:NCLUST
        cl_i=zeros(size_x,size_y);
        cl_1_index_=xx((clustInd==centInd(i,2)),:);
        %     cl_1_index_=xx((clustInd==centInd(i,2)),:);
        
        clustNum=size(cl_1_index_,1);
        for j =1:clustNum
            cl_i(cl_1_index_(j,1), cl_1_index_(j,2))=1;
        end
        CC = bwconncomp(cl_i);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        big_clust_idx = CC.PixelIdxList{idx};
        
        cl_1_index=xx(big_clust_idx,:);
        cl_i=zeros(size_x,size_y);
        clustNum=size(cl_1_index,1);
       
        clump_sum_=zeros(clustNum,1);
        for j =1:clustNum
            cl_i(cl_1_index(j,1), cl_1_index(j,2))=1;
            clump_sum_(j,1)=data(cl_1_index(j,1), cl_1_index(j,2));
        end
        
        clustsum=sum(clump_sum_);
        clump_Cen(i,:)= clump_sum_'*cl_1_index(:,:) / clustsum;
        clustVolume(i)=clustNum;
        clustSum(i)=clustsum;
        %     peak是在滤波之后的数据得到的
        %     clustPeak(i)=data_filter(cluster_points(i,1),cluster_points(i,2),cluster_points(i,3));
        clustPeak(i)=data(clump_Peak(i,1),clump_Peak(i,2));
        x_i=cl_1_index-clump_Cen(i,:);
        clustSize(i,:)=sqrt(((clump_sum_'*x_i(:,:).^2) / clustsum)-(clump_sum_'*x_i(:,:)./clustsum).^2);
        
        temp=cl_i.*data;
        temp_sort=sort(temp(:),'descend');
        
        inde_end = round(precent*size(cl_1_index,1));
        temp_mean=mean(temp_sort(inde_end));
        temp(temp<temp_mean)=0;
        temp(temp>=temp_mean)=i;
        cl_result=cl_result+temp;
    end
end
% disp(size(cluster_points))
% disp(size(clump_Cen))
% disp(size(clustSize))
% disp(size(clustPeak))
% disp(size(clustSum'))
% disp(size(clustVolume'))
outcat=[clump_Peak,clump_Cen,clustSize,clustSum',clustPeak,clustVolume'];