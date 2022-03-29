function [beta0, x_, y_] = get_beta0_mix(outcat_loc,mask,data_3d,mask_idx)
% idx 代表outcat的第idx行
% low --->12  high,l_40 --->11
peak = outcat_loc(mask_idx,12);
cen1 = outcat_loc(mask_idx,5);
cen2 = outcat_loc(mask_idx,6);
cen3 = outcat_loc(mask_idx,7);
size1 = outcat_loc(mask_idx,8);
size2 = outcat_loc(mask_idx,9);
size3 = outcat_loc(mask_idx,10);
angle = zeros(length(mask_idx),1);
beta0 = [peak,cen2,size2,cen1,size1,angle,cen3,size3];


% clump_id代表编号为clump_id的云核
clump_id = outcat_loc(mask_idx,1);
% data_1 = zeros(size(data_3d));

[size_x,size_y,size_z] = size(data_3d);
[p_i, p_j, p_k] = meshgrid(1:size_y, 1:size_x, 1:size_z);
x = [p_j(:), p_i(:), p_k(:)];
y = data_3d(:);
y_ = zeros(size(y));
x_ = zeros(size(x));

for i = 1:length(mask_idx)
    mask1 = mask;
    mask1(mask1~=clump_id(i))=0;
    % 将离群点去掉
    CC = bwconncomp(mask1);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    cluster_points_index = CC.PixelIdxList{idx};  % 得到云核主要部分的点的索引
    
    y_(cluster_points_index) = y(cluster_points_index);
    x_(cluster_points_index,:) = x(cluster_points_index,:);
%     data_temp = data_3d.*mask1./clump_id(i);
%     data_1 = data_1 + data_temp;
end

% 将clump_id编号的云核数据剥离出来

% figure
% work_plot(data_1)



index_ = find(y_==0);
x_(index_,:) = [];
y_(index_) = [];
% 
% index_1 = find(y_>0);
% y_ = y_(index_1);
% x_ = x_(index_1,:);

end