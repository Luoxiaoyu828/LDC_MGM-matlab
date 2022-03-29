function [outcat_simu, outcat_find, peak,beta0, x, y] = get_beta0_filter(outcat_s_f,mask,data_3d,idx)
% 对数据进行滤波
% idx 代表outcat的第idx行
outcat_simu = outcat_s_f(idx,1:18);
outcat_find = outcat_s_f(idx,19:31);
% peak = outcat_find(12);
cen1 = outcat_find(5);
cen2 = outcat_find(6);
cen3 = outcat_find(7);
size1 = outcat_find(8);
size2 = outcat_find(9);
size3 = outcat_find(10);
size_e = 3;
h = ones(size_e,size_e,size_e)/(size_e^3);
data_3d_filter = convn(data_3d,h,'same');

beta0 = [cen2;size2;cen1;size1;0;cen3;size3];

mask1=mask;
% clump_id代表编号为clump_id的云核
clump_id = outcat_s_f(idx,19);
mask1(mask1~=clump_id)=0;
% 将clump_id编号的云核数据剥离出来
data_1 = data_3d_filter.*mask1./clump_id;
% data_1 = data_1 + data_3d.*(16-mask);
% work_plot(data_1)
% 均值滤波后求peak
data_2 = data_3d_filter.*mask1./clump_id;
peak = max(data_2(:));

[size_x,size_y,size_z]=size(data_1);
[p_i, p_j, p_k] = meshgrid(1:size_y, 1:size_x, 1:size_z);
x=[p_j(:), p_i(:), p_k(:)];
index_x = find(data_1>0);
x = x(index_x,[1 2 3]);
y= data_1(:);
y = y(index_x);

end