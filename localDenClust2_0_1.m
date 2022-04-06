function [outcat, out, mask, Gradient] = localDenClust2_0_1(data,para,is_plot,get_out_mask)
Gradient = 0;
%   INPUT:
%       data: 3D data or 2D data
%       xx: 3D or 2D data coordinates
%       para: clustering parameters
%                 para.rhomin: Minimum density
%                 para.deltamin: Minimum delta
%                 para.v_min: Minimum volume
%                 para.rms: The noise level of the data, used for data truncation calculation
%                 para.sigma: Standard deviation of Gaussian filtering
%       is_plot: 1 denotes that plot the Decision Graph, otherwise 0.
%       get_out_mask: 1 denotes that return out and mask, otherwise 0.
%   OUTPUT:
%       outcat: description of the clump:
%               Peak1,Peak2,Peak3,Cen1,Cen2,Cen3,Size1,Size2,Size3,Sum,Peak,Volume
%       out: data obtained through mask
%       mask:  1,2,3,…
% 2021/05/25  修改  保留密度和距离最低阈值作为聚类中心的判别准则
addpath E:\local_density_clustering\model\local_density_cluster_function
tic
if length(size(data))<=2
    disp('2d ok')
    [size_x,size_y]=size(data);
    [p_i, p_j] = meshgrid(1:size_y, 1:size_x);
    xx=[p_j(:), p_i(:)];
    [NClust,centInd,clustInd,~,mask,out,~]=densityCluster_2d(data,xx,para,is_plot);
end
if length(size(data))==3
    disp('3d ok')
    [size_x,size_y,size_z]=size(data);
    [p_i, p_j, p_k] = meshgrid(1:size_y, 1:size_x, 1:size_z);
    xx=[p_j(:), p_i(:), p_k(:)];
    [NClust,centInd,clustInd,~,mask,out,~]=densityCluster_3d(data,xx,para,is_plot,get_out_mask);
end
[outcat] = extroclump_parameters(NClust,xx,clustInd,centInd,data);
toc
