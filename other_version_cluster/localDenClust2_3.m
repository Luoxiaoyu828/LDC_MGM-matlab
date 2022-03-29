function [outcat, out, mask, Gradient] = localDenClust2_3(data,para,is_plot,get_out_mask,figure_name)
%   INPUT:
%       data: 3D data or 2D data
%       xx: 3D or 2D data coordinates
%       para: clustering parameters
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

% 2021/04/22  修改
%去掉了最小密度、距离的判断准则。改为先聚类，再确定聚类边界，最后根据类成员数判断是否保留。
% 修改：对于体积滤波放在对边界确定完成后进行


addpath E:\local_density_clustering\local_density_model
tic
if length(size(data))<=2
    disp('clustering for 2d data!')
    [size_x,size_y]=size(data);
    [p_i, p_j] = meshgrid(1:size_y, 1:size_x);
    xx=[p_j(:), p_i(:)];
    [NClust,centInd,clustInd,mask,out,Gradient]=densityCluster_2d(data,xx,para,is_plot,get_out_mask);
end
if length(size(data))==3
    [size_x,size_y,size_z]=size(data);
    [p_i, p_j, p_k] = meshgrid(1:size_y, 1:size_x, 1:size_z);
    xx=[p_j(:), p_i(:), p_k(:)];
    [NClust,centInd,clustInd,~,mask,out,~]=densityCluster_3d(data,xx,para,is_plot,get_out_mask,figure_name);
end
[outcat] = findclumps(NClust,xx,clustInd,centInd,data);
toc
end
%% densityCluster 3D
function [NCLUST_,centInd,clustInd_re,cluster_info,mask,out,Gradient]=densityCluster_3d(data,xx,para,is_plot,get_out_mask,figure_name)
%   SEE the following paper published in *SCIENCE* for more details:
%       Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%       Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
%       做了部分修改,2020/12/08,vastlxy@163.com

%   INPUT:
%       data: 3D data
%       xx: 3D data coordinates
%       para: clustering parameters
%                 para.v_min: Minimum volume
%                 para.rms: The noise level of the data, used for data truncation calculation
%                 para.sigma: Standard deviation of Gaussian filtering
%       is_plot: 1 denotes that plot the Decision Graph, otherwise 0.
%   OUTPUT:
%       NCLUST: number of clusters
%       clustInd: cluster index that each point belongs to, NOTE that -1 represents no clustering assignment (haloInd points)
%       centInd:  centroid index vector
%       cluster_info: [delta(icl)',rho(icl),clustVolume]
addpath E:\local_density_clustering\local_density_model
gradtmin=para.gradtmin;
v_min=para.v_min;
rms=para.rms;

data_filter=imgaussfilt3(data,para.sigma);
[size_x,size_y,size_z] = size(data_filter);
rho=data_filter(:);

% [rho_sorted,rho_Ind]% rho_sorted是排序的结果% rho_Ind是排序的索引
[rho_sorted,rho_Ind]=sort(rho,'descend');

% 初始化
maxd=size_x+size_y+size_z;
ND=length(rho);
delta=zeros(1,ND);
IndNearNeigh=zeros(1,ND);
Gradient=zeros(1,ND);

r_region = 4;

% delta  记录距离，
% IndNearNeigh  记录：两个密度点的联系% index of nearest neighbor with higher density
delta(rho_Ind(1))=sqrt(size_x^2+size_y^2+size_z^2);   % 2021/01/18 做了修改  原来是-1
IndNearNeigh(rho_Ind(1))=ND+2;    % 2021/01/18 做了修改  原来是0
% Gradient(ordrhoInd(1))=0;

% 计算 delta, Gradient
tic
for ii=2:ND
    %密度降序排序后，即密度第ii大的索引（在rho中）
    ordrho_ii=rho_Ind(ii);
    % 密度第ii大点的密度
    rho_ii =  rho_sorted(ii);
    if rho_ii >= rms
        % 首先对delta赋一个较大的初始值
        delta(ordrho_ii)=maxd;
        %密度第ii大的点对应的坐标
        delta_ii_xy = xx(ordrho_ii,:);
        % 得到以delta_ii_xy这一点为中心的 k*k*k 的邻域内的点在原图中的坐标
        bt = kc_coord_3d(delta_ii_xy,size_z,size_y,size_x,r_region);
        for j_=1:size(bt,1)
            rho_jj= data_filter(bt(j_,1),bt(j_,2),bt(j_,3));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            if dist_i_j == 0
                continue
            end
            gradient = (rho_jj - rho_ii) / dist_i_j;
            if dist_i_j <= delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                delta(ordrho_ii)= dist_i_j/(r_region*sqrt(3)) + 0.0001;
                Gradient(ordrho_ii)= gradient;
                IndNearNeigh(ordrho_ii)=(bt(j_,3)-1)*size_x*size_y + (bt(j_,2)-1)*size_x + bt(j_,1);
            end
        end
        if delta(ordrho_ii)==maxd% 表明，在设置的"框"中没有找到比该点强度值高的点，则进行全局搜索
            for jj=1:ii-1
                rho_jj = rho_sorted(jj);
                dist_i_j=distx(ordrho_ii,rho_Ind(jj),xx);
                gradient = (rho_jj - rho_ii) / dist_i_j;
                if dist_i_j<=delta(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    IndNearNeigh(ordrho_ii)=ND+2;
                end
            end
        end
    else
        % rho_ii < para.rms则将点标记为一个特殊值
        IndNearNeigh(ordrho_ii)=ND+1;
    end
end

[detla_sorted,~] = sort(delta,'descend');
delta(rho_Ind(1))=detla_sorted(2);

toc
disp('delata, rho and Gradient are ok')

NCLUST=0;
clustInd=-1*ones(1,ND+1);

% 将局部极大值点选出来，局部极大值点：在设定的"框"中未找到比该点强度值大的点
clust_index = find(IndNearNeigh==ND+2);

clust_num = length(clust_index);
str_1 = sprintf('find local maximun points number: %d\n',clust_num);
fprintf(str_1)
% icl是用来记录第i个类中心在xx中的索引值
icl = zeros(clust_num,1);
for ii=1:clust_num
    i = clust_index(ii);  % i 表示的密度点排序后的索引，这里都均为识别为聚类中心的点
    NCLUST=NCLUST+1;
    clustInd(i)=NCLUST;
    icl(NCLUST)=i;  % icl是用来记录第i个类中心在xx中的索引值
end
% assignation  将其他非类中心分配到离它最近的类中心中去
% clustInd=-1 表示该点不是类的中心点，属于其他点，等待被分配到某个类中去
% 类的中心点的梯度Gradient被指定为-1
for i=1:ND
    ordrho_i=rho_Ind(i);
    if clustInd(ordrho_i)==-1 % not centroid
        clustInd(ordrho_i)=clustInd(IndNearNeigh(ordrho_i));
    else
        Gradient(ordrho_i)=-1;
    end
end

% 计算每个聚类的体积
clustVolume = zeros(NCLUST,1);
for i=1:NCLUST
    clustVolume(i,1)=length(clustInd(clustInd==i));
end
% the Decision Graph
if is_plot
    DecisionGraph_3D(NCLUST,rho,delta,icl)
end

cluster_info=[delta(icl)',rho(icl),clustVolume];

centInd = icl(clustVolume >= v_min);
centNum = find(clustVolume >= v_min);
% centInd [类中心点在xx坐标下的索引值， 类中心在centInd的索引值: 代表类别编号]

centInd = [centInd, centNum];
NCLUST_ = length(centNum);
str_2 = sprintf('Number of cluster(volume > %d): %d\n',v_min, NCLUST_);
fprintf(str_2)

clustInd_re = -1 * ones(size(clustInd));  % 保存最后确定下来的云核的坐标索引
mask = zeros(size_x,size_y,size_z);
out = zeros(size_x,size_y,size_z);

% 根据梯度将类的边界进一步压缩，得到云核主要部分
mask_grad = find(Gradient>gradtmin);
clust_id = 1;
for i = 1:size(centInd,1)
    
    rho_clust_i = zeros(size(rho));
    index_clust_i = find(clustInd==centInd(i,2));% 得到第i类所有点的索引
    
    index_cc =  intersect(mask_grad,index_clust_i);% 得到第i类梯度大于gradtmin的点的索引
    rho_clust_i(index_clust_i) = rho(index_clust_i); % 得到第i类所有点的密度
    
    rho_cc_mean = mean(rho(index_cc)); % 得到第i类梯度大于gradtmin的点的密度的均值
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % 得到第i类密度大于该阈值的点的索引
    
    index_clust_rho = union(index_cc, index_cc_rho); %得到两个条件都满足的点的索引
    
    clustInd_re(index_clust_rho) = centInd(i,2);
    if get_out_mask
        cl_i_point=xx(index_clust_rho,:); % 得到第i类梯度大于gradtmin的点的坐标
        %     ini = find(clustInd==1);
        %     cl_i_point=xx(ini,:); % 得到第i类梯度大于gradtmin的点的坐标
        clump_Image=zeros(size_x,size_y,size_z);
        mask_out=zeros(size_x,size_y,size_z);
        if size(cl_i_point,1)>0
            for j =1:size(cl_i_point,1)
                %         mask_clust(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=i;
                mask_out(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=1;
            end
            %     对每个类的掩模做形态学处理：闭运算，填补空洞，保留面积最大连通域
            bw = imbinarize(mask_out,0.5);
            BW2 = logical(bw);           % 默认8连通； bwlabe(bw,4);
            %     se=strel('disk',3);
            %             SE = strel('sphere', 2);
            %             BW2=imclose(BW2,SE);         % 闭运算
            BW2 = imfill(BW2,'holes');   % 填补空洞
            
            STATS = regionprops3(BW2,'All');
            % 统计上一步标记图像中的连通域的面积分布
            max_ind = find( STATS.Volume==max( STATS.Volume));
            points = STATS.VoxelList{max_ind};
            %     disp(STATS.Volume)
            for ii = 1: STATS.Volume(max_ind)
                clump_Image(points(ii,2),points(ii,1),points(ii,3))=1;
            end
            %     clump_Image = double(BW2);
            mask_clust = clump_Image .* clust_id;
            mask = mask + mask_clust;
            out = out + clump_Image.*data;
            clust_id = clust_id+1;
        end
    end
end
end

%% densityCluster 2D
function [NCLUST_,centInd,clustInd,mask,out,Gradient]=densityCluster_2d(data,xx,para,is_plot,get_out_mask)
%   SEE the following paper published in *SCIENCE* for more details:
%       Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%       Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
%       做了部分修改,2020/12/08,vastlxy@163.com
%   INPUT:
%       data: 2D data
%       xx: 2D data coordinates
%       para: clustering parameters
%                 para.v_min: Minimum volume
%                 para.rms: The noise level of the data, used for data truncation calculation
%                 para.sigma: Standard deviation of Gaussian filtering
%       is_plot: 1 denotes that plot the Decision Graph, otherwise 0.
%   OUTPUT:
%       NCLUST: number of clusters
%       clustInd: cluster index that each point belongs to, NOTE that -1 represents no clustering assignment (haloInd points)
%       centInd:  centroid index vector


gradtmin = para.gradtmin;
v_min = para.v_min;
rms = para.rms;
region = para.region;

data_filter=imgaussfilt(data,para.sigma);
[size_x,size_y,size_z] = size(data_filter);
rho=data_filter(:);

[rho_sorted,rho_Ind]=sort(rho,'descend');%  rho_sorted是排序的结果% ordrhoInd是排序的索引

% 初始化
maxd = sqrt(size_x^2+size_y^2);
ND=length(rho);
delta = zeros(1,ND);
IndNearNeigh = zeros(1,ND);
Gradient = zeros(1,ND);

% delta  记录距离，
% IndNearNeigh  记录：两个密度点的联系% index of nearest neighbor with higher density
delta(rho_Ind(1))=-1;
IndNearNeigh(rho_Ind(1))=ND+2;

% 计算 delta, Gradient
for ii=2:ND
    ordrho_ii=rho_Ind(ii);% 密度降序排序后，即密度第ii大的索引（在rho中）
    rho_ii =  rho_sorted(ii);% 当前点的密度及坐标
    if rho_ii >=rms
        delta(ordrho_ii)=maxd;
        delta_ii_xy = xx(ordrho_ii,:);
        % 得到以delta_ii_xy这一点为中心的 k*k 的邻域内的点在原图中的坐标
        bt = kc_coord_2d(delta_ii_xy, size_x, size_y, region);
        for j_=1:size(bt,1)
            rho_jj= data_filter(bt(j_,1),bt(j_,2));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            gradient = (rho_jj - rho_ii) / dist_i_j;
            if dist_i_j <= delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                delta(ordrho_ii)= dist_i_j;
                Gradient(ordrho_ii)= gradient;
                % indNearNeigh(ordRho(i)) = ordRho(j); 下面等式的右端实现的是得到ordRho(j)
                IndNearNeigh(ordrho_ii)=(bt(j_,2)-1)*size_x + bt(j_,1);
            end
        end
        if delta(ordrho_ii)==maxd
            % 表明，在3*3*3的领域中没有找到比该点高，距离最近的点，则进行全局搜索
            for jj=1:ii-1
                rho_jj = rho_sorted(jj);
                dist_i_j=distx(ordrho_ii,rho_Ind(jj),xx);
                gradient = (rho_jj - rho_ii) / dist_i_j;
                if dist_i_j<=delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    IndNearNeigh(ordrho_ii)=ND+2;
                end
            end
        end
    else     
        IndNearNeigh(ordrho_ii)=ND+1;% rho_ii < para.rms则将点标记为一个特殊值
    end
end

[detla_sorted,~] = sort(delta,'descend');
delta(rho_Ind(1))=detla_sorted(2);
disp('delata, rho and Gradient are ok')

NCLUST=0;
clustInd=-1*ones(1,ND+1);

% 将局部极大值点选出来.局部极大值点：在设定的"框"中未找到比该点强度值大的点
clust_index = find(IndNearNeigh==ND+2);

clust_num = length(clust_index);
str_1 = sprintf('find local maximun points number: %d\n',clust_num);
fprintf(str_1)

icl = zeros(clust_num,1);% icl是用来记录第i个类中心在xx中的索引值
for ii=1:clust_num
    i = clust_index(ii);
    NCLUST=NCLUST+1;
    clustInd(i)=NCLUST;
    icl(NCLUST)=i;
end

% the Decision Graph
if is_plot
    DecisionGraph(NCLUST,rho,delta,icl,clustInd,clust_index)
end

% assignation  将其他非类中心分配到离它最近的类中心中去
% clustInd=-1 表示该点不是类的中心点，属于其他点，等待被分配到某个类中去
% 类的中心点的梯度Gradient被指定为-1
for i=1:ND
    ordrho_i=rho_Ind(i);
    if clustInd(ordrho_i)==-1 % not centroid
        clustInd(ordrho_i)=clustInd(IndNearNeigh(ordrho_i));
    else
        Gradient(ordrho_i)=-1;
    end
end

clustVolume = zeros(NCLUST,1);
for i=1:NCLUST
    clustVolume(i,1)=length(clustInd(clustInd==i));
end


centInd1 = icl(clustVolume >= v_min);
centNum = find(clustVolume >= v_min);
% centInd [类中心点在xx坐标下的索引值centInd1， 类中心在centInd的索引值centNum: 代表类别编号]
centInd = [centInd1, centNum];
NCLUST_ = length(centNum);
str_2 = sprintf('Number of cluster(volume > %d): %d\n',v_min, NCLUST_);
fprintf(str_2)
clustInd_re = -1 * ones(size(clustInd));  % 保存最后确定下来的云核的坐标索引
mask = zeros(size_x,size_y,size_z);
out = zeros(size_x,size_y,size_z);

% 根据梯度将类的边界进一步压缩，得到云核主要部分
mask_grad = find(Gradient>gradtmin);

for i = 1:size(centInd,1)
    rho_clust_i = zeros(size(rho));
    index_clust_i = find(clustInd==centInd(i,2));% 得到第i类所有点的索引
    index_cc =  intersect(mask_grad,index_clust_i);% 得到第i类梯度大于gradtmin的点的索引
    rho_clust_i(index_clust_i) = rho(index_clust_i); % 得到第i类所有点的密度
    
    rho_cc_mean = mean(rho(index_cc)); % 得到第i类梯度大于gradtmin的点的密度的均值
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % 得到第i类密度大于该阈值的点的索引
    
    index_clust_rho = union(index_cc, index_cc_rho); %得到两个条件都满足的点的索引
    
    clustInd_re(index_clust_rho) = centInd(i,2); % 重新编上类别号
end

for i = 1:size(clustInd_re,1)
           
    cl_i_point=xx(index_clust_rho,:); % 得到第i类梯度大于gradtmin的点的坐标
    if get_out_mask
        mask_out=zeros(size_x,size_y,size_z);
        for j =1:size(cl_i_point,1)
            %         mask_clust(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=i;
            mask_out(cl_i_point(j,1), cl_i_point(j,2))=1;
        end
        %     对每个类的掩模做形态学处理：闭运算，填补空洞，保留面积最大连通域
        bw = imbinarize(mask_out,0.5);
        L21 = logical(bw);  % 默认8连通； bwlabe(bw,4);
        
        BW2 = imfill(L21,'holes'); % 填补空洞
        
        L = bwlabel(BW2);   % 对连通区域进行标记
        STATS = regionprops(BW2,'All');
        % 统计上一步标记图像中的连通域的面积分布
        Area = cat(1, STATS.Area);
        ind = find(Area ==max(Area));%找到最大连通区域的标号
        BW2(L~=ind)=0;%将其它区域置为0
        
        clump_Image = double(BW2);
        mask_clust = clump_Image .* i;
        mask = mask + mask_clust;
        out = out + clump_Image.*data;
    end
end
end
