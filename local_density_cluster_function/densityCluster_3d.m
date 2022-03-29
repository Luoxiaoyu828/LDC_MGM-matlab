function [NCLUST_,centInd_re,clustInd_re,cluster_info,mask,out,Gradient]=densityCluster_3d(data,xx,para,is_plot,get_out_mask)
%% densityCluster 3D
%   SEE the following paper published in *SCIENCE* for more details:
%       Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%       Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
%       做了部分修改,2020/12/08,vastlxy@163.com

%   INPUT:
%       data: 3D data
%       xx: 3D data coordinates
%       para: clustering parameters
%                 para.rhomin: Minimum density
%                 para.deltamin: Minimum delta
%                 para.v_min: Minimum volume
%                 para.rms: The noise level of the data, used for data truncation calculation
%                 para.sigma: Standard deviation of Gaussian filtering
%       is_plot: 1 denotes that plot the Decision Graph, otherwise 0.
%   OUTPUT:
%       NCLUST: number of clusters
%       clustInd: cluster index that each point belongs to, NOTE that -1 represents no clustering assignment (haloInd points)
%       centInd:  centroid index vector
%       cluster_info: [delta(icl)',rho(icl),clustVolume]

% 参数初始化
gradtmin = para.gradtmin;
rhomin = para.rhomin;
deltamin = para.deltamin;
v_min = para.v_min;
rms = para.rms;
if para.sigma == 0
    data_filter = data;
else
    data_filter = imgaussfilt3(data,para.sigma);
end
[size_x,size_y,size_z] = size(data_filter);
rho=data_filter(:);

% [rho_sorted,rho_Ind]% rho_sorted是排序的结果% rho_Ind是排序的索引
[rho_sorted,rho_Ind] = sort(rho,'descend');

% 初始化
maxd = size_x + size_y + size_z;
ND = length(rho);
delta = zeros(1,ND);
IndNearNeigh = zeros(1,ND);
Gradient = zeros(1,ND);

% delta  记录距离，
% IndNearNeigh  记录：两个密度点的联系% index of nearest neighbor with higher density
delta(rho_Ind(1)) = sqrt(size_x^2+size_y^2+size_z^2);   % 2021/01/18 做了修改  原来是-1
IndNearNeigh(rho_Ind(1)) = rho_Ind(1);    % 2021/01/18 做了修改  原来是0

% 计算 delta, Gradient
tic
for ii = 2:ND
    %密度降序排序后，即密度第ii大的索引（在rho中）
    ordrho_ii = rho_Ind(ii);
    rho_ii = rho_sorted(ii);% 第ii的密度
    if rho_ii >= rms
        % 首先对delta赋一个较大的初始值
        delta(ordrho_ii) = maxd;
        %密度第ii大的点对应的坐标
        delta_ii_xy = xx(ordrho_ii,:);
        % 得到以delta_ii_xy这一点为中心的 (2*k+1)*(2*k+1)*(2*k+1) 的邻域内的点在原图中的坐标
        k = 2;
        bt = kc_coord_3d(delta_ii_xy,size_z,size_y,size_x,k);
        for j_ = 1:size(bt,1)
            rho_jj = data_filter(bt(j_,1),bt(j_,2),bt(j_,3));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            if dist_i_j == 0
                continue
            end
            gradient = (rho_jj - rho_ii) / dist_i_j;
            % 通过距离和梯度值的控制  可以实现找到“比当前点强度大的点”中最近的点
            if dist_i_j <= delta(ordrho_ii) && gradient >= 0
                delta(ordrho_ii)= dist_i_j;
                Gradient(ordrho_ii)= gradient;
                % indNearNeigh(ordRho(i)) = ordRho(j); 下面等式的右端实现的是得到ordRho(j)
                IndNearNeigh(ordrho_ii)=(bt(j_,3)-1)*size_x*size_y + (bt(j_,2)-1)*size_x + bt(j_,1);
            end
        end
        % 表明，在(2*k+1)*(2*k+1)*(2*k+1)的领域中没有找到比该点高，距离最近的点，则进行全局搜索
        if delta(ordrho_ii)==maxd
            % 因为点的强度降序排列了  所以jj的点总是比ii的点强度大，故在下面的判断中不用梯度判断
            for jj=1:ii-1
                rho_jj = rho_sorted(jj);
                dist_i_j=distx(ordrho_ii,rho_Ind(jj),xx);  % 计算两点间的距离
                gradient = (rho_jj - rho_ii) / dist_i_j;    % 计算两点间的梯度
                if dist_i_j<=delta(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    
                    % 如果两个点的距离大于4,则标记为中心点侯选体
                    if dist_i_j <= 4
                        IndNearNeigh(ordrho_ii)=rho_Ind(jj);
                    else
                        IndNearNeigh(ordrho_ii)=ND+1;
                    end
                    
                end
            end
        end
    else
        % rho_ii < para.rms则将点标记为一个特殊值
        IndNearNeigh(ordrho_ii)=ND+1;
    end
    
end

[detla_sorted,~] = sort(delta,'descend');
% [rho_sorted,~] = sort(rho,'descend');
% delta(rho_Ind(1))=max(delta(:));
delta(rho_Ind(1))=detla_sorted(2);

toc
disp('delata, rho and Gradient are ok')

% 根据密度和距离来确定类中心
NCLUST=0;
clustInd=-1*ones(1,ND+1);

disp('find cluster center point by delata and rho')

% 选取类中心点：将满足密度，距离最低值的点选出来作为类中心点
clust_index = find(and(rho>rhomin,delta'>deltamin)>0);
clust_num = length(clust_index);

% icl是用来记录第i个类中心在xx中的索引值
icl = zeros(clust_num,1);
for ii=1:clust_num
    i = clust_index(ii);
    NCLUST=NCLUST+1;
    clustInd(i)=NCLUST;
    icl(NCLUST)=i;
end

% the Decision Graph
if is_plot
    figure
    plot(rho(:),delta(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
    title ('Decision Graph','FontSize',15.0)
    xlabel ('\rho')
    ylabel ('\delta')
    for i=1:NCLUST
        hold on
        %    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',6);
    end
end
% 
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

cluster_info=[delta(icl)',rho(icl),clustVolume];

centInd = icl(clustVolume >= v_min);
centNum = find(clustVolume >= v_min);
% centInd [类中心点在xx坐标下的索引值， 类中心在centInd的索引值: 代表类别编号]
centInd = [centInd, centNum];

clustInd_re = -1 * ones(size(clustInd));  % 保存最后确定下来的云核的坐标索引
mask = zeros(size_x,size_y,size_z);
out = zeros(size_x,size_y,size_z);

% 根据梯度将类的边界进一步压缩，得到云核主要部分
mask_grad = find(Gradient>gradtmin);

% 根据体积会去掉一部分核 clust_id为重新的编号  centInd_re为新的云核中心的索引
clust_id = 1;
centInd_re = [];

for i = 1:size(centInd,1)
    
    rho_clust_i = zeros(size(rho));
    index_clust_i = find(clustInd==centInd(i,2));% 得到第i类所有点的索引
    
    index_cc =  intersect(mask_grad,index_clust_i);% 得到第i类梯度大于gradtmin的点的索引
    rho_clust_i(index_clust_i) = rho(index_clust_i); % 得到第i类所有点的密度
    
    rho_cc_mean = 0.2 * mean(rho(index_cc)); % 得到第i类梯度大于gradtmin的点的密度的均值的0.2倍
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % 得到第i类密度大于该阈值的点的索引
    
    index_clust_rho = union(index_cc, index_cc_rho); %得到两个条件都满足的点的索引
    
    % 判断同时满足两个条件的云核最终的体积是否足够大
    if length(index_clust_rho) > v_min
        clustInd_re(index_clust_rho) = clust_id; % 对保留的云核重新编号
        centInd_re = [centInd_re; [centInd(i,1), clust_id]];
        if get_out_mask
            cl_i_point=xx(index_clust_rho,:); %得到第i类梯度大于gradtmin和强度大于原始强度均值的点的坐标
            clump_Image=zeros(size_x,size_y,size_z);
            
            mask_out=zeros(size_x,size_y,size_z);
            for j =1:size(cl_i_point,1)
                mask_out(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=1;
            end
            % 对每个类的掩模做形态学处理：填补空洞，保留面积最大连通域
            bw = imbinarize(mask_out,0.5);
            BW2 = logical(bw);           % 默认8连通； bwlabe(bw,4);
            BW2 = imfill(BW2,'holes');   % 填补空洞
            
            STATS = regionprops3(BW2,'All'); % 统计上一步标记图像中的连通域的面积分布
           
            max_ind = find( STATS.Volume==max( STATS.Volume));
            points = STATS.VoxelList{max_ind};

            for ii = 1: STATS.Volume(max_ind)
                clump_Image(points(ii,2),points(ii,1),points(ii,3))=1;
            end
            mask_clust = clump_Image .* clust_id;
            mask = mask + mask_clust; % 得到掩模
            out = out + clump_Image.*data; % 用掩模得到的数据
        end
        clust_id = clust_id+1; % 对云核的编号进行更新
    end
end

NCLUST_ = size(centInd_re,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);

