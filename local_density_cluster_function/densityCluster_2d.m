function [NCLUST_,centInd,clustInd,cluster_info,mask,out,Gradient]=densityCluster_2d(data,xx,para,is_plot)
%   SEE the following paper published in *SCIENCE* for more details:
%       Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%       Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
%       做了部分修改,2020/12/08,vastlxy@163.com

%   INPUT:
%       data: 2D data
%       xx: 2D data coordinates
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

gradtmin=para.gradtmin;
rhomin=para.rhomin;
deltamin=para.deltamin;
v_min=para.v_min;
rms=para.rms;

data_filter=imgaussfilt(data,para.sigma);
[size_x,size_y,size_z] = size(data_filter);
rho=data_filter(:);

% [rho_sorted,ordrho]% rho_sorted是排序的结果% ordrhoInd是排序的索引
[rho_sorted,ordrhoInd]=sort(rho,'descend');

% 初始化
maxd=size_x+size_y+size_z;
ND=length(rho);
delta=zeros(1,ND);
IndNearNeigh=zeros(1,ND);
Gradient=zeros(1,ND);

% delta  记录距离，
% IndNearNeigh  记录：两个密度点的联系% index of nearest neighbor with higher density
delta(ordrhoInd(1))=-1;
IndNearNeigh(ordrhoInd(1))=0;
% Gradient(ordrhoInd(1))=0;

% 计算 delta, Gradient
for ii=2:ND
    %    密度降序排序后，即密度第ii大的索引（在rho中）
    ordrho_ii=ordrhoInd(ii);
    % 当前点的密度及坐标
    rho_ii =  rho_sorted(ii);
    if rho_ii >=rms
        delta(ordrho_ii)=maxd;
        delta_ii_xy = xx(ordrho_ii,:);
        %         得到以delta_ii_xy这一点为中心的 k*k 的邻域内的点在原图中的坐标
        bt = kc_coord_2d(delta_ii_xy,size_y,size_x);
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
                dist_i_j=distx(ordrho_ii,ordrhoInd(jj),xx);
                gradient = (rho_jj - rho_ii) / dist_i_j;
                if dist_i_j<=delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    IndNearNeigh(ordrho_ii)=ordrhoInd(jj);
                end
            end
        end
    else
        % rho_ii < para.rms则将点标记为一个特殊值
        IndNearNeigh(ordrho_ii)=ND+1;
    end
end
delta(ordrhoInd(1))=max(delta(:));
disp('delata, rho and Gradient are ok')

% gama = rho'.*delta;

% [gama_sorted, ~] = sort(gama,'descend');
% 
% figure
% plot(1:300,gama_sorted,'.')

% 根据密度和距离来确定类中心
NCLUST=0;
clustInd=-1*ones(1,ND+1);

% 将满足密度，距离最低值的点选出来
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

if is_plot==1
    %     figure
    deltamin_ = deltamin / max(delta(:));
    delta = delta./max(delta(:));
    plot(rho(:),delta(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
%     title ('Decision Graph','FontSize',15.0)
%     xlabel ('\rho')
%     ylabel ('\delta')
    xlim([0.6,max(rho(:))])
    for i=1:NCLUST
        hold on
        %    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',6,'MarkerEdgeColor','m');
    end
    
    hold on
    plot([0, max(rho(:))], [deltamin_, deltamin_], 'r--')
    hold on
    plot([rhomin, rhomin], [0, 1], 'r--')
    gtext('\rho_0','Color','red','FontSize',15)
    gtext('\delta_0','Color','red','FontSize',15)
end

% assignation  将其他非类中心分配到离它最近的类中心中去
% clustInd=-1 表示该点不是类的中心点，属于其他点，等待被分配到某个类中去
% 类的中心点的梯度Gradient被指定为-1
for i=1:ND
    ordrho_i=ordrhoInd(i);
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

for i = 1:size(centInd,1)
    
    rho_clust_i = zeros(size(rho));
    
    index_clust_i = find(clustInd==centInd(i,2));% 得到第i类所有点的索引
    
    index_cc =  intersect(mask_grad,index_clust_i);% 得到第i类梯度大于gradtmin的点的索引
    rho_clust_i(index_clust_i) = rho(index_clust_i); % 得到第i类所有点的密度
    
    rho_cc_mean = mean(rho(index_cc)); % 得到第i类梯度大于gradtmin的点的密度的均值
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % 得到第i类密度大于该阈值的点的索引
    
    index_clust_rho = union(index_cc, index_cc_rho); %得到两个条件都满足的点的索引
    
    clustInd_re(index_clust_rho) = centInd(i,2);
    
    cl_i_point=xx(index_clust_rho,:); % 得到第i类梯度大于gradtmin的点的坐标
    %     ini = find(clustInd==1);
    %     cl_i_point=xx(ini,:); % 得到第i类梯度大于gradtmin的点的坐标
    %     mask_clust=zeros(size_x,size_y,size_z);
    mask_out=zeros(size_x,size_y,size_z);
    for j =1:size(cl_i_point,1)
        %         mask_clust(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=i;
        mask_out(cl_i_point(j,1), cl_i_point(j,2))=1;
    end
    %     对每个类的掩模做形态学处理：闭运算，填补空洞，保留面积最大连通域
    bw = imbinarize(mask_out,0.5);
    L21 = logical(bw);  % 默认8连通； bwlabe(bw,4);
    
    
    BW2 = imfill(L21,'holes'); % 填补空洞
    
    se=strel('disk',1);
    %     BW2=imclose(BW2,se); % 闭运算
    
    L = bwlabel(BW2);   % 对连通区域进行标记
    STATS = regionprops(BW2,'All');
    % 统计上一步标记图像中的连通域的面积分布
    Ar = cat(1, STATS.Area);
    ind = find(Ar ==max(Ar));%找到最大连通区域的标号
%     BW2(L~=ind)=0;%将其它区域置为0
    
    clump_Image = double(BW2);
    mask_clust = clump_Image .* i;
    mask = mask + mask_clust;
    out = out + clump_Image.*data;
end

NCLUST_ = size(centInd,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);