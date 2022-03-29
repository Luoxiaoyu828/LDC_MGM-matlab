function [outcat, out, mask, Gradient] = localDenClust2_1(data,para,is_plot,get_out_mask,figure_name)
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

% 2021/04/11  修改  去掉了rho 和 delata  最低阈值的两个参数

tic
if length(size(data))<=2
    disp('ok')
    [size_x,size_y]=size(data);
    [p_i, p_j] = meshgrid(1:size_y, 1:size_x);
    xx=[p_j(:), p_i(:)];
    [NClust,centInd,clustInd,~,mask,out,Gradient]=densityCluster_2d(data,xx,para,is_plot);
end
if length(size(data))==3
    [size_x,size_y,size_z]=size(data);
    [p_i, p_j, p_k] = meshgrid(1:size_y, 1:size_x, 1:size_z);
    xx=[p_j(:), p_i(:), p_k(:)];
    [NClust,centInd,clustInd,~,mask,out,~]=densityCluster_3d(data,xx,para,is_plot,get_out_mask,figure_name);
end
[outcat] = findclumps(NClust,xx,clustInd,centInd,data);
toc

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
Gravitation=zeros(1,ND);
std_time = 3;
r_region = 3;


% delta  记录距离，
% IndNearNeigh  记录：两个密度点的联系% index of nearest neighbor with higher density
delta(rho_Ind(1))=sqrt(size_x^2+size_y^2+size_z^2);   % 2021/01/18 做了修改  原来是-1
IndNearNeigh(rho_Ind(1))=rho_Ind(1);    % 2021/01/18 做了修改  原来是0
% Gradient(ordrhoInd(1))=0;

% 计算 delta, Gradient
tic
rho_temp = rho;
% 利用3 sigma  原则排除异常点
% rho_temp = find_mean_std(rho_temp);

rho_th = mean(rho_temp)+std_time*std(rho_temp);
% rho_th = mean(rho)+std_time*std(rho);
disp('***************最低的密度值')
disp(rho_th)

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
                % indNearNeigh(ordRho(i)) = ordRho(j); 下面等式的右端实现的是得到ordRho(j)
%                 if dist_i_j<=1.73*3
                    IndNearNeigh(ordrho_ii)=(bt(j_,3)-1)*size_x*size_y + (bt(j_,2)-1)*size_x + bt(j_,1);
%                 else
%                     IndNearNeigh(ordrho_ii)=ND+1;
%                 end
%                                 Gravitation(ordrho_ii) = rho_ii^2 / dist_i_j;
                e_ind = 1.*(rho_ii > rho_th)-1.*(rho_ii < rho_th & rho_ii == rho_th);
                Gravitation(ordrho_ii) = (delta(ordrho_ii)^e_ind)^2;
            end
        end
        if delta(ordrho_ii)==maxd
            % 表明，在3*3*3的领域中没有找到比该点高，距离最近的点，则进行全局搜索
            for jj=1:ii-1
                rho_jj = rho_sorted(jj);
                dist_i_j=distx(ordrho_ii,rho_Ind(jj),xx);
                gradient = (rho_jj - rho_ii) / dist_i_j;
                if dist_i_j<=delta(ordrho_ii) 
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    
                    %                     if dist_i_j<=1.73*2
                    %                         IndNearNeigh(ordrho_ii)=rho_Ind(jj);
                    %                     else
                    IndNearNeigh(ordrho_ii)=ND+1;
                    %                     end
                    e_ind = 1.*(rho_ii > rho_th)-1.*(rho_ii < rho_th & rho_ii ==rho_th);
                    Gravitation(ordrho_ii) = (delta(ordrho_ii)^e_ind)^2;
                    %                     Gravitation(ordrho_ii) = rho_ii^2 / dist_i_j;
                end
            end
        end
    else
        % rho_ii < para.rms则将点标记为一个特殊值
        IndNearNeigh(ordrho_ii)=ND+1;
        Gravitation(ordrho_ii) = 1;
    end
end

[detla_sorted,~] = sort(delta,'descend');
% [rho_sorted,~] = sort(rho,'descend');
% delta(rho_Ind(1))=max(delta(:));
delta(rho_Ind(1))=detla_sorted(2);

Gravitation(rho_Ind(1)) = detla_sorted(2)^2;

toc
disp('delata, rho and Gradient are ok')

% gama = rho'.*delta;
% 
% gama_temp = gama;
% gama_temp = find_mean_std(gama_temp);

% disp([length(gama_temp)])

% [gama_sorted, ~] = sort(gama,'descend');

% figure
% plot(1:300,gama_sorted(1:300),'.')

% 根据密度和距离来确定类中心
NCLUST=0;
clustInd=-1*ones(1,ND+1);
grav_th = mean(Gravitation)+std_time*std(Gravitation);
grav_th = log10(grav_th);
Gravitation(:) = log10(Gravitation(:));


% 将满足密度，距离最低值的点选出来
% clust_index = find(and(rho>rhomin,delta'>deltamin)>0);
% clust_index = find(and(rho>rho_th_1,Gravitation'>grav_th)>0);


rho_tt_1 = rho(delta>-0.00001 & delta<1);
rho_tf_1 = rho(delta>1.01);

% rho_tt_1 = find_mean_std(rho_tt_1);
% rho_tf_1 = find_mean_std(rho_tf_1);
rho_th_1 = mean(rho_tt_1)+std_time*std(rho_tt_1);
rho_th_2 = mean(rho_tf_1)+std_time*std(rho_tf_1);

figure
subplot(1,2,1)
hist(rho_tf_1,2000)
xlim([0,5])
subplot(1,2,2)
hist(rho_tt_1,2000)
xlim([-2,5])

disp('*****************Gravitation')
disp([mean(Gravitation),std(Gravitation),grav_th])

disp('************第二个密度阈值')
disp('框内点   噪声+峰值   全局')
disp([rho_th_1,rho_th_2,rho_th])

if rho_th_1<rho_th_2
    rho_th = rho_th_1;
else
    rho_th = rho_th_2;
end
% 将满足密度，距离最低值的点选出来
% clust_index = find(and(rho>rhomin,delta'>deltamin)>0);
clust_index = find(and(rho>rho_th,Gravitation'>grav_th)>0);

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
    set(gcf,'position',[50 50 1200 500])
    %  定义窗口到屏幕左边的距离是150，到屏幕下方的距离是150，图片width = 150 ,height = 150(向量第３个是width)
    subplot(1,4,[1 2])
    plot(rho(:),delta(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
    title ('Decision Graph','FontSize',15.0)
    xlabel ('\rho')
    ylabel ('\delta')
    for i=1:NCLUST
        hold on
        %    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',6);
    end
%     figure
    subplot(1,4,[3 4])
    plot(rho(:),Gravitation(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
    for i=1:NCLUST
        hold on
        %    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        plot(rho(icl(i)),Gravitation(icl(i)),'o','MarkerSize',6);
    end
    hold on
    aa_grad = floor(max(Gravitation));
    aa_rho = floor(max(rho));
    plot([0,aa_rho],[grav_th,grav_th],'-')
    hold on
    plot([rho_th,rho_th],[0, aa_grad,],'-')
    hold on
    plot([rho_th_2,rho_th_2],[0, aa_grad,],'--')
    saveas(gcf,figure_name)
    
end

% assignation  将其他非类中心分配到离它最近的类中心中去
% clustInd=-1 表示该点不是类的中心点，属于其他点，等待被分配到某个类中去
% 类的中心点的梯度Gradient被指定为-1
for i=1:ND
    ordrho_i=rho_Ind(i);
    if clustInd(ordrho_i)==-1 % not centroid
        %         disp(ordrho_i)
        %         disp(IndNearNeigh(ordrho_i))
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

NCLUST_ = size(centInd,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);


%% densityCluster 2D
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

data_filter=imgaussfilt3(data,para.sigma);
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

gama = rho'.*delta;

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
    delta = delta./max(delta(:));
    plot(rho(:),delta(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
    title ('Decision Graph','FontSize',15.0)
    xlabel ('\rho')
    ylabel ('\delta')
    xlim([0.6,max(rho(:))])
    for i=1:NCLUST
        hold on
        %    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',6);
    end
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
    BW2(L~=ind)=0;%将其它区域置为0
    
    clump_Image = double(BW2);
    mask_clust = clump_Image .* i;
    mask = mask + mask_clust;
    out = out + clump_Image.*data;
end

NCLUST_ = size(centInd,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);

%%
function index_value = kc_coord_3d(delta_ii_xy,xm,ym,zm,r_region)
jt = delta_ii_xy(1);
it = delta_ii_xy(2);
kt = delta_ii_xy(3);
r=r_region;
% bt=a(max(1,kt-r):min(zm,kt+r), max(1,jt-r):min(ym,jt+r), max(1,it-r):min(xm,it+r));
% bt1=a(max(1,kt-r):min(zm,kt+r), max(1,it-r):min(xm,it+r), max(1,jt-r):min(ym,jt+r));
% bt2=a(max(1,it-r):min(xm,it+r), max(1,kt-r):min(zm,kt+r), max(1,jt-r):min(ym,jt+r) );
% bt3=a(max(1,it-r):min(xm,it+r), max(1,jt-r):min(ym,jt+r), max(1,kt-r):min(zm,kt+r) );
% bt4=a(max(1,jt-r):min(ym,jt+r), max(1,kt-r):min(zm,kt+r),max(1,it-r):min(xm,it+r) );
% bt5=a(max(1,jt-r):min(ym,jt+r),max(1,kt-r):min(zm,kt+r),max(1,it-r):min(xm,it+r) );

[p_i,p_j,p_k] = meshgrid(max(1,kt-r):min(xm,kt+r), max(1,it-r):min(ym,it+r), max(1,jt-r):min(zm,jt+r));
index_value=[p_k(:),p_j(:),p_i(:)];

function index_value = kc_coord_2d(delta_ii_xy,xm,ym)
jt = delta_ii_xy(1);
it = delta_ii_xy(2);

r=3;
% bt=a(max(1,kt-r):min(zm,kt+r), max(1,jt-r):min(ym,jt+r), max(1,it-r):min(xm,it+r));
% bt1=a(max(1,kt-r):min(zm,kt+r), max(1,it-r):min(xm,it+r), max(1,jt-r):min(ym,jt+r));
% bt2=a(max(1,it-r):min(xm,it+r), max(1,kt-r):min(zm,kt+r), max(1,jt-r):min(ym,jt+r) );
% bt3=a(max(1,it-r):min(xm,it+r), max(1,jt-r):min(ym,jt+r), max(1,kt-r):min(zm,kt+r) );
% bt4=a(max(1,jt-r):min(ym,jt+r), max(1,kt-r):min(zm,kt+r),max(1,it-r):min(xm,it+r) );
% bt5=a(max(1,jt-r):min(ym,jt+r),max(1,kt-r):min(zm,kt+r),max(1,it-r):min(xm,it+r) );

[p_i,p_j] = meshgrid(max(1,it-r):min(ym,it+r), max(1,jt-r):min(xm,jt+r));
index_value=[p_j(:),p_i(:)];

%%
function distance = dist_xyz(point_a, point_b)
temp=point_a-point_b;
% distance=sqrt(temp(1).^2+temp(2).^2+temp(3).^2);
distance = sqrt(sum(temp.^2));

function distance = dist_xyz_2d(point_a, point_b)
temp=point_a-point_b;
distance=sqrt(temp(1).^2+temp(2).^2);

%%
function [distance] = distx(kc1,kc2,xx)
% distance=abs(xx(kc1,1)-xx(kc2,1))+abs(xx(kc1,2)-xx(kc2,2))+abs(xx(kc1,3)-xx(kc2,3));
if size(xx,2)==3
    distance=sqrt((xx(kc1,1)-xx(kc2,1)).^2+(xx(kc1,2)-xx(kc2,2)).^2+(xx(kc1,3)-xx(kc2,3)).^2);
else
    distance=sqrt((xx(kc1,1)-xx(kc2,1)).^2+(xx(kc1,2)-xx(kc2,2)).^2);
end

%% calculate the outcat
function [outcat] = findclumps(NCLUST,xx,clustInd,centInd,data)
% centInd 代表聚类中心点的索引
dim = size(xx,2);
if dim==3
    clustSum=zeros(1,NCLUST);
    clustVolume=zeros(1,NCLUST);
    clustPeak=zeros(NCLUST,1);
    clump_Cen = zeros(NCLUST,dim);
    clustSize = zeros(NCLUST,dim);
    [size_x,size_y,size_z]=size(data);
    
    clump_Peak=xx(centInd(:,1),:);
    precent=1;
    cl_result = zeros(size(data));
    for i =1:NCLUST
        cl_1=zeros(size_x,size_y,size_z);
        cl_1_index_=xx((clustInd==centInd(i,2)),:);
        %     cl_1_index_=xx((clustInd==centInd(i,2)),:);
        
        clustNum=size(cl_1_index_,1);
        for j =1:clustNum
            cl_1(cl_1_index_(j,1), cl_1_index_(j,2), cl_1_index_(j,3))=1;
        end
        CC = bwconncomp(cl_1);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        big_clust_idx = CC.PixelIdxList{idx};
        
        cl_1_index=xx(big_clust_idx,:);
        cl_1=zeros(size_x,size_y,size_z);
        clustNum=size(cl_1_index,1);
        %     cl_1_index=xx((clustInd==index_cl(i)),:);
        clump_sum_=zeros(clustNum,1);
        for j =1:clustNum
            cl_1(cl_1_index(j,1), cl_1_index(j,2), cl_1_index(j,3))=1;
            clump_sum_(j,1)=data(cl_1_index(j,1), cl_1_index(j,2), cl_1_index(j,3));
        end
        
        clustsum=sum(clump_sum_);
        clump_Cen(i,:)= clump_sum_'*cl_1_index(:,:) / clustsum;
        clustVolume(i)=clustNum;
        clustSum(i)=clustsum;
        %     peak是在滤波之后的数据得到的
        %     clustPeak(i)=data_filter(cluster_points(i,1),cluster_points(i,2),cluster_points(i,3));
        clustPeak(i)=data(clump_Peak(i,1),clump_Peak(i,2),clump_Peak(i,3));
        x_i=cl_1_index-clump_Cen(i,:);
        clustSize(i,:)=sqrt(((clump_sum_'*x_i(:,:).^2) / clustsum)-(clump_sum_'*x_i(:,:)./clustsum).^2);
        
        temp=cl_1.*data;
        temp_sort=sort(temp(:),'descend');
        
        inde_end = round(precent*size(cl_1_index,1));
        temp_mean=mean(temp_sort(inde_end));
        temp(temp<temp_mean)=0;
        temp(temp>=temp_mean)=i;
        cl_result=cl_result+temp;
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
        cl_1=zeros(size_x,size_y);
        cl_1_index_=xx((clustInd==centInd(i,2)),:);
        %     cl_1_index_=xx((clustInd==centInd(i,2)),:);
        
        clustNum=size(cl_1_index_,1);
        for j =1:clustNum
            cl_1(cl_1_index_(j,1), cl_1_index_(j,2))=1;
        end
        CC = bwconncomp(cl_1);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        big_clust_idx = CC.PixelIdxList{idx};
        
        cl_1_index=xx(big_clust_idx,:);
        cl_1=zeros(size_x,size_y);
        clustNum=size(cl_1_index,1);
        %     cl_1_index=xx((clustInd==index_cl(i)),:);
        clump_sum_=zeros(clustNum,1);
        for j =1:clustNum
            cl_1(cl_1_index(j,1), cl_1_index(j,2))=1;
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
        
        temp=cl_1.*data;
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
outcat=[clump_Peak,clump_Cen,clustSize,clustPeak,clustSum',clustVolume'];


function gama_temp = find_mean_std(gama_temp)
mean_gama_temp = [];
std_gama_temp = [];
i_gama_temp = 1;
mean_gama_temp = [mean_gama_temp mean(gama_temp)];
std_gama_temp = [std_gama_temp, std(gama_temp)];
while 1
   
   index_gama_temp = find(gama_temp>(mean_gama_temp(i_gama_temp)+3*std_gama_temp(i_gama_temp)));
   gama_temp(index_gama_temp)=[];
   mean_gama_temp = [mean_gama_temp mean(gama_temp)];
   std_gama_temp = [std_gama_temp, std(gama_temp)];
%    disp('***********************rho_ite')
%    disp([mean_rho_temp])
%    disp([std_rho_temp])
%    disp([length(rho_temp)])
   if abs(mean_gama_temp(i_gama_temp) - mean_gama_temp(i_gama_temp+1)) < 0.0001
       break
   end
   i_gama_temp = i_gama_temp+1;
end