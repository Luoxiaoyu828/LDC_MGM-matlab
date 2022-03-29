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
%       mask:  1,2,3,��

% 2021/04/22  �޸�
%ȥ������С�ܶȡ�������ж�׼�򡣸�Ϊ�Ⱦ��࣬��ȷ������߽磬���������Ա���ж��Ƿ�����
% �޸ģ���������˲����ڶԱ߽�ȷ����ɺ����


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
%       ���˲����޸�,2020/12/08,vastlxy@163.com

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

% [rho_sorted,rho_Ind]% rho_sorted������Ľ��% rho_Ind�����������
[rho_sorted,rho_Ind]=sort(rho,'descend');

% ��ʼ��
maxd=size_x+size_y+size_z;
ND=length(rho);
delta=zeros(1,ND);
IndNearNeigh=zeros(1,ND);
Gradient=zeros(1,ND);

r_region = 4;

% delta  ��¼���룬
% IndNearNeigh  ��¼�������ܶȵ����ϵ% index of nearest neighbor with higher density
delta(rho_Ind(1))=sqrt(size_x^2+size_y^2+size_z^2);   % 2021/01/18 �����޸�  ԭ����-1
IndNearNeigh(rho_Ind(1))=ND+2;    % 2021/01/18 �����޸�  ԭ����0
% Gradient(ordrhoInd(1))=0;

% ���� delta, Gradient
tic
for ii=2:ND
    %�ܶȽ�������󣬼��ܶȵ�ii�����������rho�У�
    ordrho_ii=rho_Ind(ii);
    % �ܶȵ�ii�����ܶ�
    rho_ii =  rho_sorted(ii);
    if rho_ii >= rms
        % ���ȶ�delta��һ���ϴ�ĳ�ʼֵ
        delta(ordrho_ii)=maxd;
        %�ܶȵ�ii��ĵ��Ӧ������
        delta_ii_xy = xx(ordrho_ii,:);
        % �õ���delta_ii_xy��һ��Ϊ���ĵ� k*k*k �������ڵĵ���ԭͼ�е�����
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
        if delta(ordrho_ii)==maxd% �����������õ�"��"��û���ҵ��ȸõ�ǿ��ֵ�ߵĵ㣬�����ȫ������
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
        % rho_ii < para.rms�򽫵���Ϊһ������ֵ
        IndNearNeigh(ordrho_ii)=ND+1;
    end
end

[detla_sorted,~] = sort(delta,'descend');
delta(rho_Ind(1))=detla_sorted(2);

toc
disp('delata, rho and Gradient are ok')

NCLUST=0;
clustInd=-1*ones(1,ND+1);

% ���ֲ�����ֵ��ѡ�������ֲ�����ֵ�㣺���趨��"��"��δ�ҵ��ȸõ�ǿ��ֵ��ĵ�
clust_index = find(IndNearNeigh==ND+2);

clust_num = length(clust_index);
str_1 = sprintf('find local maximun points number: %d\n',clust_num);
fprintf(str_1)
% icl��������¼��i����������xx�е�����ֵ
icl = zeros(clust_num,1);
for ii=1:clust_num
    i = clust_index(ii);  % i ��ʾ���ܶȵ����������������ﶼ��Ϊʶ��Ϊ�������ĵĵ�
    NCLUST=NCLUST+1;
    clustInd(i)=NCLUST;
    icl(NCLUST)=i;  % icl��������¼��i����������xx�е�����ֵ
end
% assignation  �������������ķ��䵽�����������������ȥ
% clustInd=-1 ��ʾ�õ㲻��������ĵ㣬���������㣬�ȴ������䵽ĳ������ȥ
% ������ĵ���ݶ�Gradient��ָ��Ϊ-1
for i=1:ND
    ordrho_i=rho_Ind(i);
    if clustInd(ordrho_i)==-1 % not centroid
        clustInd(ordrho_i)=clustInd(IndNearNeigh(ordrho_i));
    else
        Gradient(ordrho_i)=-1;
    end
end

% ����ÿ����������
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
% centInd [�����ĵ���xx�����µ�����ֵ�� ��������centInd������ֵ: ���������]

centInd = [centInd, centNum];
NCLUST_ = length(centNum);
str_2 = sprintf('Number of cluster(volume > %d): %d\n',v_min, NCLUST_);
fprintf(str_2)

clustInd_re = -1 * ones(size(clustInd));  % �������ȷ���������ƺ˵���������
mask = zeros(size_x,size_y,size_z);
out = zeros(size_x,size_y,size_z);

% �����ݶȽ���ı߽��һ��ѹ�����õ��ƺ���Ҫ����
mask_grad = find(Gradient>gradtmin);
clust_id = 1;
for i = 1:size(centInd,1)
    
    rho_clust_i = zeros(size(rho));
    index_clust_i = find(clustInd==centInd(i,2));% �õ���i�����е������
    
    index_cc =  intersect(mask_grad,index_clust_i);% �õ���i���ݶȴ���gradtmin�ĵ������
    rho_clust_i(index_clust_i) = rho(index_clust_i); % �õ���i�����е���ܶ�
    
    rho_cc_mean = mean(rho(index_cc)); % �õ���i���ݶȴ���gradtmin�ĵ���ܶȵľ�ֵ
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % �õ���i���ܶȴ��ڸ���ֵ�ĵ������
    
    index_clust_rho = union(index_cc, index_cc_rho); %�õ���������������ĵ������
    
    clustInd_re(index_clust_rho) = centInd(i,2);
    if get_out_mask
        cl_i_point=xx(index_clust_rho,:); % �õ���i���ݶȴ���gradtmin�ĵ������
        %     ini = find(clustInd==1);
        %     cl_i_point=xx(ini,:); % �õ���i���ݶȴ���gradtmin�ĵ������
        clump_Image=zeros(size_x,size_y,size_z);
        mask_out=zeros(size_x,size_y,size_z);
        if size(cl_i_point,1)>0
            for j =1:size(cl_i_point,1)
                %         mask_clust(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=i;
                mask_out(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=1;
            end
            %     ��ÿ�������ģ����̬ѧ���������㣬��ն���������������ͨ��
            bw = imbinarize(mask_out,0.5);
            BW2 = logical(bw);           % Ĭ��8��ͨ�� bwlabe(bw,4);
            %     se=strel('disk',3);
            %             SE = strel('sphere', 2);
            %             BW2=imclose(BW2,SE);         % ������
            BW2 = imfill(BW2,'holes');   % ��ն�
            
            STATS = regionprops3(BW2,'All');
            % ͳ����һ�����ͼ���е���ͨ�������ֲ�
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
%       ���˲����޸�,2020/12/08,vastlxy@163.com
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

[rho_sorted,rho_Ind]=sort(rho,'descend');%  rho_sorted������Ľ��% ordrhoInd�����������

% ��ʼ��
maxd = sqrt(size_x^2+size_y^2);
ND=length(rho);
delta = zeros(1,ND);
IndNearNeigh = zeros(1,ND);
Gradient = zeros(1,ND);

% delta  ��¼���룬
% IndNearNeigh  ��¼�������ܶȵ����ϵ% index of nearest neighbor with higher density
delta(rho_Ind(1))=-1;
IndNearNeigh(rho_Ind(1))=ND+2;

% ���� delta, Gradient
for ii=2:ND
    ordrho_ii=rho_Ind(ii);% �ܶȽ�������󣬼��ܶȵ�ii�����������rho�У�
    rho_ii =  rho_sorted(ii);% ��ǰ����ܶȼ�����
    if rho_ii >=rms
        delta(ordrho_ii)=maxd;
        delta_ii_xy = xx(ordrho_ii,:);
        % �õ���delta_ii_xy��һ��Ϊ���ĵ� k*k �������ڵĵ���ԭͼ�е�����
        bt = kc_coord_2d(delta_ii_xy, size_x, size_y, region);
        for j_=1:size(bt,1)
            rho_jj= data_filter(bt(j_,1),bt(j_,2));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            gradient = (rho_jj - rho_ii) / dist_i_j;
            if dist_i_j <= delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                delta(ordrho_ii)= dist_i_j;
                Gradient(ordrho_ii)= gradient;
                % indNearNeigh(ordRho(i)) = ordRho(j); �����ʽ���Ҷ�ʵ�ֵ��ǵõ�ordRho(j)
                IndNearNeigh(ordrho_ii)=(bt(j_,2)-1)*size_x + bt(j_,1);
            end
        end
        if delta(ordrho_ii)==maxd
            % ��������3*3*3��������û���ҵ��ȸõ�ߣ���������ĵ㣬�����ȫ������
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
        IndNearNeigh(ordrho_ii)=ND+1;% rho_ii < para.rms�򽫵���Ϊһ������ֵ
    end
end

[detla_sorted,~] = sort(delta,'descend');
delta(rho_Ind(1))=detla_sorted(2);
disp('delata, rho and Gradient are ok')

NCLUST=0;
clustInd=-1*ones(1,ND+1);

% ���ֲ�����ֵ��ѡ����.�ֲ�����ֵ�㣺���趨��"��"��δ�ҵ��ȸõ�ǿ��ֵ��ĵ�
clust_index = find(IndNearNeigh==ND+2);

clust_num = length(clust_index);
str_1 = sprintf('find local maximun points number: %d\n',clust_num);
fprintf(str_1)

icl = zeros(clust_num,1);% icl��������¼��i����������xx�е�����ֵ
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

% assignation  �������������ķ��䵽�����������������ȥ
% clustInd=-1 ��ʾ�õ㲻��������ĵ㣬���������㣬�ȴ������䵽ĳ������ȥ
% ������ĵ���ݶ�Gradient��ָ��Ϊ-1
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
% centInd [�����ĵ���xx�����µ�����ֵcentInd1�� ��������centInd������ֵcentNum: ���������]
centInd = [centInd1, centNum];
NCLUST_ = length(centNum);
str_2 = sprintf('Number of cluster(volume > %d): %d\n',v_min, NCLUST_);
fprintf(str_2)
clustInd_re = -1 * ones(size(clustInd));  % �������ȷ���������ƺ˵���������
mask = zeros(size_x,size_y,size_z);
out = zeros(size_x,size_y,size_z);

% �����ݶȽ���ı߽��һ��ѹ�����õ��ƺ���Ҫ����
mask_grad = find(Gradient>gradtmin);

for i = 1:size(centInd,1)
    rho_clust_i = zeros(size(rho));
    index_clust_i = find(clustInd==centInd(i,2));% �õ���i�����е������
    index_cc =  intersect(mask_grad,index_clust_i);% �õ���i���ݶȴ���gradtmin�ĵ������
    rho_clust_i(index_clust_i) = rho(index_clust_i); % �õ���i�����е���ܶ�
    
    rho_cc_mean = mean(rho(index_cc)); % �õ���i���ݶȴ���gradtmin�ĵ���ܶȵľ�ֵ
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % �õ���i���ܶȴ��ڸ���ֵ�ĵ������
    
    index_clust_rho = union(index_cc, index_cc_rho); %�õ���������������ĵ������
    
    clustInd_re(index_clust_rho) = centInd(i,2); % ���±�������
end

for i = 1:size(clustInd_re,1)
           
    cl_i_point=xx(index_clust_rho,:); % �õ���i���ݶȴ���gradtmin�ĵ������
    if get_out_mask
        mask_out=zeros(size_x,size_y,size_z);
        for j =1:size(cl_i_point,1)
            %         mask_clust(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=i;
            mask_out(cl_i_point(j,1), cl_i_point(j,2))=1;
        end
        %     ��ÿ�������ģ����̬ѧ���������㣬��ն���������������ͨ��
        bw = imbinarize(mask_out,0.5);
        L21 = logical(bw);  % Ĭ��8��ͨ�� bwlabe(bw,4);
        
        BW2 = imfill(L21,'holes'); % ��ն�
        
        L = bwlabel(BW2);   % ����ͨ������б��
        STATS = regionprops(BW2,'All');
        % ͳ����һ�����ͼ���е���ͨ�������ֲ�
        Area = cat(1, STATS.Area);
        ind = find(Area ==max(Area));%�ҵ������ͨ����ı��
        BW2(L~=ind)=0;%������������Ϊ0
        
        clump_Image = double(BW2);
        mask_clust = clump_Image .* i;
        mask = mask + mask_clust;
        out = out + clump_Image.*data;
    end
end
end
