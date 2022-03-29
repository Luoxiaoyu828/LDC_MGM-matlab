function [outcat, out, mask, Gradient] = localDenClust(data,para,is_plot,get_out_mask)
%   INPUT:
%       data: 3D data or 2D data
%       xx: 3D or 2D data coordinates
%       para: clustering parameters(Default parameters)
%                 para.rhomin: Minimum density. (2.5) 
%                 para.deltamin: Minimum delta. (4)
%                 para.v_min: Minimum volume. (64)
%                 para.rms: The noise level of the data, used for data
%                 truncation calculation. (0.23)
%                 para.sigma: Standard deviation of Gaussian filtering. (0.6)
%                 para.gradtmin: Data point gradient. (0.4);
%       is_plot: 1 denotes that plot the Decision Graph, otherwise 0.(0)
%       get_out_mask: 1 denotes that return out and mask, otherwise 0.(0)
%   OUTPUT:
%       outcat: description of the clump:
%               Peak1,Peak2,Peak3,Cen1,Cen2,Cen3,Size1,Size2,Size3,Sum,Peak,Volume
%       out: data obtained through mask
%       mask:  1,2,3,��


if nargin ==1
    % ����Ĭ�ϲ���
    para.rhomin=2.5;
    para.deltamin=4;
    para.v_min=64;
    para.rms=0.23;
    para.sigma=0.6;
    para.gradtmin=0.4;
    is_plot = 0;
    get_out_mask=0;
end


tic

if length(size(data))<=2
    disp('clustering 2d data')
    [size_x,size_y]=size(data);
    [p_i, p_j] = meshgrid(1:size_y, 1:size_x);
    xx=[p_j(:), p_i(:)];
    [NClust,centInd,clustInd,~,mask,out,Gradient]=densityCluster_2d(data,xx,para,is_plot);
end
if length(size(data))==3
    disp('clustering 3d data')
    [size_x,size_y,size_z]=size(data);
    [p_i, p_j, p_k] = meshgrid(1:size_y, 1:size_x, 1:size_z);
    xx=[p_j(:), p_i(:), p_k(:)];
    [NClust,centInd,clustInd,~,mask,out,~]=densityCluster_3d(data,xx,para,is_plot,get_out_mask);
end
[outcat] = findclumps(NClust,xx,clustInd,centInd,data);

toc

%% densityCluster 3D
function [NCLUST_,centInd,clustInd_re,cluster_info,mask,out,Gradient]=densityCluster_3d(data,xx,para,is_plot,get_out_mask)
%   SEE the following paper published in *SCIENCE* for more details:
%       Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%       Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
%       ���˲����޸�,2020/12/08,vastlxy@163.com

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

%% ��ʼ���������  ***************************************************
gradtmin=para.gradtmin;
rhomin=para.rhomin;
deltamin=para.deltamin;
v_min=para.v_min;
rms=para.rms;

data_filter=imgaussfilt3(data,para.sigma);
[size_x,size_y,size_z] = size(data_filter);
rho=data_filter(:);

% [rho_sorted,ordrho]% rho_sorted������Ľ��% ordrhoInd�����������
[rho_sorted,ordrhoInd]=sort(rho,'descend');

% ��ʼ��
maxd=size_x+size_y+size_z;
ND=length(rho);
delta=zeros(1,ND);
IndNearNeigh=zeros(1,ND);
Gradient=zeros(1,ND);

% delta  ��¼���룬
% IndNearNeigh  ��¼�������ܶȵ����ϵ% index of nearest neighbor with higher density
delta(ordrhoInd(1))=sqrt(size_x^2+size_y^2+size_z^2);   % 2021/01/18 �����޸�  ԭ����-1
IndNearNeigh(ordrhoInd(1))=ordrhoInd(1);    % 2021/01/18 �����޸�  ԭ����0
% Gradient(ordrhoInd(1))=0;

%% ���� delta, Gradient *********************************************
tic
for ii=2:ND
    %�ܶȽ�������󣬼��ܶȵ�ii�����������rho�У�
    ordrho_ii=ordrhoInd(ii);
    % ��ǰ����ܶȼ�����
    rho_ii =  rho_sorted(ii);
    if rho_ii >=rms
        % ���ȶ�delta��һ���ϴ�ĳ�ʼֵ
        delta(ordrho_ii)=maxd;
        %�ܶȵ�ii��ĵ��Ӧ������
        delta_ii_xy = xx(ordrho_ii,:);
        % �õ���delta_ii_xy��һ��Ϊ���ĵ� k*k*k �������ڵĵ���ԭͼ�е�����
        bt = kc_coord_3d(delta_ii_xy,size_z,size_y,size_x);
        for j_=1:size(bt,1)
            rho_jj= data_filter(bt(j_,1),bt(j_,2),bt(j_,3));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            gradient = (rho_jj - rho_ii) / dist_i_j;
            if dist_i_j <= delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                delta(ordrho_ii)= dist_i_j;
                Gradient(ordrho_ii)= gradient;
                % indNearNeigh(ordRho(i)) = ordRho(j); �����ʽ���Ҷ�ʵ�ֵ��ǵõ�ordRho(j)
                IndNearNeigh(ordrho_ii)=(bt(j_,3)-1)*size_x*size_y + (bt(j_,2)-1)*size_x + bt(j_,1);
            end
        end
        if delta(ordrho_ii)==maxd
            % ��������3*3*3��������û���ҵ��ȸõ�ߣ���������ĵ㣬�����ȫ������
            for jj=1:ii-1
                rho_jj = rho_sorted(jj);
                dist_i_j=distx(ordrho_ii,ordrhoInd(jj),xx);
                if dist_i_j==0
                    % ���������0ʱ������i,j����һ��������������ѭ��
                    continue
                end
                gradient = (rho_jj - rho_ii) / dist_i_j;
                if dist_i_j<=delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    IndNearNeigh(ordrho_ii)=ordrhoInd(jj);
                end
            end
        end
    else
        % rho_ii < para.rms�򽫵���Ϊһ������ֵ
        IndNearNeigh(ordrho_ii)=ND+1;
    end
end
delta(ordrhoInd(1))=max(delta(:));
toc
disp('delata, rho and Gradient are ok')

%% �����ܶȺ;�����ȷ�������� *****************************************
NCLUST=0;
clustInd=-1*ones(1,ND+1);

% �������ܶȣ��������ֵ�ĵ�ѡ����
clust_index = find(and(rho>rhomin,delta'>deltamin)>0);
clust_num = length(clust_index);

% icl��������¼��i����������xx�е�����ֵ
icl = zeros(clust_num,1);
for ii=1:clust_num
    i = clust_index(ii);
    NCLUST=NCLUST+1;
    clustInd(i)=NCLUST;
    icl(NCLUST)=i;
end

% the Decision Graph
if is_plot==1
    % ������ͼ���Ƴ���
    make_decision_plot(delta,rho,NCLUST,icl)
end

%% assignation  �������������ķ��䵽�����������������ȥ**********************
% clustInd=-1 ��ʾ�õ㲻��������ĵ㣬���������㣬�ȴ������䵽ĳ������ȥ
% ������ĵ���ݶ�Gradient��ָ��Ϊ-1
for i=1:ND
    ordrho_i=ordrhoInd(i);
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
% centInd [�����ĵ���xx�����µ�����ֵ�� ��������centInd������ֵ: ���������]
centInd = [centInd, centNum];

clustInd_re = -1 * ones(size(clustInd));  % �������ȷ���������ƺ˵���������
mask = zeros(size_x,size_y,size_z);
out = zeros(size_x,size_y,size_z);

%% �����ݶȽ���ı߽��һ��ѹ�����õ��ƺ���Ҫ���� ***********************
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
            clust_id= clust_id+1;
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
%       ���˲����޸�,2020/12/08,vastlxy@163.com

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


%% ��ʼ���������  ***************************************************
gradtmin=para.gradtmin;
rhomin=para.rhomin;
deltamin=para.deltamin;
v_min=para.v_min;
rms=para.rms;

data_filter=imgaussfilt3(data,para.sigma);
[size_x,size_y,size_z] = size(data_filter);
rho=data_filter(:);

% [rho_sorted,ordrho]% rho_sorted������Ľ��% ordrhoInd�����������
[rho_sorted,ordrhoInd]=sort(rho,'descend');

% ��ʼ��
maxd=size_x+size_y+size_z;
ND=length(rho);
delta=zeros(1,ND);
IndNearNeigh=zeros(1,ND);
Gradient=zeros(1,ND);
d_density=zeros(1,ND);

% delta  ��¼���룬
% IndNearNeigh  ��¼�������ܶȵ����ϵ% index of nearest neighbor with higher density
delta(ordrhoInd(1))=maxd;
IndNearNeigh(ordrhoInd(1))=ordrhoInd(1);  %��ǰ�ǽ�0��ֵ�����������޸�ΪordrhoInd(1)
% Gradient(ordrhoInd(1))=0;

%% ���� delta, Gradient *********************************************
tic
for ii=2:ND
    %    �ܶȽ�������󣬼��ܶȵ�ii�����������rho�У�
    ordrho_ii=ordrhoInd(ii);
    % ��ǰ����ܶȼ�����
    rho_ii =  rho_sorted(ii);
    if rho_ii >=rms
        delta(ordrho_ii)=maxd;
        delta_ii_xy = xx(ordrho_ii,:);
        %         �õ���delta_ii_xy��һ��Ϊ���ĵ� k*k �������ڵĵ���ԭͼ�е�����
        bt = kc_coord_2d(delta_ii_xy,size_y,size_x);
        for j_=1:size(bt,1)
            rho_jj= data_filter(bt(j_,1),bt(j_,2));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            if dist_i_j==0
                % ���������0ʱ������i,j����һ��������������ѭ��
                continue
            end
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
                dist_i_j=distx(ordrho_ii,ordrhoInd(jj),xx);
                d_density = rho_jj - rho_ii;
                gradient = d_density / dist_i_j;
                % ������б��������Ի����ܶȣ��������ݶ����жϻᵼ��
                if dist_i_j<=delta(ordrho_ii) && gradient>=Gradient(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    IndNearNeigh(ordrho_ii)=ordrhoInd(jj);
                end
            end
        end
    else
        % rho_ii < para.rms�򽫵���Ϊһ������ֵ
        IndNearNeigh(ordrho_ii)=ND+1;
    end
end
delta(ordrhoInd(1))=max(delta(:));
toc
disp('delata, rho and Gradient are ok')

%% �����ܶȺ;�����ȷ��������******************************************
NCLUST=0;
clustInd=-1*ones(1,ND+1);

% �������ܶȣ��������ֵ�ĵ�ѡ����
clust_index = find(and(rho>rhomin,delta'>deltamin)>0);
clust_num = length(clust_index);

% icl��������¼��i����������xx�е�����ֵ
icl = zeros(clust_num,1);
for ii=1:clust_num
    i = clust_index(ii);
    NCLUST=NCLUST+1;
    clustInd(i)=NCLUST;
    icl(NCLUST)=i;
end

if is_plot==1
    % ������ͼ���Ƴ���
    make_decision_plot(delta,rho,NCLUST,icl)
end

% assignation  �������������ķ��䵽�����������������ȥ
% clustInd=-1 ��ʾ�õ㲻��������ĵ㣬���������㣬�ȴ������䵽ĳ������ȥ
% ������ĵ���ݶ�Gradient��ָ��Ϊ-1
for i=1:ND
    ordrho_i=ordrhoInd(i);
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
% centInd [�����ĵ���xx�����µ�����ֵ�� ��������centInd������ֵ: ���������]
centInd = [centInd, centNum];

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
    
    clustInd_re(index_clust_rho) = centInd(i,2);
    
    cl_i_point=xx(index_clust_rho,:); % �õ���i���ݶȴ���gradtmin�ĵ������
    %     ini = find(clustInd==1);
    %     cl_i_point=xx(ini,:); % �õ���i���ݶȴ���gradtmin�ĵ������
    %     mask_clust=zeros(size_x,size_y,size_z);
    mask_out=zeros(size_x,size_y,size_z);
    for j =1:size(cl_i_point,1)
        %         mask_clust(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=i;
        mask_out(cl_i_point(j,1), cl_i_point(j,2))=1;
    end
    %     ��ÿ�������ģ����̬ѧ���������㣬��ն���������������ͨ��
    bw = imbinarize(mask_out,0.5);
    L21 = logical(bw);  % Ĭ��8��ͨ�� bwlabe(bw,4);
    
    
    BW2 = imfill(L21,'holes'); % ��ն�
    
    %     se=strel('disk',1);
    %     BW2=imclose(BW2,se); % ������,�ᵼ���������������ƺ˵���ģֵ���
    
    L = bwlabel(BW2);   % ����ͨ������б��
    STATS = regionprops(BW2,'All');
    % ͳ����һ�����ͼ���е���ͨ�������ֲ�
    Ar = cat(1, STATS.Area);
    ind = find(Ar ==max(Ar));%�ҵ������ͨ����ı��
    BW2(L~=ind)=0;%������������Ϊ0
    
    clump_Image = double(BW2);
    mask_clust = clump_Image .* i;
    mask = mask + mask_clust;
    out = out + clump_Image.*data;
end

NCLUST_ = size(centInd,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);

%%
function index_value = kc_coord_3d(delta_ii_xy,xm,ym,zm)
jt = delta_ii_xy(1);
it = delta_ii_xy(2);
kt = delta_ii_xy(3);
r=3;
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
distance=sqrt(sum(temp.^2));

%%
function [distance] = distx(kc1,kc2,xx)
% distance=abs(xx(kc1,1)-xx(kc2,1))+abs(xx(kc1,2)-xx(kc2,2))+abs(xx(kc1,3)-xx(kc2,3));
xx1=xx(kc1,:);
xx2=xx(kc2,:);
distance = sqrt(sum((xx1-xx2).^2));

%% calculate the outcat
function [outcat] = findclumps(NCLUST,xx,clustInd,centInd,data)
% centInd ����������ĵ������
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
        %     peak�����˲�֮������ݵõ���
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
        %     peak�����˲�֮������ݵõ���
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

function make_decision_plot(delta,rho,NCLUST,icl)
figure
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
