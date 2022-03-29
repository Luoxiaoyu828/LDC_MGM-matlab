function [NCLUST_,centInd_re,clustInd_re,cluster_info,mask,out,Gradient]=densityCluster_3d(data,xx,para,is_plot,get_out_mask)
%% densityCluster 3D
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

% ������ʼ��
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

% [rho_sorted,rho_Ind]% rho_sorted������Ľ��% rho_Ind�����������
[rho_sorted,rho_Ind] = sort(rho,'descend');

% ��ʼ��
maxd = size_x + size_y + size_z;
ND = length(rho);
delta = zeros(1,ND);
IndNearNeigh = zeros(1,ND);
Gradient = zeros(1,ND);

% delta  ��¼���룬
% IndNearNeigh  ��¼�������ܶȵ����ϵ% index of nearest neighbor with higher density
delta(rho_Ind(1)) = sqrt(size_x^2+size_y^2+size_z^2);   % 2021/01/18 �����޸�  ԭ����-1
IndNearNeigh(rho_Ind(1)) = rho_Ind(1);    % 2021/01/18 �����޸�  ԭ����0

% ���� delta, Gradient
tic
for ii = 2:ND
    %�ܶȽ�������󣬼��ܶȵ�ii�����������rho�У�
    ordrho_ii = rho_Ind(ii);
    rho_ii = rho_sorted(ii);% ��ii���ܶ�
    if rho_ii >= rms
        % ���ȶ�delta��һ���ϴ�ĳ�ʼֵ
        delta(ordrho_ii) = maxd;
        %�ܶȵ�ii��ĵ��Ӧ������
        delta_ii_xy = xx(ordrho_ii,:);
        % �õ���delta_ii_xy��һ��Ϊ���ĵ� (2*k+1)*(2*k+1)*(2*k+1) �������ڵĵ���ԭͼ�е�����
        k = 2;
        bt = kc_coord_3d(delta_ii_xy,size_z,size_y,size_x,k);
        for j_ = 1:size(bt,1)
            rho_jj = data_filter(bt(j_,1),bt(j_,2),bt(j_,3));
            dist_i_j = dist_xyz(delta_ii_xy, bt(j_,:));
            if dist_i_j == 0
                continue
            end
            gradient = (rho_jj - rho_ii) / dist_i_j;
            % ͨ��������ݶ�ֵ�Ŀ���  ����ʵ���ҵ����ȵ�ǰ��ǿ�ȴ�ĵ㡱������ĵ�
            if dist_i_j <= delta(ordrho_ii) && gradient >= 0
                delta(ordrho_ii)= dist_i_j;
                Gradient(ordrho_ii)= gradient;
                % indNearNeigh(ordRho(i)) = ordRho(j); �����ʽ���Ҷ�ʵ�ֵ��ǵõ�ordRho(j)
                IndNearNeigh(ordrho_ii)=(bt(j_,3)-1)*size_x*size_y + (bt(j_,2)-1)*size_x + bt(j_,1);
            end
        end
        % ��������(2*k+1)*(2*k+1)*(2*k+1)��������û���ҵ��ȸõ�ߣ���������ĵ㣬�����ȫ������
        if delta(ordrho_ii)==maxd
            % ��Ϊ���ǿ�Ƚ���������  ����jj�ĵ����Ǳ�ii�ĵ�ǿ�ȴ󣬹���������ж��в����ݶ��ж�
            for jj=1:ii-1
                rho_jj = rho_sorted(jj);
                dist_i_j=distx(ordrho_ii,rho_Ind(jj),xx);  % ���������ľ���
                gradient = (rho_jj - rho_ii) / dist_i_j;    % �����������ݶ�
                if dist_i_j<=delta(ordrho_ii)
                    delta(ordrho_ii)=dist_i_j;
                    Gradient(ordrho_ii)= gradient;
                    
                    % ���������ľ������4,����Ϊ���ĵ��ѡ��
                    if dist_i_j <= 4
                        IndNearNeigh(ordrho_ii)=rho_Ind(jj);
                    else
                        IndNearNeigh(ordrho_ii)=ND+1;
                    end
                    
                end
            end
        end
    else
        % rho_ii < para.rms�򽫵���Ϊһ������ֵ
        IndNearNeigh(ordrho_ii)=ND+1;
    end
    
end

[detla_sorted,~] = sort(delta,'descend');
% [rho_sorted,~] = sort(rho,'descend');
% delta(rho_Ind(1))=max(delta(:));
delta(rho_Ind(1))=detla_sorted(2);

toc
disp('delata, rho and Gradient are ok')

% �����ܶȺ;�����ȷ��������
NCLUST=0;
clustInd=-1*ones(1,ND+1);

disp('find cluster center point by delata and rho')

% ѡȡ�����ĵ㣺�������ܶȣ��������ֵ�ĵ�ѡ������Ϊ�����ĵ�
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

% ���������ȥ��һ���ֺ� clust_idΪ���µı��  centInd_reΪ�µ��ƺ����ĵ�����
clust_id = 1;
centInd_re = [];

for i = 1:size(centInd,1)
    
    rho_clust_i = zeros(size(rho));
    index_clust_i = find(clustInd==centInd(i,2));% �õ���i�����е������
    
    index_cc =  intersect(mask_grad,index_clust_i);% �õ���i���ݶȴ���gradtmin�ĵ������
    rho_clust_i(index_clust_i) = rho(index_clust_i); % �õ���i�����е���ܶ�
    
    rho_cc_mean = 0.2 * mean(rho(index_cc)); % �õ���i���ݶȴ���gradtmin�ĵ���ܶȵľ�ֵ��0.2��
    index_cc_rho = find(rho_clust_i>rho_cc_mean); % �õ���i���ܶȴ��ڸ���ֵ�ĵ������
    
    index_clust_rho = union(index_cc, index_cc_rho); %�õ���������������ĵ������
    
    % �ж�ͬʱ���������������ƺ����յ�����Ƿ��㹻��
    if length(index_clust_rho) > v_min
        clustInd_re(index_clust_rho) = clust_id; % �Ա������ƺ����±��
        centInd_re = [centInd_re; [centInd(i,1), clust_id]];
        if get_out_mask
            cl_i_point=xx(index_clust_rho,:); %�õ���i���ݶȴ���gradtmin��ǿ�ȴ���ԭʼǿ�Ⱦ�ֵ�ĵ������
            clump_Image=zeros(size_x,size_y,size_z);
            
            mask_out=zeros(size_x,size_y,size_z);
            for j =1:size(cl_i_point,1)
                mask_out(cl_i_point(j,1), cl_i_point(j,2), cl_i_point(j,3))=1;
            end
            % ��ÿ�������ģ����̬ѧ������ն���������������ͨ��
            bw = imbinarize(mask_out,0.5);
            BW2 = logical(bw);           % Ĭ��8��ͨ�� bwlabe(bw,4);
            BW2 = imfill(BW2,'holes');   % ��ն�
            
            STATS = regionprops3(BW2,'All'); % ͳ����һ�����ͼ���е���ͨ�������ֲ�
           
            max_ind = find( STATS.Volume==max( STATS.Volume));
            points = STATS.VoxelList{max_ind};

            for ii = 1: STATS.Volume(max_ind)
                clump_Image(points(ii,2),points(ii,1),points(ii,3))=1;
            end
            mask_clust = clump_Image .* clust_id;
            mask = mask + mask_clust; % �õ���ģ
            out = out + clump_Image.*data; % ����ģ�õ�������
        end
        clust_id = clust_id+1; % ���ƺ˵ı�Ž��и���
    end
end

NCLUST_ = size(centInd_re,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);

