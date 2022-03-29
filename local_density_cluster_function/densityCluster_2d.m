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

gradtmin=para.gradtmin;
rhomin=para.rhomin;
deltamin=para.deltamin;
v_min=para.v_min;
rms=para.rms;

data_filter=imgaussfilt(data,para.sigma);
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
delta(ordrhoInd(1))=-1;
IndNearNeigh(ordrhoInd(1))=0;
% Gradient(ordrhoInd(1))=0;

% ���� delta, Gradient
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
disp('delata, rho and Gradient are ok')

% gama = rho'.*delta;

% [gama_sorted, ~] = sort(gama,'descend');
% 
% figure
% plot(1:300,gama_sorted,'.')

% �����ܶȺ;�����ȷ��������
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

% assignation  �������������ķ��䵽�����������������ȥ
% clustInd=-1 ��ʾ�õ㲻��������ĵ㣬���������㣬�ȴ������䵽ĳ������ȥ
% ������ĵ���ݶ�Gradient��ָ��Ϊ-1
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
    
    se=strel('disk',1);
    %     BW2=imclose(BW2,se); % ������
    
    L = bwlabel(BW2);   % ����ͨ������б��
    STATS = regionprops(BW2,'All');
    % ͳ����һ�����ͼ���е���ͨ�������ֲ�
    Ar = cat(1, STATS.Area);
    ind = find(Ar ==max(Ar));%�ҵ������ͨ����ı��
%     BW2(L~=ind)=0;%������������Ϊ0
    
    clump_Image = double(BW2);
    mask_clust = clump_Image .* i;
    mask = mask + mask_clust;
    out = out + clump_Image.*data;
end

NCLUST_ = size(centInd,1);
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST_);