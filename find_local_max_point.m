function local_max = find_local_max_point(data)
% clear
% close all
% ����׼��
addpath E:\efficient_SVM_edge_detection\model
n=3;    % �˲����ߴ�Ĵ�С��(2*n+1)*(2*n+1)
p1=3;   % ����ʽ�Ľ���
p2=0.8;    % ������Ͷ���ʽ�˺����ı���  p2Ϊ������ı���
sig2=7;     %������˺����Ŀ��  sigma^2
gama=1-0.99;      % ���ӵ�λ����ǰ��ϵ��
step=1;        % �����ƶ��Ĳ���


[AXK,XK,XK1,XK2,XK1T,XK2T,XK11T,XK22T,XK12T,XKKT,B] = ikernel(n,p1,p2,sig2,gama,step);

% ׼���˲�����
L = reshape(XK(((2*n+1)^2+1)/2,:),2*n+1,2*n+1);
Lr = reshape(XK1T(((2*n+1)^2+1)/2,:),2*n+1,2*n+1);
Lc = reshape(XK2T(((2*n+1)^2+1)/2,:),2*n+1,2*n+1);
Lrr=reshape(XK11T(((2*n+1)^2+1)/2,:),2*n+1,2*n+1);
Lrc=reshape(XK12T(((2*n+1)^2+1)/2,:),2*n+1,2*n+1);
Lcc=reshape(XK22T(((2*n+1)^2+1)/2,:),2*n+1,2*n+1);

% �������˲����ӿ��ӻ�

 imstar = data;
Dr=imfilter(imstar,Lr,'conv');
Dc=imfilter(imstar,Lc,'conv');
D1=imfilter(imstar,Lrr,'conv');
D2=imfilter(imstar,Lrr,'conv').*imfilter(imstar,Lcc,'conv')-imfilter(imstar,Lrc,'conv').*imfilter(imstar,Lrc,'conv');

%
ThresholdD1 = 0.01;
ThresholdD2 = 0.01;
ThresholdDr = 1;

result3 = abs(Dr)<ThresholdDr;
result4 = abs(Dc)<ThresholdDr;

result1=D1<ThresholdD1;
result2=D2>ThresholdD2;

% �õ�����ֵ���λ��
result=and(result1,result2);
result = and(result,result3);
result = and(result,result4);

STATS = regionprops(result,'All');
local_max_p = size(STATS,1);
local_max = zeros(local_max_p,3);
cat = 1;
for i = 1 : local_max_p
    pointList = STATS(i).PixelList;
    pointIdxList = STATS(i).PixelIdxList;
    point_num = size(pointList,1);
    point = zeros(point_num,4);
    if STATS(i).Area > 4
        for j = 1 : point_num
            point(j,:) = [pointList(j,:),imstar(pointList(j,2),pointList(j,1)),pointIdxList(j)];
        end
        index_temp = find(point(:,3) == max(point(:,3)));
        local_max(cat,:) = [point(index_temp,[1 2]), pointIdxList(index_temp)];
        cat = cat + 1;
    end
end
local_max(cat:end,:)=[];