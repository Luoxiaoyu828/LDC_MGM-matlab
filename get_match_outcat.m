function [match_num, match_record,F1,find_outcat_points,simulate_outcat_points]=get_match_outcat(simulate_outcat_path,find_outcat_points, dc)
% 计算仿真核与聚类中心的匹配情况
% 输入：simulate_points_filepath,cluster_points_filepath, dc
% 仿真核路径，检测核表路径，匹配上核的最大距离
% 输出：match_num, match_record
% 匹配上的点数，每个点的对应情况（检测编号，仿真编号，两者间的距离）

% simulate_outcat_path='check_data\simulate_20_40_60_outcat.FIT';
% find_outcat_path='check_data\gauss_outcat.FIT';
if ~isstr(find_outcat_points)
    find_outcat_points=find_outcat_points(:,[4 5 6]);
else
    find_outcat_points = importdata(find_outcat_points);
    find_outcat_points = find_outcat_points.data(:,[6,5,7]);
end
temp = split(simulate_outcat_path,'.');
temp = temp{2};
if temp(1)=='F'
    simulate_points = fitsread(simulate_outcat_path,'binarytable');
    simulate_outcat_points = [simulate_points{6},simulate_points{5},simulate_points{7}];
    % disp(simulate_points.textdata)
    
else
    simulate_points = importdata(simulate_outcat_path);
    simulate_outcat_points = simulate_points.data(:,[6,5,7]);
end


% % 计算距离矩阵，方式1和方式2  点数比较少时方式1快  反之方式2快
% 时间已过 0.004156 秒。
% 时间已过 0.037077 秒。

% 计算距离矩阵，方式1
distance = pdist2(find_outcat_points, simulate_outcat_points);
% 计算距离矩阵，方式2
% size_i=size(find_outcat_points,1);
% size_j=size(simulate_outcat_points,1);
% distance=zeros(size_i,size_j);
% for i=1:size_i
%     for j=1:size_j
%         distance(i,j)=dist_xyz(find_outcat_points(i,:), simulate_outcat_points(j,:));
%     end
% end


match_num = 0;
match_num_ = 0;
max_d = max(distance(:)) * 1.2;
match_record = [];
while 1
    distance_min = min(distance(:));
    [index_simulated, index_find] = find(distance == distance_min);
    simulated_points = simulate_outcat_points(index_find,:);
    find_points = find_outcat_points(index_simulated,:);
    match_num_ = match_num;
    %     dc = [2 2 2];
    % size = find_outcat_points(index_simulated,:);
    % dc = min(dc, size);  注意一下  轴的对应关系
    
    
    if sum(abs(simulated_points - find_points) <= dc) == 3
        match_record = [match_record; [index_simulated, index_find, distance_min]];
        distance(index_simulated, :) = ones(1, size(distance,2)) .* max_d;
        distance(:, index_find) = ones(size(distance,1),1) .* max_d;
        match_num = match_num + 1;  
    end
    if match_num_ == match_num
       break 
    end
%     if min(distance(:)) > dc * sqrt(3) || match_num >= size(simulate_outcat_points,1)
%         break
%     end
end
size_i = size(simulate_outcat_points,1);
size_j = size(find_outcat_points,1);

F1.precision = match_num / size_j;
F1.recall = match_num / size_i;
F1.f1 = 2 * F1.precision * F1.recall / (F1.precision + F1.recall);

% simulate_50_model=fitsread('simulate_13CO\simulate_13CO_50_model.fits');
% simulate_points=fitsread('simulate_13CO\simulate_13CO_50_outcat.fits', 'binarytable');
% simulate_points=[simulate_points{5},simulate_points{6},simulate_points{7}];
% simulate_points1(:,1)=-(simulate_points(:,1)-18.2499)*60*2+1;
% simulate_points1(:,2)=(simulate_points(:,2)-(-0.0003))*60*2+1;
% simulate_points1(:,3)=(simulate_points(:,3)-15939.8)/167+1;
% figure
% work_plot(simulate_50_model,max(simulate_50_model(:)),min(simulate_50_model(:)))
% hold on
% plot3(simulate_points1(:,1),simulate_points1(:,2),simulate_points1(:,3),'*')




