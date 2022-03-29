function [result, mask_id, connect_outcat] = connect_clump(loc_outcat,mult)
% 判断仿真云核之间是否相连
% loc_outcat为局部聚类算法中的核表
% mult 为允许的倍数


outcat = loc_outcat;
% outcat = [id, Peak1,Peak2,Peak3,Cen1,Cen2,Cen3,Size1,Size2,Size3,Sum,Peak,Volume]
col_outcat = size(outcat,1);
record = cell(1,col_outcat);

for i = 1:col_outcat
    outcat_i = outcat(i,:);
    a = [];
    j_index = i: col_outcat;
    for j = j_index
        outcat_j = outcat(j,:);
        if touch_clump(outcat_i, outcat_j,mult)
            a = [a,j];
        end
    end
    record{i} = a; % 与第i号云核连接的云核编号
end

connect_outcat = record;
record_ = cell(1,col_outcat);
for i = 1:col_outcat
    cluster_i = record{i};
    for j = 1: col_outcat
        cluster_j = record{j};
       if  intersect(cluster_i,cluster_j)
           cluster_i = union(cluster_i,cluster_j);
       end
    end
    record_{i} = cluster_i;
end

% 循环一下，让相连的充分连接
for ii = 1:10
    record = record_;
    record_ = cell(1,col_outcat);
    for i = 1:col_outcat
        cluster_i = record{i};
        for j = 1: col_outcat
            cluster_j = record{j};
            if  intersect(cluster_i,cluster_j)
                cluster_i = union(cluster_i,cluster_j);
            end
        end
        record_{i} = cluster_i;
    end
end

record_1 = record_;
for i = 1:col_outcat
    cluster_i = record_{i};
    for j = i+1:col_outcat
        cluster_j = record_{j};
        if length(intersect(cluster_i,cluster_j)) > 0
            if all(intersect(cluster_i,cluster_j)==cluster_j) && all(union(cluster_i,cluster_j)==cluster_i)
                record_1{j}=[];
            end
        end
    end
end
result = {};
j=1;
for i = 1:col_outcat
    cluster_i = record_1{i};
    if length(cluster_i)>0
        result{j} = cluster_i;
        j=j+1;
    end
end



mask_id = {};
cluster_id = loc_outcat(:,1);
jj=1;
for i = 1:length(result)
    temp = [];
    for j = result{i}
       temp = [temp,cluster_id(j)];
    end
    mask_id{jj}=temp;
    jj = jj + 1;
end

   



