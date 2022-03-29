function touch_ = touch_clump(outcat_i, outcat_j,mult)
% 判断两个云核是否重叠，返回1或者0
% mult表示倍数  
% 判断条件：
% 两个云核中心点间的距离d1 
% 两个云核的轴长(sigma)之和的长度对 d2
% d1 > mult* d2 --> 0  没有重叠  否则为重叠

start_ind = 5; 
clump_1 = outcat_i(start_ind:start_ind+5);
clump_2 = outcat_j(start_ind:start_ind+5);

distance_cen = sqrt(sum((clump_1(1:3)-clump_2(1:3)).^2)); % 两个云核中心点间的距离
distance_size = sqrt(sum((clump_1(4:6)+clump_2(4:6)).^2)); % 轴长之和构成的长度

if distance_cen > distance_size * mult
    touch_ = 0;
else
    touch_ = 1;
end
