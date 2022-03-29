function err_indx = get_fit_error_index(Match_Outcat)
% 
cen1 = Match_Outcat(:,[34 36 39]);
[index_x,~] = find(cen1<0 | cen1>120);
index1 =unique(index_x);

size = Match_Outcat(:,[35 37 40]);
[index_x,~] = find(size<1 | size>15);
index2 =unique(index_x);

err_indx = unique([index1',index2']);