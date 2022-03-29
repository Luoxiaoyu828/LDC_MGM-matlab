function index_value = kc_coord_3d(delta_ii_xy,xm,ym,zm,r)
jt = delta_ii_xy(1);
it = delta_ii_xy(2);
kt = delta_ii_xy(3);
% r=2;
% bt=a(max(1,kt-r):min(zm,kt+r), max(1,jt-r):min(ym,jt+r), max(1,it-r):min(xm,it+r));
% bt1=a(max(1,kt-r):min(zm,kt+r), max(1,it-r):min(xm,it+r), max(1,jt-r):min(ym,jt+r));
% bt2=a(max(1,it-r):min(xm,it+r), max(1,kt-r):min(zm,kt+r), max(1,jt-r):min(ym,jt+r) );
% bt3=a(max(1,it-r):min(xm,it+r), max(1,jt-r):min(ym,jt+r), max(1,kt-r):min(zm,kt+r) );
% bt4=a(max(1,jt-r):min(ym,jt+r), max(1,kt-r):min(zm,kt+r),max(1,it-r):min(xm,it+r) );
% bt5=a(max(1,jt-r):min(ym,jt+r),max(1,kt-r):min(zm,kt+r),max(1,it-r):min(xm,it+r) );

[p_i,p_j,p_k] = meshgrid(max(1,kt-r):min(xm,kt+r), max(1,it-r):min(ym,it+r), max(1,jt-r):min(zm,jt+r));
index_value=[p_k(:),p_j(:),p_i(:)];