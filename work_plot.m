function work_plot(M)
% close all
% clear
% load M.mat
% M=fitsread('dfji.fits');
mMax = max(M(:));
mMin = min(M(:));
% 归一化
% M=(M-min(M(:)))/(max(M(:))-min(M(:)));
M=(M-mMin)/(mMax-mMin);
[size_x,size_y,size_z]=size(M);
% mat用来保存坐标点和强度，如[1,1,1,0.0641675483423438]
mat=zeros(size_x*size_y*size_z,4);
l=1;
for i =1:size(M,1)
    for j =1:size(M,2)
       for k =1:size(M,3)
           mat(l,:)=[i,j,k,M(i,j,k)];
           l=l+1;
       end
    end
end

% dat_sm = smooth3(dat,'box',3);
% 对3D数据进行平滑处理：有box和gaussian两种方式
dat_sm = smooth3(M,'gaussian',5);
% dat_sm =M;

% 计算三个方向的积分图，显示以供检验
M1=sum(M,1);
M1 = reshape(M1,size_y,size_z);

M2=sum(M,2);
M2 = reshape(M2,size_x,size_z);

M3=sum(M,3);
M3 = reshape(M3,size_x,size_y);


[X_s,Y_s]=meshgrid([1:size_y],[1:size_x]);

n_layer = 10;
col_mm=flipud(jet(n_layer)); 


lev = 1:1:n_layer+1;
lev = lev./11;
hold on
for i = 1: n_layer
    mm_lev=lev(i);
%     patch(X,Y,C)   C:表示颜色
    p1 = patch(isosurface(dat_sm,mm_lev),'FaceColor',col_mm(i,:),'EdgeColor','None','FaceAlpha',0.1+0.015*i);
    isonormals(dat_sm,p1);
    axis([0 size_x 0 size_y 0 size_z])
end
% view(30,30)
view(43.6 ,35.6);
% 设置颜色范围
% caxis([0 1])
% colormap(flipud(jet))
% colorbar
camlight headlight; 
lighting phong
set(gca,'GridLineStyle',':')
grid on
axis([0 size_y 0 size_x 0 size_z])
% axis equal

lighting gouraud
