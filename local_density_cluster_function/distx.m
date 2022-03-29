function [distance] = distx(kc1,kc2,xx)
% distance=abs(xx(kc1,1)-xx(kc2,1))+abs(xx(kc1,2)-xx(kc2,2))+abs(xx(kc1,3)-xx(kc2,3));
if size(xx,2)==3
    distance=sqrt((xx(kc1,1)-xx(kc2,1)).^2+(xx(kc1,2)-xx(kc2,2)).^2+(xx(kc1,3)-xx(kc2,3)).^2);
else
    distance=sqrt((xx(kc1,1)-xx(kc2,1)).^2+(xx(kc1,2)-xx(kc2,2)).^2);
end