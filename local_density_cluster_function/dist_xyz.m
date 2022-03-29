function distance = dist_xyz(point_a, point_b)
temp=point_a-point_b;
% distance=sqrt(temp(1).^2+temp(2).^2+temp(3).^2);
distance = sqrt(sum(temp.^2));