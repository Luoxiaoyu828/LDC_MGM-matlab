function distance = dist_xyz_2d(point_a, point_b)
temp=point_a-point_b;
distance=sqrt(temp(1).^2+temp(2).^2);