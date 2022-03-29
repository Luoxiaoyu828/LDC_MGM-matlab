function [angle_s_,angle_f_] = get_angle(Match_Outcat)
% 拟合结果中的角度和仿真数据中的角度对应关系
aa_result = Match_Outcat(:,[18 38 16 17 35 37 12 33 1 2]);
for aa_result_ii=1:size(aa_result,1)
    size_x = aa_result(aa_result_ii,5);
    size_y = aa_result(aa_result_ii,6);
    s_angle = aa_result(aa_result_ii,1);
    fit_angle = aa_result(aa_result_ii,2);
    if size_x>size_y
        s_angle_ = mod(s_angle,180);
        fit_angle_ = mod(fit_angle,180);
    else
        s_angle_ = mod(s_angle,180);
        fit_angle_ = mod(fit_angle+90,180);
    end
    aa_result(aa_result_ii,[9 10])= [s_angle_,fit_angle_];
end
angle_s_ = aa_result(:,9);
angle_f_ = aa_result(:,10);