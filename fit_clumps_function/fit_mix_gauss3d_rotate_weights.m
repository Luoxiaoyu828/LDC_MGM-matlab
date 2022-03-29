function [beta,MSE,R,mdl] = fit_mix_gauss3d_rotate_weights(x_tem,y_tem,beta0,is_rotate)
% 做多核高斯参数拟合的

% inputs:
% x_tem: 坐标
% y_r: 强度值
% beta0: 拟合的初始值

%   单个3维高斯分布：[A,x0,sigma_x,y0,sigma_y,theta,v0,sigma_v]
%   n个3维高斯分布：n*8
%                   [A1,x01,sigma_x1,y01,sigma_y1,theta1,v01,sigma_v1
%                    A2,x02,sigma_x2,y02,sigma_y2,theta2,v02,sigma_v2
%                      ……
%                    An,x0n,sigma_xn,y0n,sigma_yn,thetan,v0n,sigma_vn]

% output:
% beta: 拟合值  和beta0的尺寸一样
% MSE: Mean squared error
% R: Residuals for the fitted model, returned as a vector.

% 参考文献：https://en.wikipedia.org/wiki/Gaussian_function

% 带旋转的3维高斯分布
gauss_3d_rotate = @(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2)))...
    + ((x(:,3)-a(7)).^2)./(2*(a(8))^2))));

% 不带旋转的3维高斯分布
gauss_3d_no_rotate = @(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2)./(2*(a(3))^2) + ((x(:,2)-a(4)).^2)./(2*(a(5))^2)  ...
    + ((x(:,3)-a(6)).^2)./(2*(a(7))^2) ) ) );

if is_rotate
    gauss_3d = gauss_3d_rotate;
else
    gauss_3d = gauss_3d_no_rotate;
    beta0(:,6) = [];
end
% strr = sprintf('@(a, x) ( a(%d,1) * exp( -(((x(:,1)-a(%d,2)).^2).*(cos(a(%d,6))^2/(2*a(%d,3)^2) + sin(a(%d,6))^2/(2*a(%d,5)^2)) + ((x(:,2)-a(%d,4)).^2).* (sin(a(%d,6))^2/(2*a(%d,3)^2) + cos(a(%d,6))^2/(2*a(%d,5)^2)) +  (x(:,1)-a(%d,2)).*(x(:,2)-a(%d,4)).* (2*(-sin(2*a(%d,6))/(4*a(%d,3)^2) + sin(2*a(%d,6))/(4*a(%d,5)^2)))+ ((x(:,3)-a(%d,7)).^2)./(2*(a(%d,8))^2))))',i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);

% strr = '@(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2)) + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2)) +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2)))+ ((x(:,3)-a(7)).^2)./(2*(a(8))^2))))';
% strr1 = '@(a, x) ( a(1,1) * exp( -(((x(:,1)-a(1,2)).^2).*(cos(a(1,6))^2/(2*a(1,3)^2) + sin(a(1,6))^2/(2*a(1,5)^2)) + ((x(:,2)-a(1,4)).^2).* (sin(a(1,6))^2/(2*a(1,3)^2) + cos(a(1,6))^2/(2*a(1,5)^2)) +  (x(:,1)-a(1,2)).*(x(:,2)-a(1,4)).* (2*(-sin(2*a(1,6))/(4*a(1,3)^2) + sin(2*a(1,6))/(4*a(1,5)^2)))+ ((x(:,3)-a(1,7)).^2)./(2*(a(1,8))^2))))';
% strr2 = sprintf('@(a, x) ( a(%d,1) * exp( -(((x(:,1)-a(%d,2)).^2).*(cos(a(%d,6))^2/(2*a(%d,3)^2) + sin(a(%d,6))^2/(2*a(%d,5)^2)) + ((x(:,2)-a(%d,4)).^2).* (sin(a(%d,6))^2/(2*a(%d,3)^2) + cos(a(%d,6))^2/(2*a(%d,5)^2)) +  (x(:,1)-a(%d,2)).*(x(:,2)-a(%d,4)).* (2*(-sin(2*a(%d,6))/(4*a(%d,3)^2) + sin(2*a(%d,6))/(4*a(%d,5)^2)))+ ((x(:,3)-a(%d,7)).^2)./(2*(a(%d,8))^2))))', i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);
% gauss_3d = str2func(strr2);

% 多个成分的维高斯分布
mix_num = size(beta0,1);
% fprintf('fit %d clump(s)\n', mix_num)
% 生成mix_num个高斯函数的表达式
gauss_express = cell(1,mix_num);
for i = 1 : mix_num
   gauss_express{i} = sprintf('( a(%d,1) * exp( -(((x(:,1)-a(%d,2)).^2).*(cos(a(%d,6))^2/(2*a(%d,3)^2) + sin(a(%d,6))^2/(2*a(%d,5)^2)) + ((x(:,2)-a(%d,4)).^2).* (sin(a(%d,6))^2/(2*a(%d,3)^2) + cos(a(%d,6))^2/(2*a(%d,5)^2)) +  (x(:,1)-a(%d,2)).*(x(:,2)-a(%d,4)).* (2*(-sin(2*a(%d,6))/(4*a(%d,3)^2) + sin(2*a(%d,6))/(4*a(%d,5)^2)))+ ((x(:,3)-a(%d,7)).^2)./(2*(a(%d,8))^2))))', i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);
end


% 用str2func将字符串表达式转化成函数句柄 function_handle
gauss_3d_mix = str2func(['@(a,x)(', strjoin(gauss_express, ' + '), ')']);


rng('default') % for reproducibility
opts = statset('nlinfit');

% opts.RobustWgtFun = 'bisquare';
% weights = @(y_tem) y_tem.^2 ./ sum(y_tem .^ 2);

weights = gauss_3d_rotate(beta0,x_tem);
weights = @(weights) weights.^4 ./ sum(weights .^ 4);
% mdl = fitnlm(x_tem, y_tem, gauss_3d_mix,beta0);
mdl = [];
[beta, R, J, CovB, MSE, ErrorModelInfo] = nlinfit(x_tem, y_tem, gauss_3d_mix, beta0, opts,'Weights',weights);
% 将角度由弧度转化为角度
if is_rotate
    beta(:,6) = beta(:,6)./pi.*180;

else
    beta = [beta(:,1:5),zeros(size(beta,1),1),beta(:,6:7)];

end



