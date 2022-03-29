function [beta,MSE,R,Sum] = fit_gauss3d_rotate(x_tem,y_r,beta0)
% inputs:
% x_tem: ����
% y_r: ǿ��ֵ
% beta0: ��ϵĳ�ʼֵ
%   ����2ά��˹�ֲ���[A,x0,sigma_x,y0,sigma_y,theta]
%   ����3ά��˹�ֲ���[A,x0,sigma_x,y0,sigma_y,theta,v0,sigma_v]
%   ������2ά��˹�ֲ���[A,x0,sigma_x,y0,sigma_y,theta,A1,x01,sigma_x1,y01,sigma_y1,theta1]
%   ����3ά��˹�ֲ���[A,x0,sigma_x,y0,sigma_y,theta,v0,sigma_v
%                    A1,x01,sigma_x1,y01,sigma_y1,theta1,v01,sigma_v1]

% output:
% beta: ���ֵ
% MSE: Mean squared error
% R: Residuals for the fitted model, returned as a vector.

% �ο����ף�https://en.wikipedia.org/wiki/Gaussian_function

% a=[10,30,7,30,9,45/180*pi];
% A=10
% x0=30
% y0=30
% sigma_X = 7
% sigma_Y = 9
% Sita = 45��
Sum = 0;
% ����ת��2ά��˹�ֲ�
gauss_2d = @(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2)))...
    )));
gauss_2d_b = @(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2)))...
    )) + a(7));

% ����ת��2ά��˹�ֲ���������˹�ɷ֣�
gauss_2d_plus = @(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2))) ))...
    + a(7) * exp( -(((x(:,1)-a(8)).^2).*(cos(a(12))^2/(2*a(9)^2) + sin(a(12))^2/(2*a(11)^2))...
    + ((x(:,2)-a(10)).^2).* (sin(a(12))^2/(2*a(9)^2) + cos(a(12))^2/(2*a(11)^2))...
    +  (x(:,1)-a(8)).*(x(:,2)-a(10)).* (2*(-sin(2*a(12))/(4*a(9)^2) + sin(2*a(12))/(4*a(11)^2)))...
    )));

% ������ת��2ά��˹�ֲ�
% gauss_2d_0 = @(a,x) ( a(1) * exp( -(((x(:,1)-a(2)).^2)./(2*(a(3))^2) + ((x(:,2)-a(4)).^2)./(2*(a(5))^2) )));
%     a = cos(theta)^2/(2*sigma_X^2) + sin(theta)^2/(2*sigma_Y^2);
%     b = -sin(2*theta)/(4*sigma_X^2) + sin(2*theta)/(4*sigma_Y^2);
%     c = sin(theta)^2/(2*sigma_X^2) + cos(theta)^2/(2*sigma_Y^2);
%     Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2));
% ����ת��3ά��˹�ֲ�
gauss_3d = @(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2)))...
    + ((x(:,3)-a(7)).^2)./(2*(a(8))^2))));


% ����ת������3ά��˹�ֲ�
gauss_3d_plus = @(a, x) (...
    a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2)))...
    + ((x(:,3)-a(7)).^2)./(2*(a(8))^2)))...
    + a(9) * exp( -(((x(:,1)-a(10)).^2).*(cos(a(14))^2/(2*a(11)^2) + sin(a(14))^2/(2*a(13)^2))...
    + ((x(:,2)-a(12)).^2).* (sin(a(14))^2/(2*a(11)^2) + cos(a(14))^2/(2*a(13)^2))...
    +  (x(:,1)-a(10)).*(x(:,2)-a(12)).* (2*(-sin(2*a(14))/(4*a(11)^2) + sin(2*a(14))/(4*a(13)^2)))...
    + ((x(:,3)-a(15)).^2)./(2*(a(16))^2)))...
    );
gauss_2d_plus_b =@(a, x) ( a(1) * exp( -(((x(:,1)-a(2)).^2).*(cos(a(6))^2/(2*a(3)^2) + sin(a(6))^2/(2*a(5)^2))...
    + ((x(:,2)-a(4)).^2).* (sin(a(6))^2/(2*a(3)^2) + cos(a(6))^2/(2*a(5)^2))...
    +  (x(:,1)-a(2)).*(x(:,2)-a(4)).* (2*(-sin(2*a(6))/(4*a(3)^2) + sin(2*a(6))/(4*a(5)^2))) ))...
    + a(7) * exp( -(((x(:,1)-a(8)).^2).*(cos(a(12))^2/(2*a(9)^2) + sin(a(12))^2/(2*a(11)^2))...
    + ((x(:,2)-a(10)).^2).* (sin(a(12))^2/(2*a(9)^2) + cos(a(12))^2/(2*a(11)^2))...
    +  (x(:,1)-a(8)).*(x(:,2)-a(10)).* (2*(-sin(2*a(12))/(4*a(9)^2) + sin(2*a(12))/(4*a(11)^2)))...
    )) +a(13));

rng('default') % for reproducibility
opts = statset('nlinfit');

opts.RobustWgtFun = 'cauchy';

if length(beta0)==6
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x_tem,y_r,gauss_2d,beta0,opts);
    beta(6) = beta(6)/pi*180;
end
if length(beta0)==7
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x_tem,y_r,gauss_2d_b,beta0,opts);
    beta(6) = beta(6)/pi*180;
end
if length(beta0)==12
    [beta,R,J,CovB,MSE,ErrorModelInfo]  = nlinfit(x_tem,y_r,gauss_2d_plus,beta0,opts);
    beta([6,12]) = mod(beta([6,12])./pi.*180,360);
end
if length(beta0)==13
    [beta,R,J,CovB,MSE,ErrorModelInfo]  = nlinfit(x_tem,y_r,gauss_2d_plus_b,beta0,opts);
    beta([6,12]) = mod(beta([6,12])./pi.*180,360);
end
if length(beta0)==8
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x_tem,y_r,gauss_3d,beta0,opts);
%     beta(6) = mod(beta(6)/pi*180,360);
    beta(6) = beta(6)/pi*180;
end
if length(beta0)==16
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x_tem,y_r,gauss_3d_plus,beta0,opts);
    beta([6,14]) = mod(beta([6,14])./pi.*180,360);
end
% if length(beta0)==8 && weigths
%     yr_weigths = y_r.^2./sum(y_r);
%     [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x_tem,y_r,gauss_3d,beta0,opts,'Weights',yr_weigths);
%     beta(6) = mod(beta(6)/pi*180,360);
% end
end