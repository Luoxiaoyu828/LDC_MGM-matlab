function touch_ = touch_clump(outcat_i, outcat_j,mult)
% �ж������ƺ��Ƿ��ص�������1����0
% mult��ʾ����  
% �ж�������
% �����ƺ����ĵ��ľ���d1 
% �����ƺ˵��᳤(sigma)֮�͵ĳ��ȶ� d2
% d1 > mult* d2 --> 0  û���ص�  ����Ϊ�ص�

start_ind = 5; 
clump_1 = outcat_i(start_ind:start_ind+5);
clump_2 = outcat_j(start_ind:start_ind+5);

distance_cen = sqrt(sum((clump_1(1:3)-clump_2(1:3)).^2)); % �����ƺ����ĵ��ľ���
distance_size = sqrt(sum((clump_1(4:6)+clump_2(4:6)).^2)); % �᳤֮�͹��ɵĳ���

if distance_cen > distance_size * mult
    touch_ = 0;
else
    touch_ = 1;
end
