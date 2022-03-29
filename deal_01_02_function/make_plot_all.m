function make_plot_all(s_position, s_size, s_sum, s_Peak, angle_s_,f_position, f_size, f_sum, f_Peak, angle_f_)
% Î»ÖÃÆ«²î*********  1
x_range = [[0 120];[0 120];[0 120]];
y_range = [[0 120];[0 120];[0 120]];
data_x = s_position;
data_y = f_position;
make_plot_simu_fit(data_x,data_y,x_range,y_range,'simulated\_cen', 'fit\_cen')

x_range = [[3.5 4];[2 2.5];[3 5]];
y_range = [[3.5 4];[2 2.5];[3 5]];
make_plot_simu_fit(s_size,f_size,x_range,y_range,'simulated\_size', 'fit\_size')

% SumµÄÆ«²î
x_range = [1000, 6000];
y_range = [1000, 6000];
x_label = 'simulate\_Sum';
y_label = 'fit\_Sum';
data_x = s_sum;
data_y = f_sum;
make_plot_sum_peak_sita(data_x, data_y,x_range,y_range,x_label, y_label)
make_plot_sum_peak_sita(angle_s_, angle_f_,[0 180],[0 180], 'simulate\_sita','fit\_sita')
make_plot_sum_peak_sita(s_Peak, f_Peak,[2 10],[2 10], 'simulate\_peak','fit\_peak')