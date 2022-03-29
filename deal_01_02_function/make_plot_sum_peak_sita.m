function make_plot_sum_peak_sita(data_x, data_y,x_range,y_range,x_label, y_label)
figure
plot(data_x, data_y,'r*')
hold on
plot(x_range,y_range,'b-')
ylabel(x_label)
xlabel(y_label)