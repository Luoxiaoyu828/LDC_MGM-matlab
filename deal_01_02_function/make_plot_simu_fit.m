function make_plot_simu_fit(data_x,data_y,x_range,y_range,x_label, y_label)
% ps = s_position;
%     pf = f_position;
    figure
    for i = 1:3
        subplot(1,3,i)
        plot(data_x(:,i),data_y(:,i),'*')
        hold on
        plot(x_range(i,:),y_range(i,:),'b-')
        ylabel(y_label)
        xlabel(x_label)
    end
%     subplot(1,3,2)
%     plot(ps(:,2),pf(:,2),'r*')
%     hold on
%     plot(0:0.01:x_range,0:0.01:x_range,'b.')
%     ylabel('fit\_cen')
%     xlabel('simulate\_cen')
%     subplot(1,3,2)
%     plot(ps(:,2),pf(:,2),'r*')
%     hold on
%     plot(0:0.01:x_range,0:0.01:x_range,'b.')
%     ylabel('fit\_cen')
%     xlabel('simulate\_cen')
%     subplot(1,3,3)
%     plot(ps(:,3),pf(:,3),'r*')
%     hold on
%     plot(0:0.01:x_range,0:0.01:x_range,'b.')
%     ylabel('fit\_cen')
%     xlabel('simulate\_cen')