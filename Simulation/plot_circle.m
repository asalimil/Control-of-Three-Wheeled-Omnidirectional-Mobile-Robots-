function plot_circle(O,R,boundarycolor,fill_obs)

x = O(1); y = O(2);
% t = 0:pi/50:2*pi;
% circle = [];

if fill_obs == 0
    %     for t = 1:length(t)
    %         circle = [circle; O(1)+R*cos(t) O(2)+R*sin(t)];
    %     end
    %     plot(circle(:,1), circle(:,2),'.','Color',rgb./255,'LineWidth',3); hold on;
    th = 0:pi/50:2*pi;
    x_circle = R*cos(th) + x;
    y_circle = R*sin(th) + y;
    plot(x_circle, y_circle,'Color',boundarycolor/255,'LineWidth',3); hold on;
    % fill(x_circle, y_circle, c)
else
    %     for t = 1:length(t)
    %         circle = [circle; O(1)+R*cos(t) O(2)+R*sin(t)];
    %     end
    %     plot(circle(:,1), circle(:,2),'.','Color',rgb./255,'LineWidth',3); hold on;
    %     fill(circle(:,1), circle(:,2),'r'); hold on;
    th = 0:pi/50:2*pi;
    x_circle = R*cos(th)+x;
    y_circle = R*sin(th)+y;
    plot(x_circle, y_circle,'Color',boundarycolor/255,'LineWidth',3); hold on;
    fill(x_circle, y_circle, 'y')
end

drawnow;

end