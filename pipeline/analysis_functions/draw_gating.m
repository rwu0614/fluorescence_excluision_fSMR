function [gating_x, gating_y] = draw_gating(figure_number)
% Description: draw_gating function serves to choose points on the desired figure in which to use
% for live cell gating, and plots the area of gating on figure based on chosen points.
% Inputs:
%   figure_number = previous figure generated in which to draw gating on
% Outputs:
%   gating_x = array of x-values for each point selected in order of selection
%   gating_y = array of y-values for each point selected in order of selection


gating_x = [];
gating_y = [];

figure(figure_number)
hold on
exitflag = 0;
while exitflag == 0
    [selected_x, selected_y] = ginput(1);
    gating_x(end+1) = selected_x;
    gating_y(end+1) = selected_y;
    plot(gating_x,gating_y,'k','Marker','.')
    choose_points = input('Want to choose more points? (1 = Yes, 0 = No):    ');
    if choose_points == 0
        exitflag = 1;
    else
        exitflag = 0;
    end
end
end
