% Make Plot Axes Nice
% Script to improve plots taken from STANDARDIZE_FIGURE within PlottingTemplate_Short.m

% INCLUDE PLOT_STANDARDS
    
    PS = PLOT_STANDARDS();

%get axes for current figure
ax = fig.CurrentAxes;        


% ADJUST AXES PROPERTIES
        
        ax.Box = 'on';
        ax.TickDir = 'out';
        ax.TickLength = [PS.AxisTickLength, PS.AxisTickLength];
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        %ax.ZMinorTick = 'on';
        ax.XColor = PS.AxisColor;
        ax.YColor = PS.AxisColor;
        %ax.ZColor = PS.AxisColor;
        % ax.XLabel.Color = PS.AxisLabelColor;
        % ax.YLabel.Color = PS.AxisLabelColor;
        % ax.ZLabel.Color = PS.AxisLabelColor;
        ax.LineWidth = PS.DefaultLineWidth;

        ax.FontWeight = 'bold';