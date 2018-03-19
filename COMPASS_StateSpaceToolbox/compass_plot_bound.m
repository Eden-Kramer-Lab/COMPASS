function h=compass_plot_bound(mode, x ,y, y_low, y_high, x_title,y_title,fig_title,color)
% plot figures    
    if nargin < 5
        disp('Not enough input, quit');
        return;
    end
    if nargin ==5
        x_title = [];
        y_title = [];
        fig_title = [];
        color =[ 0.9 0.9 0.9];
    end
    if nargin ==6
        y_title = [];
        fig_title = [];
        color =[ 0.9 0.9 0.9];
    end
    if nargin ==7
        fig_title = [];
        color =[ 0.9 0.9 0.9];
    end
    

   % plots the data plus its bound
    x_inc = x;
    x_dec = x(end:-1:1);
    set(fill([x_inc x_dec],[y_high y_low(end:-1:1)],color),'EdgeColor',color);
    hold on;
    h=plot(x,y,'k');
    set(h,'LineWidth',3);
    if mode == 1    % single figure
        hold off;
    end
    if isempty(x_title)==0
        xlabel(x_title);
    end
    if isempty(y_title)==0
        ylabel(x_title);
    end
    if isempty(fig_title)==0
        title(fig_title);
    end

end

