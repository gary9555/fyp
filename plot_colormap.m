function [h,htxt]=plot_colormap(interval,pos,labels,txt)
% PLOT_COLORMAP
%   interval        - Interval of direction mapping to colors.
%   pos             - Position [NorthEast, NorthWest, SouthEast, SouthWest].
%                     Positions calculated for image mode.
%   labels          - Plotting labels [0 1]?
%   txt             - Text label.
if nargin<1, interval = [-pi pi]; end
if nargin<2, pos = 'NorthEast'; end
if nargin<3, labels = 0; end
if nargin<4, txt = []; end
% Select position x0, y0, and size.
XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');
x0 = XLim(1);
w = XLim(end) - x0;
y0 = YLim(1);
h = YLim(end) - y0;
l = min(h, w)/10;
pos = upper(pos);
% Default position is NorthEast.
s = (1+0.05)*l;
xp = x0 + w - s;
yp = y0 + s;
% if strcmp(pos, 'NORTHEAST'), xp = w - s - 1/8 * s; yp = s;
if strcmp(pos, 'NORTHEAST'), xp = x0 + w - s; yp = y0 + s;
elseif strcmp(pos, 'NORTHWEST'), xp = x0 + s; yp = y0 + s;
elseif strcmp(pos, 'SOUTHEAST'), xp = x0 + w - s; yp = y0 + h - s;
elseif strcmp(pos, 'SOUTHWEST'), xp = x0 + s; yp = y0 + h - s;
else error(['Unknown position ', pos, ' in function plot_colormap!']); 
end
% select
map = hsv;
interval = -interval;
rmin=3/10*s;
rmax=l;
n=length(map);
phi=linspace(interval(1),interval(end),n+1);

% Draw black background.
if abs(interval(end)-interval(1))>pi, h=patch([xp-s xp-s xp+s xp+s],[yp-s yp+s yp+s yp-s],[0 0 0]);
else h=patch([xp-s xp-s xp+s xp+s],[yp-s yp+0.05*l yp+0.05*l yp-s],[0 0 0]); 
end

for i=1:n
    [x,y]=pol2cart([phi(i) phi(i) phi(i+1) phi(i+1)],[rmin rmax rmax rmin]);
    h = [h,patch(xp+x,yp+y,map(n-i+1,:))];
    set(h,'LineStyle','none');
end
% Plot labels in interval.
if labels,
    htxt=[];
    xi=linspace(interval(end)-interval(end)/8,interval(1),8);
    ti=linspace(-interval(1),-(interval(end)-interval(end)/8),8);
    ti=circshift(ti,[0 -1]);
    for i=1:length(xi)
        [xtext,ytext]=pol2cart(xi(i),1.1*rmax);
        htxt=[htxt,text(x0+xtext,y0+ytext,[num2str(ti(i)*180/pi),'°'],'HorizontalAlignment','center')];
    end
end
if ~isempty(txt) 
    htxt=[htxt, text(xp,yp,txt,'HorizontalAlignment','center')];
end


% function [h,htxt]=plot_colormap(interval,label,x0,y0,s)
% % PLOT_COLORMAP
% %   interval        - Interval of direction mapping to colors.
% %   label           - Label in the center of the colormap.
% %   x0              - X-position.
% %   y0              - Y-position.
% %   s               - Size.
% map = hsv;
% interval = -interval;
% rmin=15;
% rmax=50;
% n=length(map);
% phi=linspace(interval(1),interval(end),n+1);
% h=[];
% for i=1:n
%     [x,y]=pol2cart([phi(i) phi(i) phi(i+1) phi(i+1)],[rmin rmax rmax rmin]);
%     h = [h,patch(x,y,map(i,:))];
% end
% 
% % interval = -interval;
% htxt=[];
% xi=linspace(interval(end)-interval(end)/8,interval(1),8);
% ti=linspace(-interval(1),-(interval(end)-interval(end)/8),8);
% ti=circshift(ti,[0 -1]);
% for i=1:length(xi)
%     [xtext,ytext]=pol2cart(xi(i),1.1*rmax);
%     htxt=[htxt,text(xtext,ytext,[num2str(ti(i)*180/pi),'°'],'HorizontalAlignment','center')];
% end
% 
% set(h,'LineStyle','none');
% axis('equal','off',[-1.2 1.2 -1.2 1.2]*rmax);
% 
% htxt=[htxt, text(0,0,label,'HorizontalAlignment','center')];

