function [hc] = create_2D_plot(xdata,ydata,n,colormapName,clims)
%CREATE_2D_PLOT rotate the solution to create a 2D plot
% n: #divisions/grid points
% colormapName: name of the colormap; if it is just a color name, we will
% interpolate from black to this color
% clims: [climMin,climMax]
% hc: handle for the colorbar

% select colormap
if nargin >= 4 && exist('colormapName','var')
    % built-in colormaps
    colormapList = {'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', ...
                'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', ...
                'lines', 'colorcube', 'prism', 'flag'};
    if ~any(contains(colormapList,colormapName))
        try
            RGBcolor = validatecolor(colormapName);
            colormapName = (linspace(0,1,256).')*RGBcolor;
        catch
            warning('Wrong colormapName')
            return;
        end
    end
end

% start time
t0 = min(xdata);
% simulation length
tsim = max(xdata);
lim = tsim-t0;%/sqrt(2);
% rotating
x = linspace(-lim,lim,n);
[X,Y] = meshgrid(x);
R = (X.^2+Y.^2).^0.5;

% interpolate for the simulation space points
C = zeros(size(X));
C(:) = interp1(xdata,ydata,t0+R(:));

surf(X,Y,C)
shading interp
view(2)
grid off
axis equal
axis off
axis([-lim lim -lim lim])
colormap(gca,colormapName)
hc = colorbar;
% hc.Ruler.TickLabelFormat = '%2.2f';
% ylabel(hc,'\bf\boldmath$\hat{G}$ / -','FontSize',18,'interpreter','latex')
if nargin >= 5
    clim(clims)
end
set(gca,'FontSize',16)
set(hc,'FontSize',16)
set(gca,'LineWidth',1.5)
set(hc,'LineWidth',1.5)
% Show tick marks above plot
set(gca,'Layer','top')
set(gca,'FontSize',16)
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% xlabel('\bf\boldmath$\hat{x}$ / -','FontSize',18,'interpreter','latex')
% ylabel('\bf\boldmath$\hat{y}$ / -','FontSize',18,'interpreter','latex')
box on


end