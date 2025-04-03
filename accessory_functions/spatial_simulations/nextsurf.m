function hs = nextsurf(c,Lx1,Lx2,nx1,nx2,name)
hs = surf(linspace(0,Lx2,nx2+1),linspace(0,Lx1,nx1+1),c);
% shading interp
shading flat
axis equal
hc = colorbar;
% title(['t = ' num2str(t,'%2.1f')])
axis([0 Lx2 0, Lx1])
% clim([0 1])
hc.Ruler.TickLabelFormat = '%.2f';
ylabel(hc,name,'FontSize',18,'interpreter','none')%,'FontWeight','bold')
xtickformat('%.1f')
ytickformat('%.1f')
set(gca,'FontSize',16)
set(hc,'FontSize',16)
set(gca,'LineWidth',1.5)
set(hc,'LineWidth',1.5)
xlabel('\bf\boldmath$x$','FontSize',18,'interpreter','latex')
ylabel('\bf\boldmath$y$','FontSize',18,'interpreter','latex')
view(2)
% Show tick marks above plot
set(gca,'Layer','top')
grid off
box on
% set(gca,'xaxisLocation','top')
end