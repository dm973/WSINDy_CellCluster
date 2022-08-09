surf(xx,yy,F_dat,'edgeColor','none')
colorbar
view([0 90])
if force_ind<3
    xlabel('$x$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    axis equal
else
    xlabel('$s$','interpreter','latex')
    ylabel('$\theta$','interpreter','latex')
end
ylim([min(y) max(y)])
xlim([min(x) max(x)])
colorbar('ticklabelinterpreter','latex','fontsize',14)
colormap(jet(100))
set(gca,'ticklabelinterpreter','latex','fontsize',14)
caxis(c_range)