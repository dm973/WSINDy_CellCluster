function plot_forces(F,cap,fignum,xx,yy,cmap)

    if ~isempty(F)
        x = xx(1,:);
        y = yy(:,1);
        F = min(max(F,-cap),cap);
        figure(fignum);
        surf(xx,yy,F,'edgecolor','none')
        view([0 90])
        axis equal
        ylim([min(y) max(y)])
        xlim([min(x) max(x)])
        colorbar('ticklabelinterpreter','latex','fontsize',14)
        colormap(cmap)
        set(gca,'ticklabelinterpreter','latex','fontsize',14)
    end

end