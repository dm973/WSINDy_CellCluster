%% plot interaction force

n = 200;
r=0.2;
v=[0 1];
cap = inf;
cmap = jet(100);
[rr,th,xx,yy] = build_polar_grid(n,r,v);
try 
    F = f_learned(rr,th);
catch
    F = [];
end
plot_forces(F,cap,1,xx,yy,cmap)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
legend({'$f^\star_{a-r}$'},'interpreter','latex','fontsize',12)
caxis([-15 15])

%% plot alignment force

n = 200;
r=0.2;
v=[0 1];
cap = inf;
cmap = jet(100);
[rr,th,xx,yy] = build_polar_grid(n,r,v);
try 
    H = h_learned(rr,th);
catch
    H = [];
end
plot_forces(H,cap,2,xx,yy,cmap)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)
caxis([-15 15])

%% plot drag force

n = 200;
r=0.05;
v=[0 1];
cap = inf;
cmap = jet(100);
[rr,th,xx,yy] = build_polar_grid(n,r,v);
try 
    D = d_learned(rr,th);
catch
    D = [];
end
plot_forces(D,cap,3,xx,yy,cmap)
xlabel('$v_x$','interpreter','latex')
ylabel('$v_y$','interpreter','latex')
legend({'$f^\star_{drag}$'},'interpreter','latex','fontsize',12)
caxis([-1 1])