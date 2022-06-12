%% plot interaction force

n = 200;
r=0.5;
v=[0 1];
cap = inf;
levs = 200;
[rr,th,xx,yy] = build_polar_grid(n,r,v);

try 
    F = f_learned(rr,th);
catch
    F = [];
end

plot_forces(F,cap,levs,1,xx,yy)

%% plot alignment force

n = 200;
r=0.5;
v=[0 1];
cap = inf;
levs = 200;
[rr,th,xx,yy] = build_polar_grid(n,r,v);

try 
    H = h_learned(rr,th);
catch
    H = [];
end
plot_forces(H,cap,levs,1,xx,yy)

%% plot drag force

n = 200;
r=0.05;
v=[0 1];
cap = inf;
levs = 200;
[rr,th,xx,yy] = build_polar_grid(n,r,v);

try 
    D = d_learned(rr,th);
catch
    D = [];
end

plot_forces(D,cap,levs,1,xx,yy)

function plot_forces(F,cap,levs,fignum,xx,yy)

    if ~isempty(F)
        x = xx(1,:);
        y = yy(:,1);
        F = min(max(F,-cap),cap);
        figure(fignum);
        contourf(xx,yy,F,levs,'edgecolor','none')
        xlabel('x')
        ylabel('y')
        ylim([min(y) max(y)])
        xlim([min(x) max(x)])
        colorbar
        colormap(jet(100))
        axis equal
    end

end

function [rr,th,xx,yy] = build_polar_grid(n,r,v)

    x = linspace(-r,r,n);
    y = linspace(-r,r,n);
    [xx,yy] = meshgrid(x,y);
    dX = [xx(:) yy(:)];
    dX = reshape([0 0],1,1,2) - reshape(dX,1,n^2,2);
    rr = reshape(squeeze(vecnorm(dX,2,3)),n,[]);

    dV = repmat(reshape(v,1,1,2),1,n^2,1);
    th = dot(-dX,dV,3);
    th = reshape(atan2(abs(dot(-flip(cat(3,-dX(:,:,1),dX(:,:,2)),3),dV,3)),th),n,[]);

end