

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