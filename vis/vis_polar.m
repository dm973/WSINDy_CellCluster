ftag=0;
figstart=1;
n = 200;rr=0.5;
x = linspace(-rr,rr,n);
y = linspace(-rr,rr,n);
[xx,yy] = meshgrid(x,y);
dX = [xx(:) yy(:)];
dX = reshape([0 0],1,1,2) - reshape(dX,1,n^2,2);
rr = reshape(squeeze(vecnorm(dX,2,3)),n,[]);
v = [0 1];
dV = repmat(reshape(v,1,1,2),1,n^2,1);
th = dot(-dX,dV,3);
th = reshape(atan2(abs(dot(-flip(cat(3,-dX(:,:,1),dX(:,:,2)),3),dV,3)),th),n,[]);
% figure(1);clf;
% figure(2);clf;
% figure(3);clf;

if isequal(ftag,'true')
    try 
        F = f_xv_true(rr,th);
    catch
        F = [];
    end
else
    try 
        F = f_learned(rr,th);
    catch
        F = [];
    end
end

if ~isempty(F)
    cap = inf;
    F = min(max(F,-cap),cap);
    levs = 200;
    figure(figstart);
    contourf(xx,yy,F,levs,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    ylim([min(y) max(y)])
    xlim([min(x) max(x)])
    colorbar
    colormap(jet(100))
    axis equal
%     saveas(gcf,['~/Desktop/',num2str(inds_cell{1}(1)),'_ftrue.png'])
else
    close(figure(figstart));
end


if isequal(ftag,'true')
    try 
        F = h_xv_true(rr,th);
    catch
        F = [];
    end
else
    try 
        F = h_learned(rr,th);
    catch
        F = [];
    end
end

if ~isempty(F)
    cap = inf;
    F = min(max(F,-cap),cap);
    levs = 200;
    figure(figstart+1);
    contourf(xx,yy,F,levs,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    ylim([min(y) max(y)])
    xlim([min(x) max(x)])
    colorbar
    colormap(jet(100))
    axis equal
%     saveas(gcf,['~/Desktop/',num2str(inds_cell{1}(1)),'_htrue.png'])
else
    close(figure(figstart+1));
end

n = 200;rr=0.05;
x = linspace(-rr,rr,n);
y = linspace(-rr,rr,n);
[xx,yy] = meshgrid(x,y);
dX = [xx(:) yy(:)];
dX = reshape([0 0],1,1,2) - reshape(dX,1,n^2,2);
rr = reshape(squeeze(vecnorm(dX,2,3)),n,[]);
v = [0 1];
dV = repmat(reshape(v,1,1,2),1,n^2,1);
th = dot(-dX,dV,3);
th = reshape(atan2(abs(dot(-flip(cat(3,-dX(:,:,1),dX(:,:,2)),3),dV,3)),th),n,[]);

if isequal(ftag,'true')
    try 
        F = d_xv_true(rr,th);
    catch
        F = [];
    end
else
    try 
        F = d_learned(rr,th);
    catch
        F = [];
    end
end

if ~isempty(F)
    cap = inf;
    F = min(max(F,-cap),cap);
    levs = 200;
    figure(figstart+2);
    contourf(xx,yy,F,levs,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    ylim([min(y) max(y)])
    xlim([min(x) max(x)])
    colorbar
    colormap(jet(100))
    axis equal
%     saveas(gcf,['~/Desktop/',num2str(inds_cell{1}(1)),'_dtrue.png'])
else
    close(figure(figstart+2));
end