
function W = get_true_weights(ftrue,tpoints,r,v,fx_cell,fv_cell,lambda)
    n = floor(sqrt(tpoints));
    [rr,th,~,~] = build_polar_grid(n,r,v);
    A = [];
    for i=1:length(fx_cell)
        for j=1:length(fv_cell)
            A = [A reshape(fx_cell{i}(rr).*fv_cell{j}(th),[],1)];
        end
    end
    b = reshape(ftrue(rr,th),[],1);
    W = sparsifyDynamics(A,b,lambda,0,ones(size(A,2),1),size(A,2));
end