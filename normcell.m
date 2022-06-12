function [hv,ev] = normcell(Xscell,numbins)
    n = length(Xscell);
    dims = zeros(n,2);
    for nn=1:n
        dims(nn,:) = size(Xscell{nn},[1 3]);
    end
    Xnorms = zeros(sum(prod(dims,2)),1);
    ntot = 0;
    for nn=1:n
        Xnorms_temp = reshape(vecnorm(Xscell{nn},2,2),[],1);
        L = length(Xnorms_temp);
        Xnorms(ntot+1:ntot+L) = Xnorms_temp;
        ntot = ntot+L;
    end
    [hv,ev] = histcounts(Xnorms,numbins,'normalization','pdf');
end