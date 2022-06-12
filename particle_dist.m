function [h,e]  = particle_dist(Xs,numbins,maxents,R)
    if ~exist('R','var')
        R=inf;
    end
    n = length(Xs);
    Xrad_cell = cell(n,1);
    Ls = zeros(n+1,1);
    for nn=1:n
        X = Xs{nn}(1:min(end,1000),:,:);
        normX =vecnorm(X,2,2);
        X = X./normX.*(min(R,normX));
        [N,~,M] = size(X);
        kk = ceil(N*(N-1)/2*M*n/maxents);
        Mk = floor((M-1)/kk)+1;
        Xrad = zeros(N*(N-1)/2,Mk);
        Ls(nn+1) = N*(N-1)/2*Mk;
        [N,~,~] = size(X);    
        ind_mat = logical(~triu(ones(N,N))');
        for m=1:kk:M
            radmat = squeeze(dist(X(:,:,m)'));
            radmat = triu(radmat);
            radmat = radmat(ind_mat);
            Xrad(:,(m-1)/kk+1) = radmat(:);
        end
        Xrad_cell{nn} = Xrad(:);
    end
    Xrad = zeros(prod(Ls),1);
    for nn=1:n
        Xrad(sum(Ls(1:nn))+1:sum(Ls(1:nn+1))) = Xrad_cell{nn};
    end
    [h,e] = histcounts(Xrad,numbins,'normalization','pdf');
    e = (e(1:end-1)+e(2:end))/2;
end