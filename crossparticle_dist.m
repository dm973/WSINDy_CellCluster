function h  = crossparticle_dist(X,Y,e)
    [Nx,~,M]=size(X);
    [Ny,~,~]=size(Y);
    logmat=~tril(ones(Nx,Ny));
    Xrad=[];
    for j=1:M
        Xt=squeeze(X(:,:,j));
        Yt=squeeze(Y(:,:,j));
        difft=squeeze(vecnorm(permute(Xt,[1 3 2])-permute(Yt,[3 1 2]),2,3));
        Xrad=[Xrad;reshape(difft(logmat),[],1)];
    end
    h = histcounts(Xrad,e,'normalization','pdf');
end
