function neighbs = get_neighbs(X, V, ind_cell, ninds, nearestKLLneighbs, alphaKL, knnp1, tobs, numbins)
    dist_mat=squeeze(vecnorm(X(ind_cell,:,1)-X(ninds,:,1),2,2));
    [~,neighbs]=sort(min(dist_mat,[],2));
    i2 = ninds(neighbs(1:min(nearestKLLneighbs,end)));
    
    Rcap=inf;
    tpts=1:length(tobs);
    [~,KLxx,~,~]=crossKL(X,X,ind_cell,i2,tpts,numbins,Rcap);
    [~,KLvv,~,~]=crossKL(V,V,ind_cell,i2,tpts,numbins,Rcap);
    [~,KLv,~,~]=crossvel(V,V,ind_cell,i2,tpts,numbins,Rcap);
    [~,Kallinds]=sort(alphaKL(1)*KLxx.^2+alphaKL(2)*KLvv.^2+alphaKL(3)*KLv.^2);
    neighbs=i2(Kallinds(1:min(knnp1,end)));
end