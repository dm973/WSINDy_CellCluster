function errs = compute_errs(ninds,simdat,Xscell_obs,Vscell_obs,neighbs,errfun)
    Ni=length(simdat);
    if ~isempty(neighbs)
    
        knnp1= length(simdat{1}{1});
        errs=zeros(Ni,knnp1);
        
        for i=1:Ni
            for j=1:knnp1
                X = simdat{i}{1}{j}{1};
                V = simdat{i}{1}{j}{2};
                M = size(X,3);
                try
                    Y = Xscell_obs{1}(ninds(neighbs{i}(j)),:,1:M);
                    W = Vscell_obs{1}(ninds(neighbs{i}(j)),:,1:M);
                catch
                    Y = Xscell_obs{1}(neighbs{i}(j),:,1:M);
                    W = Vscell_obs{1}(neighbs{i}(j),:,1:M);
    %                 if j==1
    %                     neighbs{i}(j)
    %                 end
                end
                errs(i,j) = errfun(X,Y,V,W);
            end
        end
    else
        
        errs=zeros(Ni,1);
    
        for i=1:Ni
            X = simdat{i}{1};
            V = simdat{i}{2};
            M = size(X,3);
            Y = Xscell_obs{1}(ninds(i),:,1:M);
            W = Vscell_obs{1}(ninds(i),:,1:M);
            errs(i) = errfun(X,Y,V,W);
        end
    end
end