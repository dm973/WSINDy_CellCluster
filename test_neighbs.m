function [Xspred,errs,nux,nuv] = test_neighbs(X,V,tobs,nu_learned,nufac_x,nufac_v,f_learned,h_learned,d_learned,valid_cells,opts,subdt,avg_v0,test_tinds_frac,verbose,toggle_par)

nux = nu_learned*nufac_x;
nuv = nu_learned*nufac_v;

tcap=1;
errfun = @(X,Y,V,W) norm(reshape(V(:,:,1:floor(tcap*end))-W(:,:,1:floor(tcap*end)),[],1))/...
    norm(reshape(W(:,:,1:floor(tcap*end)),[],1));

if ~isempty(f_learned)
    f_xv_pred = @(x,v,t) f_learned(x,v);
else
    f_xv_pred = [];
end
if ~isempty(h_learned)
    h_xv_pred = @(x,v,t) h_learned(x,v);
else
    h_xv_pred = [];
end
if ~isempty(d_learned)
    d_xv_pred = @(x,v,t) d_learned(x,v);
else
    d_xv_pred = [];
end

tinds_test = 1:floor(test_tinds_frac*length(tobs));
ttest = tobs(tinds_test);
X = X(:,:,tinds_test);
V = V(:,:,tinds_test);

Xspred=cell(length(valid_cells),1);
errs=zeros(length(valid_cells),1);
[N,d,~] = size(X);
RHS = @(Z,t) rhs_lean_pred(f_xv_pred, h_xv_pred, d_xv_pred, Z, d, t, opts,1);

if length(valid_cells) == 1
    toggle_par = 0;
end

if toggle_par

    parfor jj=1:length(valid_cells)
        if verbose>1
            tic;
        end
        inds_n = valid_cells(jj);
        inds_np = find(~ismember(1:N,inds_n));
        Xtemp = X(inds_n,:,1);
        Vtemp = mean(V(inds_n,:,1:avg_v0),3);
        X_off = X(inds_np,:,:);
        V_off = V(inds_np,:,:);
        Ztemp = reshape([Xtemp Vtemp],[],1);
        [~,Z,B] = simnu_split(RHS,ttest,Ztemp,subdt,nux,nuv,0,X_off,V_off);
        X_pred = zeros(1,d,length(ttest));
        V_pred = zeros(1,d,length(ttest));
        for m=1:length(ttest)
            X_pred(:,:,m) = reshape(Z(m,1:end/2)',1,d);
            V_pred(:,:,m) = reshape(Z(m,end/2+1:end)',1,d);
        end
        Xspred{jj}={X_pred,V_pred};
        errs(jj) = errfun(X_pred,X(inds_n,:,:),V_pred,V(inds_n,:,:));

        if verbose>1
            disp([toc inds_n jj errs(jj)]);
        end
    end

else

    
    for jj=1:length(valid_cells)
        if verbose>1
            tic;
        end
        inds_n = valid_cells(jj);
        inds_np = find(~ismember(1:N,inds_n));
        Xtemp = X(inds_n,:,1);
        Vtemp = mean(V(inds_n,:,1:avg_v0),3);
        X_off = X(inds_np,:,:);
        V_off = V(inds_np,:,:);
        Ztemp = reshape([Xtemp Vtemp],[],1);
        [~,Z,B] = simnu_split(RHS,ttest,Ztemp,subdt,nux,nuv,0,X_off,V_off);
        X_pred = zeros(1,d,length(ttest));
        V_pred = zeros(1,d,length(ttest));
        for m=1:length(ttest)
            X_pred(:,:,m) = reshape(Z(m,1:end/2)',1,d);
            V_pred(:,:,m) = reshape(Z(m,end/2+1:end)',1,d);
        end
        Xspred{jj}={X_pred,V_pred};
        errs(jj) = errfun(X_pred,X(inds_n,:,:),V_pred,V(inds_n,:,:));

        if verbose>1
            disp([toc inds_n jj errs(jj)]);
            plot(squeeze(X(inds_n,1,:)),squeeze(X(inds_n,2,:)),'k',squeeze(X_pred(1,1,:)),squeeze(X_pred(1,2,:)),'r--','linewidth',2)
            drawnow
        end
    end

end

if verbose
    disp(['Errors on neighbor cells: '])
    disp(errs(:)');
end

end