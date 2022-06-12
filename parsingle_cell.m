if i1<=length(ninds)

ind_cell= ninds(i1);

tic;
%%% Build linear system
[Gs,bs,Cfs] = wsindy_sde2nd_fun_nonspatial(X,V,tobs,pt,mt,1,ind_cell,neighbs_cell{expr},tinds,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);    

%%% Solve sparse regression problem
[W,its,lossvals,lambda_hat] = sparsifyDynamics_seq(lambda,gamma*norm(Gs.*(M')),[Gs.*(M');G_append],[bs;b_append],M,Aineq.*(M'),bineq,Aeq.*(M'),beq,excl_inds,const_tol,max_its_qp,disp_opt,max_its_stls,alpha);
resid = (bs-Gs*W)/norm(bs);

%%% Compute learned forces
[f_learned,h_learned,d_learned] = gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);

%%% Compute learned diffusion
nu_learned = rms(Gs*W-bs)^2/2;

%%% Compute neighbors of homing cell
valid_cells = get_neighbs(X, ind_cell, home_cell, V, expr, nearestKLLneighbs, alphaKL, knnp1, tobs, numbins);

%%% Validate model on neighbor cells
[Xspred,errs] = test_neighbs(X,V,tobs,nu_learned,nufac_x,nufac_v,f_learned,h_learned,d_learned,valid_cells,opts,subdt,avg_v0,test_tinds_frac,verbose,1);
disp(['Time to learn single-cell model: ',num2str(toc), 's; ind: ',num2str(i1)])

%%% Save single-cell model
algouttemp={f_learned,h_learned,d_learned,W,resid,lossvals,lambda_hat};

 else

   disp(['Exceeding cell limit'])

end
