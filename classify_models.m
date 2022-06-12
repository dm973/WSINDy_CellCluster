classify_inputs_txt=fileread('classify_inputs.m');

J_fx = length(fx_fcn_cell);
J_fv = length(fv_fcn_cell);
J_hx = length(hx_fcn_cell);
J_hv = length(hv_fcn_cell);
J_dx = length(dx_fcn_cell);
J_dv = length(dv_fcn_cell);

test_tinds = 1:floor(length(tobs)*test_tinds_frac);

nux = nu_learned*nufac_x;
nuv = nu_learned*nufac_v;
Xtest=Xscell_obs{1}(:,:,test_tinds);
Vtest=Vscell_obs{1}(:,:,test_tinds);
tobs_test=tobs(test_tinds);
[N,d,M]=size(Xtest);

species_inds={};
species_models={};
valpairs_cell={};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% learn models

subinds=1:length(ninds);

while ~isempty(subinds)

    %%% cluster models
    ge = compute_errs(ninds(subinds),simdat(subinds),Xscell_obs,Vscell_obs,neighbs(subinds),errfun);
    [Wsmat,Mod,inds_keep,ge_score,models,inds_pat,ninds_rep]=...
        clustermodels(alpha,beta,gamma,normGE,normALL,tol1,tol2,ninds(subinds),simdat(subinds),neighbs(subinds),algout(subinds),J_fv,J_fx,J_hv,J_hx,J_dv,J_dx,ge);

    %%% compute top model
    if ~isempty(inds_pat)
        W=mean(Wsmat(:,inds_pat{1}{1}),2);
        %%% could use median here
        W(log10(abs(W))<max(log10(abs(W)))-logcutoff)=0;

        [f_learned,h_learned,d_learned] = gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);

        %%% get cells to validate
        valid_cells=ninds(subinds);
        num_gen=length(valid_cells);
        
        %%% Validate model on remaining cells
        [valpairs,~] = test_neighbs(Xtest,Vtest,tobs,nu_learned,nufac_x,nufac_v,f_learned,h_learned,d_learned,valid_cells,opts,subdt,avg_v0,test_tinds_frac,2,1);

        valpairs_cell=[valpairs_cell,{valpairs}];
        
        %%% compute log val. errs
        relerr = compute_errs(valid_cells,valpairs,{Xtest},{Vtest},[],errfun);
        frelerr = log10(relerr);
        
        %%% if >100(halt_prob)% of cells have error > 100(accept_err)%, check for 2-species
        if all([sum(relerr>accept_err)/length(relerr)>halt_prob length(species_inds)<max_species length(relerr)>halt_num])
            %%% fit to gaussian 
            gm1=fitgmdist(frelerr(:),1);
            idk1 = cluster(gm1,frelerr(:));

            %%% if 2-guassian mixture has lower mean BIC, 2-species, if not,
            %%% forms single species.        
            idks=repmat(idk1*0,1,num_gm_tries);
            bic = zeros(1,num_gm_tries);
            for j=1:num_gm_tries
                gm=fitgmdist(frelerr(:),2,'RegularizationValue',10^-6);
                idk = cluster(gm,frelerr(:));
                bic(j)=gm.BIC;
                [~,ii]=min(gm.mu);
                idks(:,j)= idk==ii;
            end
            if mean(bic) < gm1.BIC
                disp('multi-species')
                subinds=subinds(mean(idks,2)>0.5);
            else
                disp('mono-species')
            end
            species_inds=[species_inds,{subinds}];
            species_models=[species_models,{{W,f_learned,h_learned,d_learned}}];
            subinds = oppinds(cell2mat(species_inds),length(ninds));
        else
            low_err_inds = find(relerr<accept_err)';
            if ~isempty(low_err_inds)
                species_inds=[species_inds,{low_err_inds}];
                species_models=[species_models,{{W,f_learned,h_learned,d_learned}}];
            end
            species_inds=[species_inds,{oppinds(cell2mat(species_inds),length(ninds))}];
            species_models=[species_models,{{}}];
            subinds = [];
        end
    else
        species_inds=[species_inds,{oppinds(cell2mat(species_inds),length(ninds))}];
        species_models=[species_models,{{}}];
        subinds = [];
    end            
end