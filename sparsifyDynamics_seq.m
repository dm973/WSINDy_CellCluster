function [W,its_all,lossvals,lambda_hat] = sparsifyDynamics_seq(lambdas,gamma,G,b,M,Aineq,bineq,Aeq,beq,excl_inds,const_tol,max_its_qp,disp_opt,max_its_stls,alpha)

if ~isequal(length(lambdas),1)
    W_ls = [G;gamma*eye(size(G,2))] \ [b;zeros(size(G,2),1)];
    GW_ls = norm(G*W_ls);
    its_all = 1;

    if isempty(lambdas)
        lam_max = max(max(abs(G'*b),[],2)./vecnorm(G).^2');
        lam_min = min(vecnorm(G*W_ls))/size(G,2)/max(vecnorm(G));
        lambdas = 10.^linspace(log10(lam_min), log10(lam_max),50);
    end

    proj_cost = zeros(1,length(lambdas));
    overfit_cost = proj_cost;
    lossvals = proj_cost;

    for l=1:length(lambdas)
        lambda = lambdas(l);
        [W,its,~] = sparsifyDynamics_qp(G,b,lambda,gamma,M,Aineq,bineq,Aeq,beq,excl_inds,const_tol,max_its_qp,disp_opt,max_its_stls);
        if ~isempty(M)
            proj_cost(l) = alpha*norm(G*(W./M-W_ls))/GW_ls;
        else
            proj_cost(l) = alpha*norm(G*(W-W_ls))/GW_ls;
        end
        overfit_cost(l) = (1-alpha)*length(find(W~=0))/length(find(W_ls(:)));
        lossvals(l) = proj_cost(l) + overfit_cost(l);
        its_all = its_all +its;
    end

    l = find(lossvals == min(lossvals),1);

    lambda_hat= lambdas(l);
else
    lambda_hat= lambdas;
    proj_cost=0; overfit_cost=1-alpha; lossvals = 1-alpha;l=1;
    its_all=0;
end

[W,its,~] = sparsifyDynamics_qp(G,b,lambda_hat,gamma,M,Aineq,bineq,Aeq,beq,excl_inds,const_tol,max_its_qp,disp_opt,size(G,2));
its_all = its_all + its;
lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];

end
