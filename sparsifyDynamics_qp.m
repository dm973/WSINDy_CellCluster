function [Xi,its,thrs_EL] = sparsifyDynamics_qp(Theta,dXdt,lambda,gamma,M,Aineq,bineq,Aeq,beq,excl_inds,const_tol,max_its,disp_opt,max_its_stls)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% modified by Daniel A. Messenger, 2020 to prevent return of zero vector
% and include regularization
%
% compute Sparse regression: sequential least squares

if ~isempty(Aineq)
    options = optimoptions('quadprog','Display',disp_opt,'ConstraintTolerance',const_tol,'MaxIterations',max_its,'OptimalityTolerance',10^-10);
    n = min(size(dXdt));
    nn = size(Theta,2);

    dXdt_conj = Theta'*dXdt;
    Theta_conj = Theta'*Theta;

    if  gamma ~= 0
%         if ~isempty(M)
%             Theta_conj = Theta_conj+gamma^2*diag(M.^2);
%         else
            Theta_conj = Theta_conj+gamma^2*eye(nn);
%         end
    end

    Xi = quadprog(Theta_conj,-dXdt_conj,Aineq,bineq,Aeq,beq,[],[],[],options);  % initial guess: Least-squares
    if ~isempty(M)
        Xi = M.*Xi;
    end

    if isempty(M)
        thrs_EL = [];
    else
        bnds = norm(dXdt)./vecnorm(Theta)'.*M;
        LBs = lambda*max(1,bnds);
        UBs = 1/lambda*min(1,bnds);
        LBs = lambda*bnds;
        UBs = 1/lambda*bnds;
        thrs_EL = [LBs bnds UBs];
    end

    smallinds = 0*Xi;
    for j=1:max_its_stls
        if ~isempty(M)
            smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
%             smallinds_new = (abs(Xi)./M)<lambda;
            smallinds_new(excl_inds) = 0;
            if or(all(smallinds_new(:)==smallinds(:)),length(find(smallinds_new))==length(Xi))
                its = j;
                return
            else
                smallinds = smallinds_new;
                Xi(smallinds)=0;    
                for ind=1:n
                    biginds = ~smallinds(:,ind);
                    Xi(biginds,ind) = M(biginds).*quadprog(Theta_conj(biginds,biginds),-dXdt_conj(biginds,ind),Aineq(:,biginds),bineq,Aeq(:,biginds),beq,[],[],[],options);
                end
            end
        else
            smallinds_new = (abs(Xi)<lambda);
            smallinds_new(excl_inds) = 0;
            if or(all(smallinds_new(:)==smallinds(:)),length(find(smallinds_new))==length(Xi))
                its = j;
                return
            else
                smallinds = smallinds_new;
                Xi(smallinds)=0;
                for ind = 1:n        
                    biginds = ~smallinds(:,ind);
                    Xi(biginds,ind) = quadprog(Theta_conj(biginds,biginds),-dXdt_conj(biginds,ind),Aineq(:,biginds),bineq,Aeq(:,biginds),beq,[],[],[],options);
                end
            end
        end
    end
    its = nn;
else
    [Xi,its,thrs_EL] = sparsifyDynamics(Theta,dXdt,lambda,gamma,M,max_its_stls);
end
end