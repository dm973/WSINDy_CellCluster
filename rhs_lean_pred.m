function dZ = rhs_lean_pred(f_xv_true, h_xv_true, d_xv_true, Ztemp, d, t, opts,n)

    if isempty(opts)
        opts = zeros(3,1);
    end

    Ztemp = reshape(Ztemp,[],2*d);
    agents = Ztemp(:,1:2);
    vels = Ztemp(:,3:4);
        
    N = size(Ztemp,1);
    if ~exist('n','var')
        n = [];
    end
    if isempty(n)
        n = N;
    end
    dX = reshape(agents(1:n,:),n,1,d) - reshape(agents,1,N,d);
    rr = squeeze(vecnorm(dX,2,3));

    dV = repmat(reshape(vels(1:n,:),n,1,d),1,N,1);
    vv = dot(dV,-dX,3);
    vv = atan2(squeeze(dV(:,:,2).*dX(:,:,1)-dV(:,:,1).*dX(:,:,2)),vv);

    updateV = vels(1:n,:)*0;
    
    if ~isempty(f_xv_true)
        forces = f_xv_true(rr,vv,t);
        if opts(1)==1
            updateV = updateV + reshape(mean((forces./max(squeeze(vecnorm(dX,2,3)),eps)).*dX,2),n,d);
        else
            updateV = updateV + reshape(mean(dX.*forces,2),n,d);
        end
    end
    
    if ~isempty(h_xv_true)
        forces = h_xv_true(rr,vv,t);
        dV = reshape(vels(1:n,:),n,1,d) - reshape(vels,1,N,d);
        if opts(2)==1
            updateV = updateV + reshape(mean((forces./max(squeeze(vecnorm(dV,2,3)),eps)).*dV,2),n,d);
        else
            updateV = updateV + reshape(mean(forces.*dV,2),n,d);
        end
    end
    
    if ~isempty(d_xv_true)
        rr = repmat(vecnorm(vels(1:n,:),2,2),1,N);
        forces = d_xv_true(rr,vv,t);
        dV = reshape(vels(1:n,:),[],1,d);
        if opts(3)==1
            updateV = updateV + reshape(mean(forces.*dV./max(squeeze(vecnorm(vels,2,3)),eps),2),n,d);
        else
            updateV = updateV + reshape(mean(forces.*dV,2),n,d);
        end
    end
    
    dZ = reshape([vels [updateV;zeros(N-n,d)]],[],1);
    
end
