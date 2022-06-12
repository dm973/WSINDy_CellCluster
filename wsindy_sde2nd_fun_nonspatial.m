function [G,b,Cfs_t] = wsindy_sde2nd_fun_nonspatial(X,V,t,pt,mt,phi_t_class,ninds,neighbinds,tinds,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell)

Lt = length(tinds);
J_fx = length(fx_fcn_cell);
J_fv = length(fv_fcn_cell);
J_hx = length(hx_fcn_cell);
J_hv = length(hv_fcn_cell);
J_dx = length(dx_fcn_cell);
J_dv = length(dv_fcn_cell);

[N,d,NtV] = size(V);
if isempty(ninds)
    ninds = 1:N;
end
Ni = length(ninds);
if isempty(neighbinds)
    neighbinds=repmat({(1:N)'},1,NtV);
end

dt = t(2)-t(1);
Cfs_t = phi_int_weights(mt,2,-pt,phi_t_class);
norml = norm(Cfs_t(1,:),2)*sqrt(dt);
Cfs_t = Cfs_t/norml;

b = zeros(Lt*Ni*d,1);
for dd = 1:d
    Btemp = conv2(1,Cfs_t(3,:)/(mt*dt)^2,reshape(X(ninds,dd,:),length(ninds),NtV),'valid');
    b(Lt*Ni*(dd-1)+(1:Lt*Ni))  = reshape(Btemp(:,tinds)*dt,[],1);
end

f_mat_cell = repmat({zeros(Ni,d,NtV)},J_fx*J_fv,1);
h_mat_cell = repmat({zeros(Ni,d,NtV)},J_hx*J_hv,1);
d_mat_cell = repmat({zeros(Ni,d,NtV)},J_dx*J_dv,1);

for tt = 1:NtV
    agents = squeeze(X(neighbinds{tt},:,tt));
    vels = squeeze(V(neighbinds{tt},:,tt));
    agents_sub = squeeze(X(ninds,:,tt));
    vels_sub =  squeeze(V(ninds,:,tt));
    vels = vels(agents(:,1)~=0,:);
    agents = agents(agents(:,1)~=0,:);
    N = size(agents,1);

    dX = reshape(agents_sub,Ni,1,d) - reshape(agents,1,N,d);
    rr = reshape(vecnorm(dX,2,3),Ni,N);
    vnorms = vecnorm(vels_sub,2,2);  
%     vdist_mat = max(squeeze(vecnorm(Vdiffs,2,3)),eps);
    vxdist_mat = dot(repmat(reshape(vels_sub,Ni,1,d),1,N,1),-dX,3);
    vxdist_mat = atan2(abs(dot(repmat(reshape(vels_sub,Ni,1,d),1,N,1),-flip(cat(3,-dX(:,:,1),dX(:,:,2)),3),3)),vxdist_mat);
        
    fx_temp = cellfun(@(x)x(rr),fx_fcn_cell,'uniformoutput',0);
    fv_temp = cellfun(@(x)x(vxdist_mat),fv_fcn_cell,'uniformoutput',0);
    hx_temp = cellfun(@(x)x(rr),hx_fcn_cell,'uniformoutput',0);
    hv_temp = cellfun(@(x)x(vxdist_mat),hv_fcn_cell,'uniformoutput',0);
    dx_temp = cellfun(@(x)x(vnorms),dx_fcn_cell,'uniformoutput',0);
    dv_temp = cellfun(@(x)x(vxdist_mat),dv_fcn_cell,'uniformoutput',0);
    
    for i=1:J_fx 
        for j=1:J_fv
%         fj_mat_cell{j}(:,:,tt) = squeeze(sum(fj_fcn_cell{j}(dist_mat)./dist_denom.*Xdiffs,2));
            f_mat_cell{(i-1)*J_fv+j}(:,:,tt) = squeeze(mean(fx_temp{i}.*fv_temp{j}.*(dX./max(rr,eps)),2));
%        fj_mat_cell{j}(:,:,tt) = squeeze(mean(fj_fcn_cell{j}(dist_mat).*Xdiffs,2));
        end
    end
    for i=1:J_hx 
        for j=1:J_hv
%         fj_mat_cell{j}(:,:,tt) = squeeze(sum(fj_fcn_cell{j}(dist_mat)./dist_denom.*Xdiffs,2));
            h_mat_cell{(i-1)*J_hv+j}(:,:,tt) = squeeze(mean(hx_temp{i}.*hv_temp{j}.*(reshape(vels_sub,Ni,1,d) - reshape(vels,1,N,d)),2));
%        fj_mat_cell{j}(:,:,tt) = squeeze(mean(fj_fcn_cell{j}(dist_mat).*Xdiffs,2));
        end
    end
    for i=1:J_dx 
        for j=1:J_dv
%         fj_mat_cell{j}(:,:,tt) = squeeze(sum(fj_fcn_cell{j}(dist_mat)./dist_denom.*Xdiffs,2));
            d_mat_cell{(i-1)*J_dv+j}(:,:,tt) = squeeze(mean(dv_temp{j},2)).*dx_temp{i}.*(vels_sub);
%        fj_mat_cell{j}(:,:,tt) = squeeze(mean(fj_fcn_cell{j}(dist_mat).*Xdiffs,2));
        end
    end
%     for i=1:J_hx 
%         h_mat_cell{i}(:,:,tt) = squeeze(mean(hx_temp{i}.*Vdiffs,2));
% %            fv_mat_cell{j}(:,:,tt) = squeeze(mean(fv_fcn_cell{j}(dist_mat)./max(vecnorm(Vdiffs,2,3),eps).*Vdiffs,2));
%     end
end

% d_mat_cell = cellfun(@(x)x(Vnrms(ninds,1,:)).*V(ninds,:,:),dx_fcn_cell,'uniformoutput',0);

G = zeros(Lt*Ni*d,J_fx*J_fv+J_hx*J_hv+J_dx*J_dv);
for dd=1:d
    for j=1:J_fx*J_fv
        FJ_temp = conv2(1,Cfs_t(1,:),reshape(f_mat_cell{j}(:,dd,:),Ni,[]),'valid');
        G(Lt*Ni*(dd-1)+(1:Lt*Ni),j) = reshape(FJ_temp(:,tinds)*dt,[],1);
    end
    for j=J_fx*J_fv+1:J_fx*J_fv+J_hx*J_hv
        FJ_temp = conv2(1,Cfs_t(1,:),reshape(h_mat_cell{j-J_fx*J_fv}(:,dd,:),Ni,[]),'valid');
        G(Lt*Ni*(dd-1)+(1:Lt*Ni),j) = reshape(FJ_temp(:,tinds)*dt,[],1);
    end
    for j=J_fx*J_fv+J_hx*J_hv+1:J_fx*J_fv+J_hx*J_hv
        FJ_temp = conv2(1,Cfs_t(1,:),reshape(V_mat_cell{j-J_fx*J_fv-J_hx*J_hv}(:,dd,:),Ni,[]),'valid');
        G(Lt*Ni*(dd-1)+(1:Lt*Ni),j) = reshape(FJ_temp(:,tinds)*dt,[],1);
    end
    for j=J_fx*J_fv+J_hx*J_hv+1:J_fx*J_fv+J_hx*J_hv+J_dx*J_dv
        FJ_temp = conv2(1,Cfs_t(1,:),reshape(d_mat_cell{j-J_fx*J_fv-J_hx*J_hv}(:,dd,:),Ni,[]),'valid');
        G(Lt*Ni*(dd-1)+(1:Lt*Ni),j) = reshape(FJ_temp(:,tinds)*dt,[],1);
    end
end

end