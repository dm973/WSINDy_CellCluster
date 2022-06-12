tic;
%% velocity, subsampling
%%% inputs: t, subt, velmeth, subruns, Xscell, Vscell
%%% ouputs: Xscell_obs, Vscell_obs, tobs, runs, numt

runs = length(subruns);
Xscell_obs = cellfun(@(x) x(:,:,1+floor(end*subt(1)):subt(2):floor(end*subt(3))),Xscell(subruns),'uni',0);
tobs = t(1+floor(end*subt(1)):subt(2):floor(end*subt(3)));

if exist('Vscell','var')
    Vscell_obs = cellfun(@(x) x(:,:,1+floor(end*subt(1)):subt(2):floor(end*subt(3))),Vscell(subruns),'uni',0);
else
    Vscell_obs = Xscell_obs;
end

dt = mean(diff(tobs));
if or(isequal(velmeth,'FD'),and(isequal(velmeth,[]),~exist('Vscell','var')))
    ker = [1 0 -1]'/2/dt;
elseif isequal(velmeth,'Gauss')
    sig = 2;
    numsigs = 5;
    kerf = @(x,sig) -x/sig^2.*exp(-(x/sig).^2/2)/sqrt(2*pi*sig^2);
    sigx = dt*(-sig*numsigs:numsigs*sig)';
    ker = kerf(sigx,sig*dt)*dt;
end

for j=1:runs
    if ~isempty(velmeth)
        Vscell_obs{j} = permute(convn(permute(Xscell_obs{j},[3 2 1]),ker,'valid'),[3 2 1]);
        Xscell_obs{j} = Xscell_obs{j}(:,:,(length(ker)-1)/2+1:end-(length(ker)-1)/2);
        tobs = tobs(ceil(length(ker)/2):end-floor(length(ker)/2));
    end
end

X = Xscell_obs{expr};
V = Vscell_obs{expr};

%% Compute distributions: interparticle distance, interparticle relative speed, individual speed 
%%% inputs:  Xscell_obs, Vscell_obs
%%% outputs: xx_grid, hxx, hxx_dist, vv_grid, hvv, hvv_dist, hv, v_grid, vdist

[hxx,xx_grid]  = particle_dist(Xscell_obs,250,10^7);
hxx_dist = cumsum(hxx); hxx_dist = hxx_dist/max(hxx_dist);

[hvv,vv_grid]  = particle_dist(Vscell_obs,250,10^7);
hvv_dist = cumsum(hvv); hvv_dist = hvv_dist/max(hvv_dist);

[hv,v_grid] = normcell(Vscell_obs,250);
v_grid = v_grid(1:end-1);
vdist = cumsum(hv); vdist = vdist/max(vdist);

%% compute neighbors cells: vel_min_n, w_neighb_frac
%%% inputs: vel_min_n, w_neighb_frac, Xscell_obs, Vscell_obs
%%% outputs: w_neighb, neighbs_cell, interact_rad

neighbs_cell = cell(runs,1);
for nn=1:runs
    X = Xscell_obs{nn};
    V = Vscell_obs{nn};
    [N,~,numt] = size(X);
    
    muX=[mean(reshape(X(:,1,:),[],1)) mean(reshape(X(:,2,:),[],1))];
    maxX=[max(reshape(X(:,1,:),[],1)) max(reshape(X(:,2,:),[],1))];
    minX=[min(reshape(X(:,1,:),[],1)) min(reshape(X(:,2,:),[],1))];
    
    w_neighb=[ [muX(1)-w_neighb_frac*(muX(1)-minX(1)) muX(1)+w_neighb_frac*(maxX(1)-muX(1))];... 
        [muX(2)-w_neighb_frac*(muX(2)-minX(2)) muX(2)+w_neighb_frac*(maxX(2)-muX(2))]];

    for tt=1:numt
        neighbs_cell_temp = {};
        neighbs_cell_temp{end+1} = find(and(X(:,1,tt)>w_neighb(1,1),X(:,1,tt)<w_neighb(1,2)));
        neighbs_cell_temp{end+1} = find(and(X(:,2,tt)>w_neighb(2,1),X(:,2,tt)<w_neighb(2,2)));
        neighbs_cell_temp{end+1} = find(min(vecnorm(V(:,:,tt),2,2),[],2)>=vel_min_n);
        neighbs_cell{nn}{tt} = 1:N;
        for i=1:length(neighbs_cell_temp)
            neighbs_cell{nn}{tt} = intersect(neighbs_cell{nn}{tt},neighbs_cell_temp{i});
        end
    end
end

%% compute homing cells
%%% inputs: as, cutoff_n, w_home_frac, interact_rad_frac, std_vel_cap, outliermeth, vel_min_h, Xscell_obs, Vscell_obs
%%% outputs: w_home, home_cell

home_cell = cell(runs,1);
for nn=1:runs
    X = Xscell_obs{nn};
    V = Vscell_obs{nn};
    [N,~,numt] = size(X);
    
    muX=[mean(reshape(X(:,1,:),[],1)) mean(reshape(X(:,2,:),[],1))];
    maxX=[max(reshape(X(:,1,:),[],1)) max(reshape(X(:,2,:),[],1))];
    minX=[min(reshape(X(:,1,:),[],1)) min(reshape(X(:,2,:),[],1))];
    
    w_home=[ [muX(1)-w_home_frac*(muX(1)-minX(1)) muX(1)+w_home_frac*(maxX(1)-muX(1))];... 
        [muX(2)-w_home_frac*(muX(2)-minX(2)) muX(2)+w_home_frac*(maxX(2)-muX(2))]];
    
    interact_rad= interact_rad_frac*min([muX-minX maxX-muX]);

    if ~exist('as','var')
        as{nn} = ones(N,1);
    end
    home_cell_temp = {};    
    home_cell_temp{end+1} = find(as{nn}>=cutoff_n);
    home_cell_temp{end+1} = find(and(min(squeeze(X(:,1,:)),[],2)>w_home(1,1),max(squeeze(X(:,1,:)),[],2)<w_home(1,2)));
    home_cell_temp{end+1} = find(and(min(squeeze(X(:,2,:)),[],2)>w_home(2,1),max(squeeze(X(:,2,:)),[],2)<w_home(2,2)));
    home_cell_temp{end+1} = find(min(squeeze(vecnorm(V,2,2)),[],2)>=vel_min_h); 
    if isequal(outliermeth,'mahal')
        [~,maxV] = max(squeeze(vecnorm(V,2,2)),[],2);
        Vmax = zeros(N,2);
        for j=1:N
            Vmax(j,:) = V(j,:,maxV(j));
        end
        Vstack = reshape(permute(V,[1 3 2]),[],2);
        Vstack=Vstack(vecnorm(Vstack,2,2)>vel_min_h,:);
        datasampV=datasample(Vstack,10^4);
        [stdV,meanV]=robustcov_dam(datasampV);
        stdVinv = inv(stdV);
        Vmaxmu = Vmax-meanV; 
        home_cell_temp{end+1} = find(dot(Vmaxmu,Vmaxmu*stdVinv,2)<std_vel_cap^2);
    elseif isequal(outliermeth,'trad')
        [Vnormmax,~] = max(squeeze(vecnorm(V,2,2)),[],2);
        Vnormcol = reshape(vecnorm(V,2,2),[],1);
        home_cell_temp{end+1} = find(Vnormmax < mean(Vnormcol)+std_vel_cap*std(Vnormcol));
    end
    
    if rand_subsample<inf
        home_cell{nn} = randperm(N,rand_subsample);
    else
        home_cell{nn} = 1:N;
    end
    for i=1:length(home_cell_temp)
        home_cell{nn} = sort(intersect(home_cell{nn},home_cell_temp{i}));
    end
end

%% Compute trial functions
%%% inputs: {xbasis_f, J_fx, rad_cfac_fx, tol_fx_prob (x6 for f,h,d x/v)}, hxx, hvv. hv, xx_grid, vv_grid, v_grid
%%% outputs: fx_fcn_cell (x6)

tol_fx = 1;xx_grid(find(hxx_dist>=tol_fx_prob,1));
tol_hx = 1;xx_grid(find(hxx_dist>=tol_hx_prob,1));
tol_dx = v_grid(find(vdist>=tol_dx_prob,1));

[fx_fcn_cell,fx_tags,radcJ_fx] = build_fj_fcn_cell(J_fx,[hxx;xx_grid],xbasis_f,rad_cfac_fx,tol_fx);
[fv_fcn_cell,fv_tags,~] = build_fj_fcn_cell(J_fv,[hvv;vv_grid],vbasis_f,rad_cfac_fv,tol_fv);

[hx_fcn_cell,~,radcJ_hx] = build_fj_fcn_cell(J_hx,[hxx;xx_grid],xbasis_h,rad_cfac_hx,tol_hx);
[hv_fcn_cell,~,radcJ_hv] = build_fj_fcn_cell(J_hv,[hvv;vv_grid],vbasis_h,rad_cfac_hv,tol_hv);

[dx_fcn_cell,~,~] = build_fj_fcn_cell(J_dx,[hv;v_grid],xbasis_d,rad_cfac_dx,tol_dx);
[dv_fcn_cell,~,~] = build_fj_fcn_cell(J_dv,[hvv;vv_grid],vbasis_d,rad_cfac_dv,tol_dv);

%% Get test functions, allocate space for linear system 
%%% inputs:
%%% ouputs:

mps = cell(runs,1);
for nn=1:runs
    if isequal(tfmeth,'FFT')
        [mps{nn}{1},mps{nn}{2},sig_est,corners] = findcorners_sde(Xscell_obs{nn}(home_cell{nn},:,:),tobs,tau,tauhat,max_dt);
    elseif isequal(tfmeth,'custom')
        mps{nn}{1} = floor((length(tobs)-1)/mtfac);
        mps{nn}{2} = pt;
    end
end

tinds = 1:s:length(tobs)-2*mps{expr}{1};
[mt,pt] = mps{expr}{:};

%% Build constraints

%%% constraints for f
fj_const = [];

% near field repulsion
core_rad = xx_grid(find(hxx_dist>core_rad_prob,1));
core = linspace(eps,core_rad,num_core_pts);
angs_core = linspace(0,pi,2*(J_fv-1)+1);
for i=1:length(angs_core)
    fj_const = [fj_const;-kron(cell2mat(cellfun(@(x)x(core),fx_fcn_cell,'uniformoutput',false))',cellfun(@(x) x(angs_core(i)),fv_fcn_cell)')];
end

% farfield (non-positively) decay 
if interact_rad<xx_grid(end)
    farfield = linspace(interact_rad,xx_grid(end),num_ff_pts);
    angs_ff = linspace(0,pi,2*(J_fv-1)+1);
    for i=1:length(angs_ff)
        fj_const = [fj_const;kron(cell2mat(cellfun(@(x)x(farfield),fx_fcn_cell,'uniformoutput',false))',cellfun(@(x) x(angs_ff(i)),fv_fcn_cell)')];
    end
end

L1 = size(fj_const,1);
b_fj = zeros(L1,1);

%%% constraints for h
fv_const = eye(J_hx*J_hv);
L2 = size(fv_const,1);
b_fv = zeros(L2,1);

%%% constraints for d
fD_const = eye(J_dx*J_dv);
L3 = size(fD_const,1);
b_fD = zeros(L3,1);

Aineq = zeros(L1+L2+L3,J_fx*J_fv+J_hx*J_hv+J_dx*J_dv);
Aineq(1:L1,1:J_fx*J_fv) = fj_const;
Aineq(L1+1:L1+L2,J_fx*J_fv+1:J_fx*J_fv+J_hx*J_hv) = fv_const;
Aineq(L1+L2+1:end,J_fx*J_fv+J_hx*J_hv+1:end) = fD_const; 
bineq = [b_fj;b_fv;b_fD];
% bineq=0;Aineq=Aineq(1,:)*0;
% Aeq = zeros(1,size(Aineq,2)); 
% Aeq(J_fx*J_fv+J_hx*J_hv)=1;
% beq = -1;
% excl_inds =size(J_fx*J_fv+J_hx*J_hv,2);
Aeq = zeros(1,size(Aineq,2)); 
beq = 0;
excl_inds =[];

M = ones(J_fx*J_fv+J_hx*J_hv+J_dx*J_dv,1);
% M(1:J_fx*J_fv) = 1./trapz(xx_grid,hxx.*xx_grid);
% M(J_fx*J_fv+1:J_fx*J_fv+J_hx*J_hv) = 1./trapz(vv_grid,hvv.*vv_grid);
% M(J_fx*J_fv+J_hx*J_hv+1:end) = 1./trapz(v_grid,hv.*v_grid);

G_append = [];
b_append = [];