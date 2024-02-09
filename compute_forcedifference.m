%% choose experiment and load variables
clear all

exper = '111';

data_dr = ['~/Desktop/data_dr/',exper,'/'];
save_dr = data_dr;
input_data = findfilestrloc(save_dr,'sim',1);
load([save_dr,findfile(save_dr,'classify_',[])]);
load([save_dr,findfile(save_dr,'singlecell_',[])],'algout','neighbs','Xscell_obs','Vscell_obs','fx_fcn_cell','fv_fcn_cell','hx_fcn_cell','hv_fcn_cell','dx_fcn_cell','dv_fcn_cell');
load([save_dr,input_data],'f_xv_true_cell','h_xv_true_cell','d_xv_true_cell');
ftrue_cell = {f_xv_true_cell,h_xv_true_cell,d_xv_true_cell};

%% choose cluster and force to view

for cluster_num = 1:length(species_models);
for force_ind = 1:3;

inds = species_inds{cluster_num};
Mod = species_indsInModels{cluster_num}{3};

%%% set plotting arguments

n = 1000;
if force_ind<3
    r = 2;
    v =[0 1];
    [rr,th,xx,yy] = build_polar_grid(n,r,v);
else
    r = 0.05;
    [rr,th] = meshgrid(linspace(0,r,n),linspace(0,pi,n));
    xx=rr;
    yy=th;
end
x = xx(1,:);
y = yy(:,1);

%%% compute 'true' force, associating cluster with its majority species

forceDNA = cellfun(@(x) ~isempty(x), species_models{cluster_num}(2:4));
if isequal(forceDNA,[1 1 1])
    learned_spec = 1;
elseif isequal(forceDNA,[1 0 1])
    learned_spec = 2;
elseif isequal(forceDNA,[0 1 1])
    learned_spec = 3;
end

F_true = ftrue_cell{force_ind}{learned_spec};
if isempty(F_true)
    F_true = rr*0;
else
    F_true = F_true(rr,th,0);
end

%%% compute learned force

W = species_models{cluster_num}{1};

flearned_cell = cell(3,1);
[flearned_cell{:}] = gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);
f_learned = flearned_cell{force_ind};
if ~isempty(f_learned) 
    F_cube_Mod_thresh = f_learned(rr,th);
    disp(exper)
    disp([cluster_num force_ind])
    norm(F_cube_Mod_thresh(:)-F_true(:))/norm(F_true(:))
end

end
end