%%% data sampling
subruns = 1;
expr = 1;
subt = [0 10 1];

%%% velocity method
velmeth='FD';

%%% neighbor cell sampling
vel_min_n = 0;
w_neighb_frac = inf;

%%% homing cell samping
rand_subsample = 100; %%% set to inf to use all cells, otherwise a random subset is chosen
cutoff_n = 0.999;
w_home_frac = inf;
interact_rad_frac = inf;
std_vel_cap = inf;
outliermeth='trad';
vel_min_h = 0;

%%% trial bases
xbasis_f = 'laguerre';
J_fx = 18;
rad_cfac_fx = 1;
tol_fx_prob = 0.36;

vbasis_f = 'cos2';
J_fv = 3;
rad_cfac_fv = [];
tol_fv = [];

xbasis_h = 'MO';
J_hx = 8;
rad_cfac_hx = 2.^(-3:5);
tol_hx_prob = tol_fx_prob;

vbasis_h = 'cos';
J_hv = 3;
rad_cfac_hv = [];
tol_hv = [];

xbasis_d = 'mon';
J_dx = 4;
rad_cfac_dx = 1;
tol_dx_prob = 0.999;

vbasis_d = 'cos';
J_dv = 2;
rad_cfac_dv = [];
tol_dv = [];

%%% temporal test function 
max_dt=2; s=1;
tau = 10^-10;tauhat = 3;
mtfac=25; pt=20;
tfmeth='FFT';

%%% regression parameters
lambda = 10.^(linspace(-4,-0.5,8*8));
gamma = 10^-inf;
delta = 0;
max_its_stls = 10;
alpha = 0.5;
const_tol = 10^-14; 
max_its_qp = 2000;
disp_opt = 'none';

%%%% f constraints
core_rad_prob=0.001;
num_core_pts = 5;
num_ff_pts = 10;

%%% sim params
test_tinds_frac=0.25;
subdt = 32;
nufac_x = 0;
nufac_v = 0;
avg_v0 = 1;
knnp1 = 8;
numbins = 50;
nearestKLLneighbs=200;
alphaKL=[1 1 1];
opts = [1 0 0];
verbose=1;