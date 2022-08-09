%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set params

%%% threshold coefficients in aggregate model that are outside of the top
%%% 'logcutoff' orders of magnitude among all aggregate model coefficients

logcutoff = 4;

%%% define validation error metrics 

alpha=1/2; % coefficient of VE
beta=alpha; % coefficient of RES
gamma=0; % coefficient of GE
normGE=0; % norm for vector of GE
normALL=1; % norm for components of p-score
tol1=2/3; % accept neighbor model if error reduced by at least a factor of tol1
tol2=0.25; % accept neighbor model if neighbors p-score at most tol2, 
           % include models in average if p-score at most tol2
           % include class of models if at least one member has p-score at most tol2

%%% define error function 

tcap=1;
errfun = @(X,Y,V,W) norm(reshape(V(:,:,1:floor(tcap*end))-W(:,:,1:floor(tcap*end)),[],1))/...
    norm(reshape(W(:,:,1:floor(tcap*end)),[],1));
%     errfun = @(X,Y,V,W) norm(reshape(X(:,:,1:floor(tcap*end))-Y(:,:,1:floor(tcap*end)),[],1))/...
%         norm(reshape(Y(:,:,1:floor(tcap*end))-Y(:,:,1),[],1));
%     errfun = @(X,Y,V,W) sum(vecnorm(X(:,:,1:floor(tcap*end))-Y(:,:,1:floor(tcap*end)),2,2)>0.025)/size(Y(:,:,1:floor(tcap*end)),3);

%%% define validation simulation parameters 

test_tinds_frac=0.25;
nu_learned = 0;
subdt = 32;
nufac_x = 0;
nufac_v = 0;
avg_v0 = 1;
expr = 1;
opts = [1 0 0];
verbose =0;

%%% define halting criteria

num_gm_tries = 20;
max_species = 10;
accept_err = 0.05;
halt_prob = 0.01;
halt_num = 2;