exper = '111';
data_dr = ['~/Desktop/data_dr/',exper,'/'];
save_dr = data_dr;
input_data = findfilestrloc(save_dr,'sim',1);
consol_data = [save_dr,'singlecell_',input_data];

load(consol_data,'fx_fcn_cell','fv_fcn_cell','hx_fcn_cell','hv_fcn_cell','dx_fcn_cell','dv_fcn_cell','J_fx','J_fv','J_hx','J_hv','J_dx','J_dv');

f_ar_true = @(r,t) (15+10*cos(2*t)).*(exp(-20*r) - 0.25*exp(-10*r));
f_align_true = @(r,t) -(8+8*cos(t)).*exp(-8*r);
f_drag_true = @(s,t) -5*s;

tpoints = 10^5;
r = 2;
v = [0 1];
lambda = 10^-8;

W_ar = get_true_weights(f_ar_true,tpoints,r,v,fx_fcn_cell,fv_fcn_cell,lambda);
W_align = get_true_weights(f_align_true,tpoints,r,v,hx_fcn_cell,hv_fcn_cell,lambda);
W_drag = get_true_weights(f_drag_true,tpoints,r,v,dx_fcn_cell,dv_fcn_cell,lambda);

sl = [1 1 1];
W_true = [sl(1)*W_ar;sl(2)*W_align;sl(3)*W_drag];
%% compare coefficients

P = [];
for j=1:J_fv
    P = [P j:J_fv:J_fx*J_fv];
end
for j=1:J_hv
    P = [P J_fx*J_fv+(j:J_hv:J_hx*J_hv)];
end
for j=1:J_dv
    P = [P J_fx*J_fv+J_hx*J_hv+(j:J_dv:J_dx*J_dv)];
end

load(consol_data,'algout')

figure(1); clf
bar(W_true(P)); hold on
inds = 1:333;
Ws = [];

for i=inds
    W = algout{i}{4}(P);
    Ws = [Ws W];
	plot(W,'.')
	hold on
end

figure(2); clf
x = 1:length(W_true);

classified_data = [save_dr,'classify_',input_data];
load(classified_data,'species_models');

dat = [W_true(P)'; mean(Ws,2)'; species_models{3}{1}'];
bar(x,dat)
modes = [1 19 37 55 63 71 79 83]-0.5;
for j=1:length(modes)
    hold on
    plot([modes(j) modes(j)],[min(dat(:)) max(dat(:))],'k--')
end
ylim([min(dat(:)) max(dat(:))])
legend({'True','Average A cells'},'location','bestoutside')
