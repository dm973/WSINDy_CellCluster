%% precompute and save quantities for remote single-cell learning
%%%% use save_str

save_str = [data_dr,'precomp_',input_data];

singlecell_inputs;
precomp_learningenvironment;

singlecell_inputs_txt=fileread('singlecell_inputs.m');

ninds = home_cell{expr};
X = Xscell_obs{expr};
V = Vscell_obs{expr};

clear Xscell Vscell Vnormcol
save(save_str)
