%% Generate trajectories from learned models

specs=1:sum(cellfun(@(x) ~isempty(x),species_models));
X_preds = cell(length(specs),1);

for spec = 1:length(specs)

f_learned = species_models{specs(spec)}{2};
h_learned = species_models{specs(spec)}{3};
d_learned = species_models{specs(spec)}{4};

[~,a]=sort(cellfun(@(x) x{2}(1),simdat(species_inds{spec})));
a=a(min(2,end));
M=dist(squeeze(Xscell_obs{1}(ninds(species_inds{specs(spec)}),:,1))');                                
    [b,cell_inds]=sort(M(a,:));                                                                 
    valid_cells = ninds(species_inds{specs(spec)}(cell_inds));
    [X_pred,~] = test_neighbs(Xscell_obs{1},Vscell_obs{1},tobs,nu_learned,nufac_x,nufac_v,f_learned,h_learned,d_learned,valid_cells,opts,subdt,avg_v0,1,2,1);        
    X_preds{spec} = {X_pred,valid_cells};                                                                                                          
end                                                                                                       
save([data_dr,'learnedtraj_',input_data])

