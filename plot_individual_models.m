%% view forces from individual models

clear all

exper = '111';
force_ind = 1; % 1,2,3 = a-r,align,drag
num_cells_per_species = 4;

data_dr = ['~/Desktop/JRSI_data/WSINDy_CellCluster_data/data_dr/',exper,'/'];
save_dr=data_dr;
input_data = findfilestrloc(save_dr,'sim',1);
consol_data = [save_dr,'singlecell_0_',input_data];
load(consol_data,'algout')
load([data_dr,input_data],'inds_cell_true')
species_inds_true=find(cellfun(@(x) ~isempty(x),inds_cell_true));
species_str='ABC';

dims = [num_cells_per_species length(species_inds_true)];
inds = zeros(dims);
rng('shuffle')
for spec_ind=1:length(species_inds_true)
    
    %%% choose cells to display, or random selection
    new_inds = randperm(length(inds_cell_true{species_inds_true(spec_ind)}),num_cells_per_species);
%     new_inds = 1:num_cells_per_species;

    inds(:,spec_ind) = inds_cell_true{species_inds_true(spec_ind)}(new_inds);

end

n = 100;
r = 0.1;
v =[0 1];
[rr,th,xx,yy] = build_polar_grid(n,r,v);
x = xx(1,:);
y = yy(:,1);
c_range = [-15 15];

for j=1:size(inds,2)
    for i=1:size(inds,1)
        f_learned = algout{inds(i,j)}{force_ind};
        if ~isempty(f_learned) 
            F_dat = f_learned(rr,th);
        else
            F_dat = 0*xx;
        end    
        subplot(dims(1),dims(2),(i-1)*length(species_inds_true)+j);
        surfplot;
        title(species_str(species_inds_true(j)),'interpreter','latex')
    end
end
