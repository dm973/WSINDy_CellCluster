%% Consolidate data

algout=cell(length(ninds),1);
simdat=cell(length(ninds),1);
neighbs=cell(length(ninds),1);
keepinds=[];

for j=1:length(ninds)
    try
        load([homecell_dr,'home_cell_',num2str(j),'.mat'],'algouttemp','Xspred','errs','valid_cells','nu_learned')
        algout{j}=algouttemp;
        simdat{j}{1}=Xspred;
        simdat{j}{2}=errs;
        neighbs{j}=valid_cells(1:knnp1);
        disp([j rms(errs)])
        keepinds=[keepinds j];
disp(['found ind ',num2str(j)]);
    catch
    disp(['missing ',num2str(j)]);
    end
end

algout=algout(keepinds);
simdat=simdat(keepinds);
neighbs=neighbs(keepinds);
ninds=ninds(keepinds);
home_cell{expr}=home_cell{expr}(keepinds);

consol_data = [save_dr,'singlecell_',input_data];
save(consol_data)
