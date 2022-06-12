load(single_cell_data)
classify_inputs;
classify_models;
classified_data = [data_dr,'classify_',input_data];
save(classified_data,'species_inds','species_models','valpairs_cell');
