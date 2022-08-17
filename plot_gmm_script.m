%% choose experiment and load variables\

clear all; close all

exper = '111';
data_dr = ['~/Desktop/data_dr/',exper,'/'];
save_dr = data_dr;

toggle_save = 1;

for spec = 1:5

    plot_gmm;

end

