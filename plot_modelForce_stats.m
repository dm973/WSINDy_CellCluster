%% choose experiment and load variables
clear all

exper = '101';

data_dr = ['~/Desktop/data_dr/',exper,'/'];
save_dr = data_dr;
input_data = findfilestrloc(save_dr,'sim',1);
load([save_dr,findfile(save_dr,'classify_',[])]);
load([save_dr,findfile(save_dr,'singlecell_',[])],'algout','neighbs','Xscell_obs','Vscell_obs','fx_fcn_cell','fv_fcn_cell','hx_fcn_cell','hv_fcn_cell','dx_fcn_cell','dv_fcn_cell');
load([save_dr,input_data],'f_xv_true_cell','h_xv_true_cell','d_xv_true_cell');
ftrue_cell = {f_xv_true_cell,h_xv_true_cell,d_xv_true_cell};

%% choose cluster and force to view

cluster_num = 2;
force_ind = 2;

inds = species_inds{cluster_num};
Mod = species_indsInModels{cluster_num}{3};

%% set plotting arguments

n = 1000;
if force_ind<3
    r = 2;
    v =[0 1];
    [rr,th,xx,yy] = build_polar_grid(n,r,v);
else
    r = 0.2;
    [rr,th] = meshgrid(linspace(0,r,n),linspace(0,pi,n));
    xx=rr;
    yy=th;
end
x = xx(1,:);
y = yy(:,1);

%% compute 'true' force, associating cluster with its majority species

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

%% compute force from raw cluster average 

F_cube = repmat(rr*0,1,1,length(inds));
for i=1:length(inds)
    f_learned = algout{inds(i)}{force_ind};
    if ~isempty(f_learned) 
        F_cube(:,:,i) = f_learned(rr,th);
    end
end
F_mean = mean(F_cube,3);
F_std = std(F_cube-F_true,[],3);

%% compute force from model replacement average + threshold

Wmat = species_indsInModels{cluster_num}{1};
F_cube_Mod_thresh = repmat(rr*0,1,1,size(Wmat,2));
for i=1:size(Wmat,2)
    W = Wmat(:,i);
    W(species_models{cluster_num}{1}==0)=0;
    flearned_cell = cell(3,1);
    [flearned_cell{:}] = gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);
    f_learned = flearned_cell{force_ind};
    if ~isempty(f_learned) 
        F_cube_Mod_thresh(:,:,i) = f_learned(rr,th);
    end
end
F_mean_Mod_thresh = mean(F_cube_Mod_thresh,3);
F_std_Mod_thresh = std(F_cube_Mod_thresh-F_true,[],3);

%% view as subplots

toggle_save_plots = '';%~/Desktop/';

cmap = jet(100);
if force_ind<3
    c_range_mean = [-15 15];
else
    c_range_mean = [-inf inf];
end
c_range_diff = [-inf inf];
c_range_std = [-inf inf];

subplot(2,3,1);
c_range = c_range_mean;
F_dat = F_mean; surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)

subplot(2,3,2);
c_range = c_range_diff;
F_dat = F_mean-F_true; surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)

subplot(2,3,3);
c_range = c_range_std;
F_dat = F_std; surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)

subplot(2,3,4);
c_range = c_range_mean;
F_dat = F_mean_Mod_thresh;surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)

subplot(2,3,5);
c_range = c_range_diff;
F_dat = F_mean_Mod_thresh-F_true;surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)

subplot(2,3,6);
c_range = c_range_std;
F_dat = F_std_Mod_thresh;surfplot;
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)

if ~isempty(toggle_save_plots)

    c_range = c_range_mean;
    figure(1); close; figure(1); 
    F_dat = F_mean; surfplot;
    saveas(gcf,[toggle_save_plots,'Fmean_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
    figure(2); close; figure(2); 
    F_dat = F_mean_Mod_thresh;surfplot;
    saveas(gcf,[toggle_save_plots,'FmeanModcut_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
    
    c_range = c_range_diff;
    figure(1); close; figure(1); 
    F_dat = F_mean-F_true; surfplot;
    saveas(gcf,[toggle_save_plots,'Fdiff_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
    figure(2); close; figure(2); 
    F_dat = F_mean_Mod_thresh-F_true;surfplot;
    saveas(gcf,[toggle_save_plots,'FdiffModcut_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
    
    c_range = c_range_std;
    figure(1); close; figure(1); 
    F_dat = F_std; surfplot;
    saveas(gcf,[toggle_save_plots,'Fstd_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
    figure(2); close; figure(2); 
    F_dat = F_std_Mod_thresh;surfplot;
    saveas(gcf,[toggle_save_plots,'FstdModcut_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])

end

%%%%%%%%%%
% F_cube = repmat(rr*0,1,1,length(inds));
% F_cube_Mod = repmat(rr*0,1,1,length(inds));
% for i=1:length(inds)
%     f_learned = algout{Mod(inds(i))}{force_ind};
%     try 
%         F_cube_Mod(:,:,i) = f_learned(rr,th);
%     catch
%     end
%     f_learned = algout{inds(i)}{force_ind};
%     try 
%         F_cube(:,:,i) = f_learned(rr,th);
%     catch
%     end
% end
% Wmat = species_indsInModels{cluster_num}{1};
% F_cube_Mod_thresh = repmat(rr*0,1,1,size(Wmat,2));
% for i=1:size(Wmat,2)
%     W = Wmat(:,i);
%     W(species_models{cluster_num}{1}==0)=0;
%     flearned_cell = cell(3,1);
%     [flearned_cell{:}]=gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell);
%     f_learned = flearned_cell{force_ind};
%     try 
%         F_cube_Mod_thresh(:,:,i) = f_learned(rr,th);
%     catch
%     end
% end
% F_mean = mean(F_cube,3);
% F_std = std(F_cube-F_true,[],3);
% F_mean_Mod = mean(F_cube_Mod,3);
% F_std_Mod = std(F_cube_Mod-F_true,[],3);
% F_mean_Mod_thresh = mean(F_cube_Mod_thresh,3);
% F_std_Mod_thresh = std(F_cube_Mod_thresh-F_true,[],3);
% figure(2); close; figure(2); 
% F_dat = F_mean_Mod;surfplot;
% saveas(gcf,['~/Desktop/FmeanMod_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
% figure(2); close; figure(2); 
% F_dat = F_mean_Mod-F_true;surfplot;
% saveas(gcf,['~/Desktop/FdiffMod_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
% figure(2); close; figure(2); 
% F_dat = F_std_Mod;surfplot;
% saveas(gcf,['~/Desktop/FstdMod_cluster',num2str(cluster_num)','_force',num2str(force_ind),'.png'])
% 
% subplot(3,3,4);
% surf(xx,yy,F_mean_Mod,'edgeColor','none')
% colorbar
% view([0 90])
% xlabel('x')
% ylabel('y')
% axis equal
% ylim([min(y) max(y)])
% xlim([min(x) max(x)])
% colorbar('ticklabelinterpreter','latex','fontsize',14)
% colormap(cmap)
% set(gca,'ticklabelinterpreter','latex','fontsize',14)
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)
% caxis(c_range)
% 
% subplot(3,3,5);
% surf(xx,yy,F_mean_Mod-F_true,'edgeColor','none')
% colorbar
% view([0 90])
% xlabel('x')
% ylabel('y')
% axis equal
% ylim([min(y) max(y)])
% xlim([min(x) max(x)])
% colorbar('ticklabelinterpreter','latex','fontsize',14)
% colormap(cmap)
% set(gca,'ticklabelinterpreter','latex','fontsize',14)
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)
% caxis(c_range_diff)
% 
% subplot(3,3,6);
% surf(xx,yy,F_std_Mod,'edgeColor','none')
% colorbar
% view([0 90])
% xlabel('x')
% ylabel('y')
% axis equal
% ylim([min(y) max(y)])
% xlim([min(x) max(x)])
% colorbar('ticklabelinterpreter','latex','fontsize',14)
% colormap(cmap)
% set(gca,'ticklabelinterpreter','latex','fontsize',14)
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)
% caxis(c_range_std)