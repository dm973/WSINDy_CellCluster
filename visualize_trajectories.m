%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load data
clear all
load('~/Desktop/blah/save_dr/learnedtraj_sim08-Feb-2022_000-111_1_comb4000.mat')

%%% Indicate species and trajectories to plot
species_ind=2;
traj_inds=1:100;

%%% Set line width and line style
lw=3; lwbg=2;
lslearned='-.';
widen_view=0.02;

%%% Set colors
Cbgcells=[0.6784    0.6784    0.6784]+0.15; % gray
CspeciesA_o='black';
CspeciesB_o='black';
CspeciesC_o='black';
CspeciesA=[0.7608 0.2784 1.0000];
CspeciesB=[0 1 1];
CspeciesC=[1 0 0];
CspeciesX=[0.85 0.3275 0.1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot background trajectories
clf;hold on
h1 = plot(squeeze(Xscell_obs{1}(:,1,:))',squeeze(Xscell_obs{1}(:,2,:))');
set(h1, 'color',Cbgcells,'linewidth',lwbg);

%%% plot learned and original trajectories
if species_ind>length(species_models)
    disp(['total number of identified species exceeded'])
else
    if isempty(species_models{species_ind})
        disp(['no model identified for given species'])
    else    
        forceDNA=cellfun(@(x) ~isempty(x), species_models{species_ind}(2:4));
        [f_learned,h_learned,d_learned]=species_models{species_ind}{2:4};
        if isequal(forceDNA,[1 1 1])
            Cspec=CspeciesA_o;
            Clearned=CspeciesA;
        elseif isequal(forceDNA,[1 0 1])
            Cspec=CspeciesB_o;
            Clearned=CspeciesB;
        elseif isequal(forceDNA,[0 1 1])
            Cspec=CspeciesC_o;
            Clearned=CspeciesC;
        else
            Cspec=CspeciesC_o;
            Clearned=CspeciesX;
        end

        X_pred = X_preds{species_ind}{1}(unique(min(traj_inds,size(X_preds{species_ind}{1},1))));
        inds_n = X_preds{species_ind}{2}(unique(min(traj_inds,size(X_preds{species_ind}{1},1))));
        X_p = cell2mat(cellfun(@(x) reshape(squeeze(x{1})',[],1),X_pred,'uni',0)');
        X_p = permute(cat(3,X_p(1:end/2,:),X_p(end/2+1:end,:)),[2 3 1]);
    
        hold on
        for k=1:length(inds_n)
            h2 = plot(squeeze(Xscell_obs{1}(inds_n(k),1,:))',squeeze(Xscell_obs{1}(inds_n(k),2,:))','-');
            h3 = plot(squeeze(X_p(k,1,:))',squeeze(X_p(k,2,:))',lslearned);
            h4 = plot(squeeze(Xscell_obs{1}(inds_n(k),1,1)),squeeze(Xscell_obs{1}(inds_n(k),2,1)),'g.','markersize',30);
            h5 = plot(squeeze(Xscell_obs{1}(inds_n(k),1,end)),squeeze(Xscell_obs{1}(inds_n(k),2,end)),'k.',...
                'markersize',30,'linewidth',2);
            set(h2, 'color',Cspec,'linewidth',lw);
            set(h3, 'color',Clearned,'linewidth',lw);
        end
        legend([h2(1);h3(1);h4(1);h5(1);h1(1)],{'original','learned','start','end','neighbors'},'location','bestoutside','interpreter','latex','fontsize',11)
        xlabel('$x$','interpreter','latex')
        ylabel('$y$','interpreter','latex')
        set(gca,'ticklabelinterpreter','latex','fontsize',14)
    
        axis equal
        xlims = [min(min(X_p(:,1,:)))-widen_view max(max(X_p(:,1,:)))+widen_view];
        ylims = [min(min(X_p(:,2,:)))-widen_view max(max(X_p(:,2,:)))+widen_view];        
        xlim(xlims)
        ylim(ylims)
    end
end