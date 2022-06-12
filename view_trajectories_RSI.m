%% Load data 
clear all
dr='/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_hetcells/';
% ls(dr)
tag='110';
% load([dr,'data',tag,'.mat'])
load('save_dr/singlecell_0_sim08-Feb-2022_000-111_1_comb4000_clusterloop_11-Jun-2022_learnedX.mat')

% %% Examine individual trajectories
% 
% spec=2;
% 
% X_pred = X_preds{spec}{1};
% inds_n = X_preds{spec}{2};
% X_pred_mat = cell2mat(cellfun(@(x) x{1}, X_pred, 'uni',0));
% V_pred_mat = cell2mat(cellfun(@(x) x{2}, X_pred, 'uni',0));
% 
% tcap=50;
% relerr = cellfun(@(i)errfun(X_pred_mat(i,:,50),Xscell_obs{1}(inds_n(i),:,1:50),V_pred_mat(i,:,1:50),Vscell_obs{1}(inds_n(i),:,1:50)),num2cell(1:length(inds_n))');
% 
% j=1;s=[];clf;
% while isempty(s)
%     plot(squeeze(Xscell_obs{1}(inds_n(j),1,:)),squeeze(Xscell_obs{1}(inds_n(j),2,:)),squeeze(X_pred{j}{1}(1,1,:)),squeeze(X_pred{j}{1}(1,2,:)),'k--','LineWidth',2)
%     hold on; scatter(squeeze(Xscell_obs{1}(inds_n(j),1,tcap)),squeeze(Xscell_obs{1}(inds_n(j),2,tcap)),100,'r','filled'); hold off
%     axis equal
%     s=input('proceed?');
%     title(['j=',num2str(j),'; VE=',num2str(relerr(j))])
%     drawnow
%     j=j+1;
% end
% 
% %% plot background
bgc = 'gray'; %red % 0.1

subspec=4;
subtraj=1:100;
opt=3;

lw=3;
lwbg=2;
lsbg='-.';
lslearned='-.';

xL = 3;
xU = 4.2;
yL = 3;
yU = 4.2;

CspeciesA_o='black';
CspeciesB_o='black';
CspeciesC_o='black';
CspeciesA=[0.7608 0.2784 1.0000];
CspeciesB=[0 1 1];'cyan';
CspeciesC=[1 0 0];'red';

CspeciesX=[0.85 0.3275 0.1];

Clearned={};%[0 0 0];'black';
specs=1:sum(cellfun(@(x) ~isempty(x),species_models));

Cspec={};
for k=1:length(specs)
    forceDNA=cellfun(@(x) ~isempty(x), species_models{k}(2:4));
    if isequal(forceDNA,[1 1 1])
        Cspec{k}=CspeciesA_o;
        Clearned{k}=CspeciesA;
    elseif isequal(forceDNA,[1 0 1])
        Cspec{k}=CspeciesB_o;
        Clearned{k}=CspeciesB;
    elseif isequal(forceDNA,[0 1 1])
        Cspec{k}=CspeciesC_o;
        Clearned{k}=CspeciesC;
    else
        Cspec{k}=CspeciesC_o;
        Clearned{k}=CspeciesX;
    end
end

clf;hold on
if isequal(bgc,'gray')
    Cbackgroundcells=[0.6784    0.6784    0.6784]+0.15; % gray
    h1=plot(squeeze(Xscell_obs{1}(:,1,:))',squeeze(Xscell_obs{1}(:,2,:))');
    set(h1, 'color',Cbackgroundcells,'linewidth',lwbg);
elseif isequal(bgc,'red')
    Cbackgroundcells=[1.0000    0.6784    0.6784]; % faded red
    h1=plot(squeeze(Xscell_obs{1}(:,1,:))',squeeze(Xscell_obs{1}(:,2,:))');
    set(h1, 'color',Cbackgroundcells,'linewidth',lwbg);
elseif and(length(bgc)==3,class(bgc)=='double')
    Cbackgroundcells=bgc;
    h1=plot(squeeze(Xscell_obs{1}(:,1,:))',squeeze(Xscell_obs{1}(:,2,:))');
    set(h1, 'color',Cbackgroundcells,'linewidth',lwbg);
elseif isequal(class(bgc),'double')
    bgcells = {};
    for k=1:length(specs)
        h1 = plot(squeeze(Xscell_obs{1}(species_inds{k},1,:))',squeeze(Xscell_obs{1}(species_inds{k},2,:))');
        set(h1, 'color',min(Cspec{k}+bgc,1),'linewidth',lwbg, 'linestyle',lsbg);
    end
end

for spec=min(subspec,max(specs))

    [f_learned,h_learned,d_learned]=species_models{spec}{2:4};

X_pred = X_preds{spec}{1};
inds_n = X_preds{spec}{2};

if opt==1
    box_inds = cellfun(@(x) and(and(all(x{1}(:,1,:)>xL,3),all(x{1}(:,1,:)<xU,3)),and(all(x{1}(:,2,:)>yL,3),all(x{1}(:,2,:)<yU,3))),X_pred);
    X_pred = X_pred(box_inds);
    inds_n = inds_n(box_inds);
    X_pred_mat = cell2mat(cellfun(@(x) x{1}, X_pred, 'uni',0));
    V_pred_mat = cell2mat(cellfun(@(x) x{2}, X_pred, 'uni',0));
    relerr = cellfun(@(i)errfun(X_pred_mat(i,:,:),Xscell_obs{1}(inds_n(i),:,:),V_pred_mat(i,:,:),Vscell_obs{1}(inds_n(i),:,:)),num2cell(1:length(inds_n))');
    [~,foo] = sort(relerr);
    inds_n = inds_n(unique(foo(min(subtraj,length(foo)))));
    X_pred = X_pred(unique(foo(min(subtraj,length(foo)))));
    xlims = [xL xU];
    ylims = [yL yU];
    X_p = cell2mat(cellfun(@(x) reshape(squeeze(x{1})',[],1),X_pred,'uni',0)');
    X_p = permute(cat(3,X_p(1:end/2,:),X_p(end/2+1:end,:)),[2 3 1]);
    inds_n
elseif opt==2
    box_inds = cellfun(@(x) and(and(all(x{1}(:,1,:)>xL,3),all(x{1}(:,1,:)<xU,3)),and(all(x{1}(:,2,:)>yL,3),all(x{1}(:,2,:)<yU,3))),X_pred);
    X_pred = X_pred(box_inds);
    inds_n = inds_n(box_inds);
    bend_inds = cellfun(@(x) max(range(x{2}(:,1,:),3),range(x{2}(:,2,:),3)),X_pred);
    bend_inds = cellfun(@(x) max(max(diff(x{1},2,3),[],3),[],2),X_pred);
    [~,foo] = sort(bend_inds,'descend');
    inds_n = inds_n(foo(min(subtraj,length(foo))));
    X_pred = X_pred(foo(min(subtraj,length(foo))));
    xlims = [xL xU];
    ylims = [yL yU];
    X_p = cell2mat(cellfun(@(x) reshape(squeeze(x{1})',[],1),X_pred,'uni',0)');
    X_p = permute(cat(3,X_p(1:end/2,:),X_p(end/2+1:end,:)),[2 3 1]);
elseif opt ==3
    X_pred = X_preds{spec}{1}(unique(min(subtraj,size(X_preds{spec}{1},1))));
    inds_n = X_preds{spec}{2}(unique(min(subtraj,size(X_preds{spec}{1},1))));
    X_p = cell2mat(cellfun(@(x) reshape(squeeze(x{1})',[],1),X_pred,'uni',0)');
    X_p = permute(cat(3,X_p(1:end/2,:),X_p(end/2+1:end,:)),[2 3 1]);
    xlims = [min(min(X_p(:,1,:)))-0.02 max(max(X_p(:,1,:)))+0.02];
    ylims = [min(min(X_p(:,2,:)))-0.02 max(max(X_p(:,2,:)))+0.02];
end


hold on
for k=1:length(inds_n) 
    h2 = plot(squeeze(Xscell_obs{1}(inds_n(k),1,:))',squeeze(Xscell_obs{1}(inds_n(k),2,:))','-');
    h3 = plot(squeeze(X_p(k,1,:))',squeeze(X_p(k,2,:))',lslearned);
    h4 = plot(squeeze(Xscell_obs{1}(inds_n(k),1,1)),squeeze(Xscell_obs{1}(inds_n(k),2,1)),'g.','markersize',30);
    h5 = plot(squeeze(Xscell_obs{1}(inds_n(k),1,end)),squeeze(Xscell_obs{1}(inds_n(k),2,end)),'k.',...
        'markersize',30,'linewidth',2);
    set(h2, 'color',Cspec{spec},'linewidth',lw);
%     set(h5, 'color',  [1 0 1]);
    set(h3, 'color',Clearned{spec},'linewidth',lw);
end
legend([h2(1);h3(1);h4(1);h5(1);h1(1)],{'original','learned','start','end','neighbors'},'location','bestoutside','interpreter','latex','fontsize',11)
% legend([h2(1);h3(1);h1(1)],{'original','learned','neighbors'},'location','bestoutside','interpreter','latex','fontsize',11)
% axis equal
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',14)
% set(gcf,'position',[1377         544         544         538])


end    

axis equal

xlim(xlims)
ylim(ylims)

% saveas(gcf,['~/Desktop/traj_',fl(strfind(fl,'_000')+(5:7)),'_',num2str(spec),'.png'])
