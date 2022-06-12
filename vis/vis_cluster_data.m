clear all; addpath(genpath('~/repos/WSINDy_SDE_kitchen'));addpath(genpath('~/Dropbox/Boulder/research/data/wound_healing/'))
dr='~/Dropbox/Boulder/research/data/wound_healing/single-cell/';
S=dir(dr);
fl_names = {S(cellfun(@(x) contains(x,'clusterloop') , {S.name})).name}';
fl_singlecell = cellfun(@(x) [x(1:strfind(x,'_clusterloop')-1),'.mat'], fl_names,'uni',0);

[mat2cell((1:length(fl_names))',ones(length(fl_names),1),1),fl_names]

%% load data
clc

%[2 1 7 4 23 5 14] %%% [100 010 001 101 011 110 111]
i=17;
fl=fl_singlecell{i}
load(fl_names{i},'species_models','species_inds','cluster_loop_txt','errfun','valpairs_cell')
load(fl,'input_data','ninds','Xscell_obs','Vscell_obs','tobs','simdat','tol_dx','tol_fx','J_fv','J_fx','J_hv','J_hx','J_dv','J_dx')
disp(species_inds)
[N,d,M]=size(Xscell_obs{1});

try
    load(input_data,'inds_cell_true','inds_cell','f_xv_true_cell','h_xv_true_cell','d_xv_true_cell','f_xv_true','h_xv_true','d_xv_true')
    if ~exist('inds_cell_true','var')
        try
            inds_cell_true=inds_cell;
        catch
            inds_cell_true={1:N};
            f_xv_true_cell={f_xv_true};
            h_xv_true_cell={h_xv_true};
            d_xv_true_cell={d_xv_true};
        end
    end
catch
    f_xv_true = [];
    h_xv_true = [];
    d_xv_true = [];
    inds_cell_true = {};
end

foo1=strfind(cluster_loop_txt,'tcap=');
foo2=strfind(cluster_loop_txt,'Xtest=X')-1;
eval(cluster_loop_txt(foo1:foo2))

%%
% f,h,d
disp(fl)
numbins = 2000; maxents= 10^7; R = inf;
[h,e] = particle_dist(Xscell_obs, numbins,maxents,R);
f_true = f_xv_true_cell{1};
h_true = h_xv_true_cell{1};
d_true = d_xv_true_cell{1};

for spec=1:length(find(cellfun(@(x) ~isempty(x),species_models)))
    Wdirmodes = dirmodes({species_models{spec}{1}},J_fx,J_fv,J_hx,J_hv,J_dx,J_dv);
    spec_frac = cellfun(@(x) length(intersect(x,species_inds{spec}))/length(intersect(x,ninds)),inds_cell_true);
    spec_frac(isnan(spec_frac))=0;
    [cs,spec_majority] = max(spec_frac);
    disp(spec_frac)
    spec_frac = cellfun(@(x) length(intersect(x,species_inds{spec})),inds_cell_true);
    spec_frac(isnan(spec_frac))=0;
    disp(spec_frac)
    try
        f_learned = species_models{spec}{2};
    catch
        f_learned = [];
    end
    if ~isempty(f_learned)

        figure(1);clf
        pimults = [0 pi/4 pi/2];
        F_true = vis_polar_line(f_true,e,pimults,'$f_{a-r}$','$f_{a-r}^\star$','$r$','-');
        F_learned = vis_polar_line(f_learned,e,pimults,'$f_{a-r}$','$\widehat{f}_{a-r}$','$r$','--');
        area(e,h*max(abs(F_true(:)))/max(h)/4,'faceAlpha',0.5,'DisplayName','$\rho_{rr}$')
        h1=legend('show'); set(h1,'Interpreter','latex','box','off','location','northeast')
        ylim([-1 18])

        err = (sum(sum(h'.*(F_learned-F_true).^2))/sum(sum(h'.*F_true.^2)))^(1/2);
        disp(['f ',num2str(spec), 'err = ',num2str(err)])
        if and(err<0.1,err>0.01)
            axes('Position',[.2 .5 .3 .36]);box on
            F_true = vis_polar_line(f_xv_true_cell{spec_majority},e,pimults,[],[],'$r$','-');
            F_learned = vis_polar_line(f_learned,e,pimults,[],[],[],'--');
            area(e,h*max(abs(F_true(:)))/max(h)/4,'faceAlpha',0.5,'DisplayName','$\rho_{rr}$')
            xlim([0 0.1])
            ylim([-1 18])
        end
        drawnow
%         saveas(gcf,['~/Desktop/f_',fl(strfind(fl,'_000')+(5:7)),'_',num2str(spec),'.png'])
    end
    
    %%% h

     try
        h_learned = species_models{spec}{3};
     catch
        h_learned = [];
     end
     if ~isempty(h_learned)

        figure(2);clf
        pimults = [0 pi/4 pi/2];
        F_true = vis_polar_line(h_true,e,pimults,'$f_{align}$','$f_{align}^\star$','$r$','-');
        F_learned = vis_polar_line(h_learned,e,pimults,'$f_{align}$','$\widehat{f}_{align}$','$r$','--');
        area(e,h*max(abs(F_true(:)))/max(h)/4,'faceAlpha',0.5,'DisplayName','$\rho_{rr}$')
        h1=legend('show'); set(h1,'Interpreter','latex','box','off','location','southeast')
        ylim([-16 8])
        
        err = (sum(sum(h'.*(F_learned-F_true).^2))/sum(sum(h'.*F_true.^2)))^(1/2);
        disp(['h ',num2str(spec), 'err = ',num2str(err)])
        if and(err<0.1,err>0.01)
            axes('Position',[.22 .2 .3 .36]);box on
            F_true = vis_polar_line(h_xv_true_cell{spec_majority},e,pimults,[],[],'$r$','-');
            F_learned = vis_polar_line(h_learned,e,pimults,[],[],[],'--');
            area(e,h*max(abs(F_true(:)))/max(h)/4,'faceAlpha',0.5,'DisplayName','$\rho_{rr}$')
            xlim([0 0.2])
            ylim([-15 0])
        end
%         saveas(gcf,['~/Desktop/h_',fl(strfind(fl,'_000')+(5:7)),'_',num2str(spec),'.png'])
    end
     try
        d_learned = species_models{spec}{4};
     catch
        d_learned = [];
     end
    if ~isempty(d_learned)

        figure(3);clf
        pimults = [0];
        numbins = 2000; maxents= 10^7; R = inf;
        [hv,ev] = histcounts(reshape(vecnorm(Vscell_obs{1}(species_inds{spec},:,:),2,2),[],1), numbins,'normalization','probability');
        ev = ev(cumsum(hv)/sum(hv)<0.99);
        hv = hv(cumsum(hv)/sum(hv)<0.99);
        F_true = vis_polar_line(d_xv_true_cell{spec_majority},ev,pimults,[],'$f^\star_{drag}$',[],'-');
        F_learned = vis_polar_line(d_learned,ev,pimults,'$f_{drag}$','$\widehat{f}_{drag}$','$|v|$','--');
        area(ev,hv*max(abs(F_true(:)))/max(hv)/4,'faceAlpha',0.5,'DisplayName','$\rho_{|v|}$')
        h1=legend('show'); set(h1,'Interpreter','latex','box','off','location','southwest')
        
        err = (sum(sum(hv'.*(F_learned-F_true).^2))/sum(sum(hv'.*F_true.^2)))^(1/2);
        disp(['d ',num2str(spec), 'err = ',num2str(err)])
        if and(err<0,err>0.01)
            axes('Position',[.2 .1 .3 .36]);box on
            F_true = vis_polar_line(d_true,ev,pimults,[],[],'$r$','-');
            F_learned = vis_polar_line(d_learned,ev,pimults,[],[],[],'--');
            area(ev,hv*max(abs(F_true(:)))/max(hv)/4,'faceAlpha',0.5,'DisplayName','$\rho_{rr}$')
            xlim([0 0.3])
            ylim([-1 3])
        end
%         saveas(gcf,['~/Desktop/d_',fl(strfind(fl,'_000')+(5:7)),'_',num2str(spec),'.png'])
    end
    
    %%%%% view gaussian mixtures

    valpairs=valpairs_cell{spec};
    valid_cells =cellfun(@(x)x{3},valpairs);

    %%% compute log val. errs
    frelerr = log10(compute_errs(valid_cells,valpairs,{Xscell_obs{1}(:,:,test_tinds)},...
        {Vscell_obs{1}(:,:,test_tinds)},[],errfun));

    gm1=fitgmdist(frelerr(:),1);
    idk1 = cluster(gm1,frelerr(:));

    if length(species_inds{spec})>2
        idks=repmat(idk1*0,1,num_gm_tries);
        bic = zeros(1,num_gm_tries);
        for j=1:num_gm_tries
            gm=fitgmdist(frelerr(:),2,'RegularizationValue',10^-6);
            idk = cluster(gm,frelerr(:));
            bic(j)=gm.BIC;
            [~,ii]=min(gm.mu);
            idks(:,j)= idk==ii;
        end
    end
    figure(4);clf
    if and(mean(bic) < gm1.BIC, all([sum(frelerr>log10(0.05))/length(frelerr)>0.01 length(frelerr)>1]))
        x = linspace(min(gm.mu)-sqrt(max(gm.Sigma))*3,max(gm.mu)+sqrt(max(gm.Sigma))*3,1000);
        disp('multi-species')
        clustbar=mean([max(frelerr(ismember(valid_cells,species_inds{spec}))) min(frelerr(~ismember(valid_cells,species_inds{spec})))]);
        y = pdf(gm,x');
        plot([clustbar clustbar],[0 1.8],'--', 'markersize',10,'linewidth',3,'DisplayName','Cluster division');
        disp(['err = ', num2str(mean(10.^frelerr(frelerr<clustbar)))])
    else
        x = linspace(min(gm1.mu)-sqrt(max(gm1.Sigma))*3,max(gm1.mu)+sqrt(max(gm1.Sigma))*3,1000);
        disp('mono-species')
        clusterbar=[];
        y = pdf(gm1,x');
        disp(['err = ', num2str(mean(10.^frelerr))])
    end
    hold on
    [hgm,egm]=histcounts(frelerr,50,'normalization','pdf');
    bar(egm(1:end-1),hgm,'DisplayName','$VE$ histogram');
    plot(x,y,'--', 'markersize',10,'linewidth',3,'DisplayName','GM fit')
    h1=legend('show');
    set(h1,'interpreter','latex','fontsize',14);
    xlabel('$\log_{10}(VE)$','interpreter','latex')
    set(gca,'ticklabelinterpreter','latex','fontsize',14)
%     legend({'Cluster division','GM fit','$VE$ histogram'},'interpreter','latex','fontsize',14)

    xlim([min(egm) max(egm)])

    saveas(gcf,['~/Desktop/gm_',fl(strfind(fl,'_000')+(5:7)),'_',num2str(spec),'.png'])
end
disp('-----------------------------------------------------')
% end
% 
%% view true forces

figure(1);clf
n=500;rr=0.2;cap=inf;levs=200;
cmap=jet(100);
F_true = vis_polar_fcn(f_xv_true_cell{1},n,rr,cap,levs,cmap);
legend({'$f^\star_{a-r}$'},'interpreter','latex','fontsize',12)
caxis([-15 15])
saveas(gcf,['~/Desktop/ftrue.png'])

figure(2);clf
F_true = vis_polar_fcn(h_xv_true_cell{1},n,rr,cap,levs,cmap);
legend({'$f^\star_{align}$'},'interpreter','latex','fontsize',12)
caxis([-15 15])
saveas(gcf,['~/Desktop/htrue.png'])

% figure(3);clf
% pimults = [0 pi/4 pi/2 3*pi/4 pi];
% x_grid=linspace(0,tol_dx,500);
% D_true = vis_polar_line(d_xv_true_cell{1},x_grid,pimults,'$f_{drag}$','$f^\star_{drag}$');
% saveas(gcf,['~/Desktop/dtrue_art.png'])


%% View trajectories

for spec=[1]%1:length(species_models);
% spec=1:4;
tic,
f_learned = species_models{spec}{2};
h_learned = species_models{spec}{3};
d_learned = species_models{spec}{4};

[~,a]=sort(cellfun(@(x) x{2}(1),simdat(species_inds{spec})));
a=a(2);
M=dist(squeeze(Xscell_obs{1}(ninds(species_inds{spec}),:,1))');
[b,a]=sort(M(a,:));cell_inds=a(1:min(8,end));

num_traj=length(cell_inds);
X_pred = cell(num_traj,1);zeros(num_traj,d,length(tobs));
V_pred = cell(num_traj,1);zeros(num_traj,d,length(tobs));
inds_n = ninds(species_inds{spec}(cell_inds));

parfor i=1:num_traj
    RHS = @(Z,t) rhs_lean_pred(f_learned, h_learned, d_learned, Z, d, t, opts, 1);
    inds_np = oppinds(inds_n(i),N);
    Xtemp = Xscell_obs{1}(inds_n(i),:,1);
    Vtemp = mean(Vscell_obs{1}(inds_n(i),:,1:avg_v0),3);
    Ztemp = reshape([Xtemp Vtemp],[],1);
    [~,Z,B] = simnu_split(RHS,tobs,Ztemp,subdt,nux,nuv,verbose,Xscell_obs{1}(inds_np,:,:),Vscell_obs{1}(inds_np,:,:));
    X_pred_temp = zeros(1,d,length(tobs));
    V_pred_temp = zeros(1,d,length(tobs));
    for m=1:length(tobs)
        X_pred_temp(:,:,m) = reshape(Z(m,1:end/2)',1,d);
        V_pred_temp(:,:,m) = reshape(Z(m,end/2+1:end)',1,d);
    end
    X_pred{i}={X_pred_temp,V_pred_temp};
end
toc

Cbackgroundcells=[0.6784    0.6784    0.6784];
CspeciesA=[0.7608    0.2784    1.0000];
CspeciesB='cyan';
CspeciesC='red';
Clearned='black';

forceDNA=cellfun(@(x) ~isempty(x), species_models{spec}(2:4));
if isequal(forceDNA,[1 1 1])
    Cspec=CspeciesA;
elseif isequal(forceDNA,[1 0 1])
    Cspec=CspeciesB;
elseif isequal(forceDNA,[0 1 1])
    Cspec=CspeciesC;
end

X_p = cell2mat(cellfun(@(x) reshape(squeeze(x{1})',[],1),X_pred,'uni',0)');
X_p = permute(cat(3,X_p(1:end/2,:),X_p(end/2+1:end,:)),[2 3 1]);
Y_p = cell2mat(cellfun(@(x) reshape(squeeze(x{1}{1}{1})',[],1),simdat( ismember(ninds,inds_n)),'uni',0)');
Y_p = permute(cat(3,Y_p(1:end/2,:),Y_p(end/2+1:end,:)),[2 3 1]);

% relerr = errfun(X_pred,Xscell_obs{1}(inds_n,:,:),V_pred,Vscell_obs{1}(inds_n,:,:))
clf
h1=plot(squeeze(Xscell_obs{1}(:,1,:))',squeeze(Xscell_obs{1}(:,2,:))','-','markersize',10,'linewidth',0.25);
hold on
h2 = plot(squeeze(Xscell_obs{1}(inds_n,1,:))',squeeze(Xscell_obs{1}(inds_n,2,:))','-');
h3 = plot(squeeze(X_p(:,1,:))',squeeze(X_p(:,2,:))','--');
h4 = plot(squeeze(Xscell_obs{1}(inds_n,1,1)),squeeze(Xscell_obs{1}(inds_n,2,1)),'gx',squeeze(Xscell_obs{1}(inds_n,1,end)),squeeze(Xscell_obs{1}(inds_n,2,end)),'kx',...
    'markersize',10,'linewidth',2);
set(h1, 'color',Cbackgroundcells,'markersize',10,'linewidth',2);
set(h2, 'color',Cspec,'markersize',10,'linewidth',2);
set(h3, 'color',Clearned,'markersize',10,'linewidth',2);
%     squeeze(Y_p(:,1,:))',squeeze(Y_p(:,2,:))','k--.',
legend([h2(1);h3(1);h4;h1(1)],{'original','learned','start','end','neighbors'},'location','bestoutside','interpreter','latex','fontsize',11)
axis equal

xlims = [min(min(X_p(:,1,:)))-0.02 max(max(X_p(:,1,:)))+0.02];
ylims = [min(min(X_p(:,2,:)))-0.02 max(max(X_p(:,2,:)))+0.02];

xlim(xlims)
ylim(ylims)

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',14)
set(gcf,'position',[1377         544         544         538])

% saveas(gcf,['~/Desktop/traj_',fl(strfind(fl,'_000')+(5:7)),'_',num2str(spec),'.png'])
end
%%

dirmodes(species_models{4}(1),J_fx,J_fv,J_hx,J_hv,J_dx,J_dv)

%%
i2=intersect(find(idk==2),1:501);
i=3;
inds_n=valid_cells(i2(i));
X_pred = valpairs{i2(i)}{1};
plot(squeeze(Xscell_obs{1}(inds_n,1,test_tinds)),squeeze(Xscell_obs{1}(inds_n,2,test_tinds)),'r-',squeeze(X_pred(1,1,:)),squeeze(X_pred(1,2,:)),'b--',...
    squeeze(Xscell_obs{1}(inds_n,1,1)),squeeze(Xscell_obs{1}(inds_n,2,1)),'gx',squeeze(Xscell_obs{1}(inds_n,1,test_tinds(end))),squeeze(Xscell_obs{1}(inds_n,2,test_tinds(end))),'kx','markersize',10,'linewidth',2)


function F=vis_polar_fcn(f,n,rr,cap,levs,cmap)

    x = linspace(-rr,rr,n);
    y = linspace(-rr,rr,n);
    [xx,yy] = meshgrid(x,y);
    dX = [xx(:) yy(:)];
    dX = reshape([0 0],1,1,2) - reshape(dX,1,n^2,2);
    rr = reshape(squeeze(vecnorm(dX,2,3)),n,[]);
    v = [0 1];
    dV = repmat(reshape(v,1,1,2),1,n^2,1);
    th = dot(-dX,dV,3);
    th = reshape(atan2(abs(dot(-flip(cat(3,-dX(:,:,1),dX(:,:,2)),3),dV,3)),th),n,[]);

    if ~isempty(f)
        F = min(max(f(rr,th),-cap),cap);
    else
        F = 0*rr;
    end
%     contourf(xx,yy,F,levs,'edgecolor','none')
    surf(xx,yy,F,'edgecolor','none')
    view([0 90])
    axis equal
    xlabel('$x$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    ylim([min(y) max(y)])
    xlim([min(x) max(x)])
    colorbar('ticklabelinterpreter','latex','fontsize',14)
    colormap(cmap)
    set(gca,'ticklabelinterpreter','latex','fontsize',14)
    

end

function F=vis_polar_line(f,x_grid,v_snaps,ylab,leglab,xlab,stl)
    vs = length(v_snaps);
    pimult=v_snaps/pi;
    F = zeros(length(x_grid),vs);
    for i=1:vs
        F(:,i)=f(x_grid,v_snaps(i))';
        if pimult(i)==0
            plot(x_grid,F(:,i),stl,'linewidth',2, 'DisplayName',[leglab, ', $\theta=0$']); hold on
        elseif pimult(i)==1
            plot(x_grid,F(:,i),stl,'linewidth',2, 'DisplayName',[leglab, ', $\theta=\pi$']); hold on
        else
            plot(x_grid,F(:,i),stl,'linewidth',2, 'DisplayName',[leglab, ', $\theta=(',num2str(pimult(i)),')\pi$']); hold on
        end
    end
    xlabel(xlab,'interpreter','latex')
    ylabel(ylab,'interpreter','latex')
    xlim([min(x_grid) max(x_grid)])
    set(gca,'ticklabelinterpreter','latex','fontsize',14)
    if isempty(leglab)
        legend off
    end
end

function G = dirmodes(Ws,J_fx,J_fv,J_hx,J_hv,J_dx,J_dv)
    angleinds_f=cell(J_fv,1);
    for i=1:J_fv
        angleinds_f{i}=i:J_fv:J_fx*J_fv;
    end
    angleinds_h=cell(J_hv,1);
    for i=1:J_hv
        angleinds_h{i}=(i:J_hv:J_hx*J_hv)+J_fx*J_fv;
    end
    angleinds_d=cell(J_dv,1);
    for i=1:J_dv
        angleinds_d{i}=(i:J_dv:J_dx*J_dv)+J_fx*J_fv+J_hx*J_hv;
    end
    G=cell2mat(cellfun(@(y)double([cellfun(@(x)any(y(x)~=0),angleinds_f);cellfun(@(x)any(y(x)~=0),angleinds_h);cellfun(@(x)any(y(x)~=0),angleinds_d)]),Ws,'uni',0));

    num_angle=J_fv+J_hv+J_dv;
    inds_all=de2bi(1:2^num_angle-1);
    ncat = 2^num_angle-1;
% 
%     inds_keep=[];ge_score=[];models={};inds_pat={};ninds_rep=[];
%     for i=1:ncat
%         inds=find(all(G==inds_all(i,:)',1));
end