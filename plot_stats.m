%% average speed
clear all;

close all

tic,
expers = {'100','010','001'};
specs = [1 2 3];

numbins = 500;
Rcap = inf;
tpoints = 200;
    
hxxs = zeros(numbins,length(expers));
hvvs = zeros(numbins,length(expers));
hvs = zeros(numbins,length(expers));
polarXs = zeros(200,length(expers));
angMoms = zeros(200,length(expers));

for jj=1:length(expers)
    exper = expers{jj};
    spec = specs(jj);

%     if isequal(exper,'001')
% 
%         dr = '/home/danielmessenger/Dropbox/Boulder/research/data/wound_healing/artificial_cells/multi-species/';
%         subdr = 'multi-species-02-03-22-000-001/';
%         sim = 'sim03-Feb-2022_000_001_2.mat';
%         load([dr,subdr,sim],'Xscell','Vscell','inds_cell_true')

%     else


        data_dr = ['~/Desktop/data_dr/',exper,'/'];
        save_dr = data_dr;
        input_data = findfilestrloc(save_dr,'sim',1);
        load([data_dr,input_data],'Xscell','Vscell','inds_cell_true')
% 
%     end
    
    Xall = Xscell{1}(:,:,1:10:tpoints*10);
    X = Xall(inds_cell_true{spec},:,:);
    Vall = Vscell{1}(:,:,1:10:tpoints*10);
    V = Vall(inds_cell_true{spec},:,:);

    e=linspace(0,3.5,numbins+1);
    xx_grid=(e(1:end-1)+e(2:end))/2;
    hxx = crossparticle_dist(X,Xall,e);
    hxxs(:,jj) = hxx;
  
    e=linspace(0,0.1,numbins+1);
    vv_grid=(e(1:end-1)+e(2:end))/2;
    hvv = crossparticle_dist(V,Vall,e);
    hvvs(:,jj) = hvv; 
    
    X = Vscell{1}(inds_cell_true{spec},:,1:10:end);
    v_grid=vv_grid;
    hv = histcounts(vecnorm(X,2,2),e,'normalization','pdf');
    hvs(:,jj) = hv; 

    polarX = zeros(tpoints,1);
    angMom = zeros(tpoints,1);
    for tt = 1:tpoints
        cm = mean(Xall(:,:,tt));
        Xtt = X(:,:,tt)-cm;    
        Vtt = V(:,:,tt);
        Xmean = mean(Vtt);
        Xabsmean = mean(vecnorm(Vtt,2,2));
        polarX(tt) = norm(Xmean/Xabsmean);
        angMom(tt) = abs(mean(dot(Xtt,[Vtt(:,2) -Vtt(:,1)]))/dot(vecnorm(Xtt,2,2),vecnorm(Vtt,2,2)));
%         plot(Xall(:,1,tt),Xall(:,2,tt),'o');drawnow
    end
    polarXs(:,jj) = polarX;
    angMoms(:,jj) = angMom;

end
toc

%%

legs = {'$\mathbf{X}_{A}$','$\mathbf{X}_{B}$','$\mathbf{X}_{C}$'};
colorz = {[0.7608,0.2784,1.0000],[0,1,1],[1,0,0]};

figure(1)
h=plot(xx_grid,hxxs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$r$','Interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
legend(legs,'interpreter','latex','location','ne','fontsize',30)
xlim([0 3.5])
ylim([0 1])
saveas(gcf,['~/Desktop/prr_homog.png'])

h=plot(vv_grid,hvvs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$s$','Interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
% legend(legs,'interpreter','latex','location','ne','fontsize',30)
xlim([0 0.1])
ylim([0 80])
saveas(gcf,['~/Desktop/pvv_homog.png'])

h=plot(v_grid,hvs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$s$','Interpreter','latex')
% legend('$d\rho_{v}/dr$','interpreter','latex','location','ne','fontsize',30)
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 0.1])
ylim([0 270])
saveas(gcf,['~/Desktop/pv_homog.png'])

h=semilogy(1:size(polarXs,1),polarXs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$t$','Interpreter','latex')
% legend('$P(t)$','interpreter','latex','location','best')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 200])
ylim([0.03 1])
saveas(gcf,['~/Desktop/polar_homog.png'])

h=semilogy(1:size(polarXs,1),angMoms,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$t$','Interpreter','latex')
% legend('$M_{ang}(t)$','interpreter','latex','location','best')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 200])
ylim([0.005 1])
saveas(gcf,['~/Desktop/angMom_homog.png'])

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% average speed
clear all;

close all

tic,
expers = {'110'};
cellspec = {'$A$','$B$','$C$'};
specs_cell = cellfun(@(x) find(str2num(x')),expers,'uni',0);
legs_cell = cellfun(@(x) cellspec(find(str2num(x'))),expers,'uni',0);

numbins = 500;
Rcap = inf;
tpoints = 400;
    
hxxs = [];
hvvs = [];
hvs = [];
polarXs = [];
angMoms = [];

for jj=1:length(expers)
    legs = legs_cell{jj};
    exper = expers{jj};
    specs = specs_cell{jj};
    for kk=1:length(specs)
        spec = specs(kk);

%     if isequal(exper,'001')
% 
%         dr = '/home/danielmessenger/Dropbox/Boulder/research/data/wound_healing/artificial_cells/multi-species/';
%         subdr = 'multi-species-02-03-22-000-001/';
%         sim = 'sim03-Feb-2022_000_001_2.mat';
%         load([dr,subdr,sim],'Xscell','Vscell','inds_cell_true')

%     else


        data_dr = ['~/Desktop/data_dr/',exper,'/'];
        save_dr = data_dr;
        input_data = findfilestrloc(save_dr,'sim',1);
        load([data_dr,input_data],'Xscell','Vscell','inds_cell_true')
% 
%     end
    
    Xall = Xscell{1}(:,:,1:10:tpoints*10);
    X = Xall(inds_cell_true{spec},:,:);
    Vall = Vscell{1}(:,:,1:10:tpoints*10);
    V = Vall(inds_cell_true{spec},:,:);

    e=linspace(0,3.5,numbins+1);
    xx_grid=(e(1:end-1)+e(2:end))/2;
    hxx = crossparticle_dist(X,Xall,e);
    hxxs = [hxxs hxx(:)];
  
    e=linspace(0,0.1,numbins+1);
    vv_grid=(e(1:end-1)+e(2:end))/2;
    hvv = crossparticle_dist(V,Vall,e);
    hvvs = [hvvs hvv(:)]; 
    
    X = Vscell{1}(inds_cell_true{spec},:,1:10:end);
    v_grid=vv_grid;
    hv = histcounts(vecnorm(X,2,2),e,'normalization','pdf');
    hvs = [hvs hv(:)]; 

    polarX = zeros(tpoints,1);
    angMom = zeros(tpoints,1);
    for tt = 1:tpoints
        cm = mean(Xall(:,:,tt));
        Xtt = X(:,:,tt)-cm;    
        Vtt = V(:,:,tt);
        Xmean = mean(Vtt);
        Xabsmean = mean(vecnorm(Vtt,2,2));
        polarX(tt) = norm(Xmean/Xabsmean);
        angMom(tt) = abs(mean(dot(Xtt,[Vtt(:,2) -Vtt(:,1)]))/dot(vecnorm(Xtt,2,2),vecnorm(Vtt,2,2)));
%         plot(Xall(:,1,tt),Xall(:,2,tt),'o');drawnow
    end
    polarXs = [polarXs polarX(:)];
    angMoms = [angMoms angMom(:)];
    end
end
toc

%%

colorz = {[0.7608,0.2784,1.0000],[0,1,1],[1,0,0]};

figure(1)
h=plot(xx_grid,hxxs,'LineWidth',4);
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$r$','Interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
legend(legs,'interpreter','latex','location','ne','fontsize',30)
xlim([0 3.5])
ylim([0 0.8])
saveas(gcf,['~/Desktop/prr_',exper,'.png'])

%%

h=plot(vv_grid,hvvs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$s$','Interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
% legend(legs,'interpreter','latex','location','ne','fontsize',30)
xlim([0 0.1])
ylim([0 55])
saveas(gcf,['~/Desktop/pvv_',exper,'.png'])

%%

h=plot(v_grid,hvs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$s$','Interpreter','latex')
% legend('$d\rho_{v}/dr$','interpreter','latex','location','ne','fontsize',30)
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 0.1])
ylim([0 150])
saveas(gcf,['~/Desktop/pv_',exper,'.png'])

%%

h= semilogy(1:size(polarXs,1),polarXs,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$t$','Interpreter','latex')
% legend('$P(t)$','interpreter','latex','location','best')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 400])
ylim([0.005 1])
saveas(gcf,['~/Desktop/polar_',exper,'.png'])

%%

h=semilogy(1:size(polarXs,1),angMoms,'LineWidth',4)
for j=1:length(specs)
    set(h(j),'color',colorz{j})
end
xlabel('$t$','Interpreter','latex')
% legend('$M_{ang}(t)$','interpreter','latex','location','best')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 400])
ylim([0.00001 1])
saveas(gcf,['~/Desktop/angMom_',exper,'.png'])

