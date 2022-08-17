%% average speed
% clear all;

close all

tic,
expers = {'100','110','101','111'};
legs = {'$\mathbf{X}_{A}$','$\mathbf{X}_{A,B}$(long)','$\mathbf{X}_{A,C}$','$\mathbf{X}_{A,B,C}$(long)'};
spec=1;

expers = {'010','110','011','111'};
legs = {'$\mathbf{X}_{B}$','$\mathbf{X}_{A,B}$(long)','$\mathbf{X}_{B,C}$','$\mathbf{X}_{A,B,C}$(long)'};
spec=2;

expers = {'001','101','011','111'};
legs = {'$\mathbf{X}_{C}$','$\mathbf{X}_{A,C}$','$\mathbf{X}_{B,C}$','$\mathbf{X}_{A,B,C}$(long)'};
spec=3;


numbins=250;
Rcap = inf;
    
hxxs = zeros(numbins,length(expers));
hvvs = zeros(numbins,length(expers));
hvs = zeros(numbins,length(expers));

for jj=1:length(expers)
    exper = expers{jj};
    data_dr = ['~/Desktop/data_dr/',exper,'/'];
    save_dr = data_dr;
    input_data = findfilestrloc(save_dr,'sim',1);
    load([data_dr,input_data],'Xscell','Vscell','inds_cell_true')
    
    X = Xscell{1}(:,:,1:10:end);
    e=linspace(0,3.5,numbins+1);
    xx_grid=(e(1:end-1)+e(2:end))/2;
    hxx = crossparticle_dist(X(inds_cell_true{spec},:,:),X,e);
    hxxs(:,jj) = hxx;
    
    X = Vscell{1}(:,:,1:10:end);
    e=linspace(0,0.1,numbins+1);
    vv_grid=(e(1:end-1)+e(2:end))/2;
    hvv = crossparticle_dist(X(inds_cell_true{spec},:,:),X,e);
    hvvs(:,jj) = hvv; 
    
    X = Vscell{1}(inds_cell_true{spec},:,1:10:end);
    v_grid=vv_grid;
    hv = histcounts(vecnorm(X,2,2),e,'normalization','pdf');
    hvs(:,jj) = hv; 
end
toc

figure(1)
plot(xx_grid,hxxs,'LineWidth',4)
xlabel('$r$','Interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
% legend(legs,'interpreter','latex','location','ne','fontsize',30)
xlim([0 3.5])
ylim([0 1])
saveas(gcf,['~/Desktop/prr_',num2str(spec),'.png'])

figure(2)
plot(vv_grid,hvvs,'LineWidth',4)
xlabel('$s$','Interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
legend(legs,'interpreter','latex','location','ne','fontsize',30)
xlim([0 0.1])
ylim([0 170])
saveas(gcf,['~/Desktop/pvv_',num2str(spec),'.png'])

figure(3)
plot(v_grid,hvs,'LineWidth',4)
xlabel('$s$','Interpreter','latex')
% legend(legs,'interpreter','latex')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 0.1])
ylim([0 250])
saveas(gcf,['~/Desktop/pv_',num2str(spec),'.png'])

%% average speed

clear all; close all;

expers = {'100','110','101','111'};
legs = {'$\mathbf{X}_{A}$','$\mathbf{X}_{A,B}$(long)','$\mathbf{X}_{A,C}$','$\mathbf{X}_{A,B,C}$(long)'};
spec=1;

expers = {'010','110','011','111'};
legs = {'$\mathbf{X}_{B}$','$\mathbf{X}_{A,B}$(long)','$\mathbf{X}_{B,C}$','$\mathbf{X}_{A,B,C}$(long)'};
spec=2;

expers = {'001','101','011','111'};
legs = {'$\mathbf{X}_{C}$','$\mathbf{X}_{A,C}$','$\mathbf{X}_{B,C}$','$\mathbf{X}_{A,B,C}$(long)'};
spec=3;

polarXs = [];
angMoms = [];
for exper = expers
    data_dr = ['~/Desktop/data_dr/',exper{1},'/'];
    save_dr = data_dr;
    input_data = findfilestrloc(save_dr,'sim',1);
    load([data_dr,input_data],'Xscell','Vscell','inds_cell_true')
    
    Xall = Xscell{1}(:,:,1:10:2000);
    X = Xscell{1}(inds_cell_true{spec},:,1:10:2000);
    V = Vscell{1}(inds_cell_true{spec},:,1:10:2000);
    polarX = zeros(size(X,3),1);
    angMom = zeros(size(X,3),1);
    for tt = 1:size(V,3)
        cm = mean(Xall(:,:,tt));
        Xtt = X(:,:,tt);    
        Vtt = V(:,:,tt);
        Xmean = mean(Vtt);
        Xabsmean = mean(vecnorm(Vtt,2,2));
        polarX(tt) = norm(Xmean/Xabsmean);
        angMom(tt) = abs(mean(dot(Xtt-cm,[Vtt(:,2) -Vtt(:,1)]))/dot(vecnorm(Xtt-cm,2,2),vecnorm(Vtt,2,2)));
        plot(Xall(:,1,tt),Xall(:,2,tt),'o');drawnow
    end
    polarXs = [polarXs polarX(:)];
    angMoms = [angMoms angMom(:)];
end    
    
figure(1)
semilogy(1:size(polarXs,1),polarXs,'LineWidth',4)
xlabel('$t$','Interpreter','latex')
% legend(legs,'interpreter','latex','location','best')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 200])
ylim([0.005 1])
saveas(gcf,['~/Desktop/polar_',num2str(spec),'.png'])

figure(1)
semilogy(1:size(polarXs,1),angMoms,'LineWidth',4)
xlabel('$t$','Interpreter','latex')
% legend(legs,'interpreter','latex','location','best')
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0 200])
ylim([0.000001 1])
saveas(gcf,['~/Desktop/angMom_',num2str(spec),'.png'])
