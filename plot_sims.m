clear all; close all
expers = {'100','010','001','110','101','011','111'};

figure(1)
for k=1:length(expers)
    exper = expers{k};
    data_dr = ['~/Desktop/data_dr/',exper,'/'];
    save_dr=data_dr;
    input_data = findfilestrloc(save_dr,'sim',1);
    load([data_dr,input_data],'Xscell','Vscell','xmus','xsigs','vmus','vsigs','t','inds_cell_true')
    
    colorz={[0.7608 0.2784 1.0000],[0 1 1], [1 0 0]};
    for tt=[1 2000]
        clf
        for i=1:length(inds_cell_true)
            quiver(Xscell{1}(inds_cell_true{i},1,tt),Xscell{1}(inds_cell_true{i},2,tt),Vscell{1}(inds_cell_true{i},1,tt),Vscell{1}(inds_cell_true{i},2,tt),2,'k','linewidth',1.25); 
            hold on;
            h=scatter(Xscell{1}(inds_cell_true{i},1,tt),Xscell{1}(inds_cell_true{i},2,tt),'filled');
            set(h,'markeredgecolor','k','markerfacecolor',colorz{i})
            hold on
        end
        axis equal
        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        set(gca,'FontSize',14,'TickLabelInterpreter','latex')
        % legend(legs,'interpreter','latex','location','ne','fontsize',30)
        xlim([2 6])
        ylim([2 6])
        saveas(gcf,['~/Desktop/sim_t=',num2str(tt),'_',exper,'.png'])
    end
end