input_data = findfilestrloc(save_dr,'sim',1);
load([save_dr,findfile(save_dr,'classify_',[])]);
load([save_dr,findfile(save_dr,'singlecell_',[])],'tobs','ninds','algout','neighbs','Xscell_obs','Vscell_obs');
tcap=1;
test_tinds_frac= 0.25;
errfun = @(X,Y,V,W) norm(reshape(V(:,:,1:floor(tcap*end))-W(:,:,1:floor(tcap*end)),[],1))/...
    norm(reshape(W(:,:,1:floor(tcap*end)),[],1));
num_gm_tries = 20;
test_tinds = 1:floor(length(tobs)*test_tinds_frac);

if spec<=length(species_inds)
    
    subinds = oppinds(cell2mat(species_inds(1:spec-1)),length(ninds));
    valid_cells=ninds(subinds);
    valpairs = valpairs_cell{spec};
    frelerr = log10(compute_errs(valid_cells,valpairs,{Xscell_obs{1}(:,:,test_tinds)},{Vscell_obs{1}(:,:,test_tinds)},[],errfun));
    
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
    
    xlim([min(egm) max(egm)])

end