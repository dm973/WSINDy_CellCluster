% clear all
% 
% dr='~/Desktop/data/wound_healing/single-cell/';
% fl='singlecell_artificial_1-21-22_N1000';
% singlecellfile=[dr,fl,'.mat'];
% topmodval=[dr,fl,'_topmodelval.mat'];
% 
% load(topmodval,'valid_cells','relerr');
% load(singlecellfile,'ninds','algout','simdat');
% 
% errtol=0.1;
% outliertol=0.01;
% inds=valid_cells(relerr>errtol);
% topcluster = find(~ismember(ninds,valid_cells));find(ismember(ninds,valid_cells(relerr<=errtol)))];
% alpha=1/sqrt(3);beta=1/sqrt(3);gamma=1/sqrt(3);pp=0;tol=0.5;
% 
% if length(inds)<outliertol*length(relerr)
%     disp(['unimodal'])
% else
%     disp(['multimodal'])
%     subinds=find(ismember(ninds,inds));
%     [Mod,inds_pat]=clustermodels(singlecellfile,alpha,beta,gamma,pp,tol,subinds);
%     if sum(cellfun(@(x)length(x{1}),inds_pat))<outliertol*length(relerr)
%         disp(['second mode unclassifiable'])
%     end
% end

function [Wsmat,Mod,inds_keep,ge_score,models,inds_pat,ninds_rep]=clustermodels(alpha,beta,gamma,normGE,normALL,tol1,tol2,ninds,simdat,neighbs,algout,J_fv,J_fx,J_hv,J_hx,J_dv,J_dx,ge)

    Mod=model_exchanges(ninds,neighbs,ge,algout,alpha,beta,gamma,normGE,normALL,tol1,tol2);
    algoutMod=algout(Mod);
    simdat=simdat(Mod);
  
    Ws ={};Wsmat=[];res=[];
    for j=1:length(ninds)
        Ws{j}=algoutMod{j}{4};
        Wsmat=[Wsmat Ws{j}]; 
        res(j,1)=norm(algoutMod{j}{5});
    end
    
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

    inds_keep=[];ge_score=[];models={};inds_pat={};ninds_rep=[];
    for i=1:ncat
        inds=find(all(G==inds_all(i,:)',1));
        if ~isempty(inds)
            if normGE==0 
                pscore=vecnorm([alpha*ge(inds,1) beta*res(inds) gamma*min(ge(inds,2:end),[],2)],normALL,2);
            else
                pscore=vecnorm([alpha*ge(inds,1) beta*res(inds) gamma*vecnorm(ge(inds,2:end),normGE,2)],normALL,2);
            end
            if min(pscore)<tol2
                [~,a]=max(pscore);
                models=[models,algoutMod(inds(a))];
                ge_score=[ge_score [ge(inds(a),1);median(ge(inds(a),2:end),2);res(inds(a));length(inds)]];
                inds_keep=[inds_keep inds_all(i,:)'];
                ninds_rep=[ninds_rep [inds(a);ninds(Mod(inds(a)))]];
                inds_pat=[inds_pat,{{inds(pscore<tol2),inds_all(i,:)',inds}}];
            end
        end
    end
    if ~isempty(ge_score)
        [~,a]=sort(ge_score(end,:),'descend');        
        ge_score=ge_score(:,a);
        inds_keep=inds_keep(:,a);
        models=models(a);
        ninds_rep=ninds_rep(:,a);
        inds_pat=inds_pat(a);
    else
        disp('no model found')
    end
    
end