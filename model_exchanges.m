% ninds is a list of particle indices in the data Xscell_obs. 
% ith entry of Mod is an index j such that particle ninds(i)'s model is
% replaced by ninds(j)'s model

function Mod=model_exchanges(ninds,neighbs,ge,algout,alpha,beta,gamma,normGE,normALL,tol1,tol2)
    Mod=(1:length(ninds))';
    if any([alpha>0 beta>0 gamma>0])
        for i=1:length(ninds)
            in=ninds(i);
            inds=cellfun(@(x)find(in==x),neighbs,'uni',0);
            inds_j=cellfun(@(x)~isempty(x),inds);
            inds_j(i)=0;
            inds_j = find(inds_j);
            inds_ji=cell2mat(inds(inds_j));
            
            if normGE==0
                ge_fun = @(L) min(L);
            else
                ge_fun = @(L) norm(L,normGE);
            end
    
            errii=norm([alpha*ge(i,1) beta*norm(algout{i}{5}) gamma*ge_fun(ge(i,2:end))],normALL);
            replacement_stats=[];
            for j=1:length(inds_j)
                errjj=norm([alpha*ge(inds_j(j),1) beta*norm(algout{inds_j(j)}{5}) gamma*ge_fun(ge(inds_j(j),2:end))],normALL);
                errji=norm([alpha*ge(inds_j(j),inds_ji(j)) beta*norm(algout{inds_j(j)}{5}) gamma*ge_fun(ge(inds_j(j),2:end))],normALL);
                errij=inf;%norm([alpha*ge(i,neighbs{i}==ninds(inds_j(j))) beta*norm(algout{i}{5}) gamma*ge_fun(ge(i,2:end))],normALL);
                if all([errii*tol1>errji errij*tol1>errjj max(errji,errjj)<tol2])
                    replacement_stats=[replacement_stats;[inds_j(j) errii errji errjj errij]];
                end
            end
    
            if ~isempty(replacement_stats)
                [~,a]=min(vecnorm(replacement_stats(:,3:4),2,2));
                Mod(i)=replacement_stats(a,1);
                Mod(Mod(1:i-1)==i)=replacement_stats(a,1);
            end
    
        end
    end
end