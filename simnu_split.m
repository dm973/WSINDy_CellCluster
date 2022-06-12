% function [t,Z,B] = simnu_split(RHS,t,Ztemp,m,nux,nuv,verbose,X,V)
%     M = length(t);
%     dt = t(2)-t(1);
%     Ntot = length(Ztemp);
%     B = zeros(Ntot,M);
%     Z = zeros(M,Ntot);
%     Z(1,:) = Ztemp;
%     for mm=2:M
%         t_temp = t(mm-1)+dt/m:dt/m:t(mm);
% %         [fcell,vcell]=hermite_cell(X{mm-1},X{mm},V{mm-1},V{mm},dt);
% %         [fcell,vcell]=hermite_cell(squeeze(X(:,:,mm-1)),squeeze(X(:,:,mm)),squeeze(V(:,:,mm-1)),squeeze(V(:,:,mm)),dt);
%         for j=1:length(t_temp)
%            Xtemp = j/length(t_temp)*[X(:,:,mm-1) V(:,:,mm-1)]+(1-j/length(t_temp))*[X(:,:,mm) V(:,:,mm)];
% %            Xtemp = j/length(t_temp)*[X{mm-1} V{mm-1}]+(1-j/length(t_temp))*[X{mm} V{mm}];
% %             Xtemp = [cell2mat(cellfun(@(x)x(t_temp(j)),fcell(:,1),'UniformOutput',false)),...
% %                 cell2mat(cellfun(@(x)x(t_temp(j)),vcell(:,1),'UniformOutput',false))];
%             Ztemp = reshape([reshape(Ztemp,[],4);Xtemp],[],1);            
%             Ztemp = Ztemp+dt/m*RHS(Ztemp,t_temp(j));
%             Ztemp = reshape(Ztemp,[],4);
%             Ztemp = reshape(Ztemp(1:Ntot/4,:),[],1);
%             if or(nux>0,nuv>0)
%                 Btemp = [sqrt(2*nux*(dt/m))*randn(Ntot/2,1);sqrt(2*nuv*(dt/m))*randn(Ntot/2,1)];
%                 Ztemp(1:end/2) = Ztemp(1:end/2) + Btemp(1:end/2);
%                 Ztemp(end/2+1:end) = Ztemp(end/2+1:end) + Btemp(end/2+1:end);
%                 B(:,mm) = B(:,mm)+Btemp;
%             end
%         end
%         Z(mm,:) = Ztemp;
%         if verbose
%             disp(mm)
%         end
%     end
% end    

function [t,Z,B] = simnu_split(RHS,t,Ztemp,m,nux,nuv,verbose,X,V)
    M = length(t);
    Ntot = length(Ztemp);
    B = zeros(Ntot,M);
    Z = zeros(M,Ntot);
    Z(1,:) = Ztemp;
    
    tq=linspace(t(1),t(2),m+1);
    for mm=1:M-1
        ttemp=linspace(t(mm),t(mm+1),m+1);
        tq=[tq,ttemp(2:end)];
    end
    
    Xinterp=zeros(size(X,1),size(X,2),length(tq));
    Vinterp=zeros(size(V,1),size(V,2),length(tq));
    Xinterp(:,1,:)=pchip(t,X(:,1,:),tq);
    Xinterp(:,2,:)=pchip(t,X(:,2,:),tq);
    Vinterp(:,1,:)=pchip(t,V(:,1,:),tq);
    Vinterp(:,2,:)=pchip(t,V(:,2,:),tq);

    for i=2:length(tq)
        dtq=(tq(i)-tq(i-1));
        Xtemp = [Xinterp(:,:,i) Vinterp(:,:,i)];
        Ztemp = reshape([reshape(Ztemp,[],4);Xtemp],[],1);            
        Ztemp = Ztemp+dtq*RHS(Ztemp,tq(i));
        Ztemp = reshape(Ztemp,[],4);
        Ztemp = reshape(Ztemp(1:Ntot/4,:),[],1);
        if or(nux>0,nuv>0)
            Btemp = [sqrt(2*nux*dtq)*randn(Ntot/2,1);sqrt(2*nuv*dtq)*randn(Ntot/2,1)];
            Ztemp(1:end/2) = Ztemp(1:end/2) + Btemp(1:end/2);
            Ztemp(end/2+1:end) = Ztemp(end/2+1:end) + Btemp(end/2+1:end);
        end
        if mod(i,m)==0
            Z(i/m+1,:) = Ztemp;
            if verbose
                disp(i/m);
            end
        end
    end
end