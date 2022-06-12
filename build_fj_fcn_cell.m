function [fj_fcn_cell,fj_tags,rad_c] = build_fj_fcn_cell(J,pdist,fj_basis,rad_cfac,tol)
    if J>0
        fj_fcn_cell = cell(J,1);
        if isequal(fj_basis,'NH') % Newton-Hermite
            rad_c = tol*max(pdist(2,:));
            rmin = 10^-3;
            fj_tags = [[[(0:ceil(J/2)-1)';floor(J/2)*ones(floor(J/2),1)] [zeros(ceil(J/2),1);(1:floor(J/2))']]];
            for j=1:J
                pow = fj_tags(j,1);
                if pow<0
                    fj = @(x) (x~=0).*max(rmin,x).^pow;
                else
                    pow2 = fj_tags(j,2);
                if and(pow2 <= 0,pow==0)
                    fj = @(x) x*0+1;
                elseif pow2>=0
                    fj = @(x) (x.^pow).*(x-rad_c).^pow2;
                end
                end
                fj_fcn_cell{j} = @(x) fj(x).*(x<=rad_c)/sqrt(trapz(pdist(2,:),pdist(1,:).*fj(pdist(2,:)).^2));
            end
        elseif isequal(fj_basis,'BP') % bernstein polynomial
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:(J-1); fj_tags = [fj_tags' -fj_tags']+[0 J-1];
            for j=1:J
                pows = fj_tags(j,:);
                fj_fcn_cell{j} = @(x) nchoosek(J-1,pows(1))*((x/rad_c).^pows(1)).*(1-(x/rad_c)).^pows(2)/(1/2)^(J-1)/nchoosek(J-1,floor((J-1)/2)).*(x<=rad_c);
    %            fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/sqrt(trapz(pdist(2,:),pdist(1,:).*fj_fcn_cell{j}(pdist(2,:))));
            end
        elseif isequal(fj_basis,'MO') % morse                
            rad_c = tol;
%             alph = 0.95;
%             f = @(lam) (1-exp(-lam*rad_c)).^2-alph^2*rad_c/2*lam.*(1-exp(-2*lam*rad_c));
%             lam0 = fzero(f,[0.00001 2/rad_c]);
%             lambdas = [0 lam0.*rad_cfac.^(0:J-2)];
            if length(rad_cfac)==1
                lambdas = [0:rad_cfac:(J-1)*rad_cfac];
            else
                lambdas=rad_cfac;
            end
            fj_fcn_cell{1} = @(x) x<=tol;
%             fj_fcn_cell{1} = @(x) fj_fcn_cell{1}(x)/sqrt(trapz(pdist(2,:),pdist(1,:)));
            fj_tags = lambdas;
            for j=2:J
                fj_fcn_cell{j} = @(x) exp(-x*lambdas(j))*lambdas(j).*(x<=tol);
            end
            for j=1:J
                fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/sqrt(trapz(pdist(2,:),pdist(1,:).*fj_fcn_cell{j}(pdist(2,:)).^2));
            end
        elseif isequal(fj_basis,'B2') % bernstein polynomial, smooth, compact support on [0 rad_c)
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:J; fj_tags = [fj_tags' -fj_tags']+[0 J];
            fj_tags = fj_tags(1:J,:);
            for j=1:J
                pows = fj_tags(j,:);
                fj_fcn_cell{j} = @(x) nchoosek(J,pows(1))*((x/rad_c).^pows(1)).*(1-(x/rad_c)).^pows(2)/(1/2)^J/nchoosek(J,floor(J/2));
    %            fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/abs(sqrt(trapz(pdist(2,:),pdist(1,:).*fj_fcn_cell{j}(pdist(2,:)))));
            end
    %     elseif fj_basis == 'PD' % particle distance
    %         R = zeros(2*J-1,1);
    %         for j=1:2*J-1
    %             R(j) = mean(reshape(Xrad.^(j-1),[],1));
    %         end
    %         for j=1:J
    %             fj_fcn_cell{j} = @(x) x.^(j-1);
    %             for k=1:j-1
    %                 fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)-R(j+k-1)/sqrt(R(2*k-1))*fj_fcn_cell{k}(x);
    %             end
    %             fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/sqrt(R(2*j-1));
    %         end 
    %         fj_tags = [];
        elseif isequal(fj_basis,'PD') % particle distance
            rad_c = tol*max(pdist(2,:));
            fj_fcn_cell{1} = @(x) x*0+1;
            fj_fcn_cell{2} = @(x) x-1;
            U = sqrt(trapz(pdist(2,:),fj_fcn_cell{2}(pdist(2,:)).^2.*pdist(1,:)));
            fj_fcn_cell{2} = @(x) fj_fcn_cell{2}(x)/U;
            for j=3:J
                fj_fcn_cell{j} = @(x) x.^(j-1);
                for k=1:j-1
                    U = trapz(pdist(2,:),fj_fcn_cell{k}(pdist(2,:)).*fj_fcn_cell{j}(pdist(2,:)).*pdist(1,:));
                    fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)-U*fj_fcn_cell{k}(x);
                end
                U = sqrt(trapz(pdist(2,:),fj_fcn_cell{j}(pdist(2,:)).^2.*pdist(1,:)));
                fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/U;
            end 
            fj_tags = [];
        elseif isequal(fj_basis,'HF') % particle distance
            ee = pdist(2,:);
            vv = cumsum(pdist(1,:));
            vv = vv/vv(end);
            eemax = ee(find(vv>=tol,1));
            if isempty(eemax)
                eemax = max(ee);
            end
            if rad_cfac == 0
                rad_c = eemax/J;
                phi = @(r) and(r>=-rad_c/2,r<rad_c/2);
                c = linspace(ee(1),eemax,J+1);
                cL = c(1:end-1);
                cR = c(2:end);
            elseif rad_cfac == 1
                rad_c = eemax/J;
                phi = @(r) max(1-abs(r/rad_c),0);
                c = linspace(ee(1),eemax,J+1);
                cL = c(1:end-1);
                cR = c(2:end);
            else
                pow = rad_cfac-1;
                rad_c = eemax/(2*(1+J*sqrt(1-2^(-1/pow))));
                t = 2*rad_c*sqrt(1-2^(-1/pow));
                phi = @(r) max((1-(r/rad_c).^2),0).^pow;
                cL = linspace(ee(1),(J-1)*t,J);
                cL = 0:t:(J-1)*t;
                cR = cL+2*rad_c;
            end
            I = integral(@(r)phi(r).^2,-rad_c,rad_c);
            minp = min(pdist(1,pdist(1,:)>0));
            indsL = cL*0;
            indsR = cR*0;
            for j=1:J
                [~,indsL(j)] = min(abs(ee-cL(j)));
                [~,indsR(j)] = min(abs(ee-cR(j)));
    %            indsR(j) = find(ee>=cR(j),1);
            end
            for j=1:J
                fj_fcn_cell{j} = @(r) phi(r-(cL(j)+cR(j))/2)/sqrt(max(mean(pdist(1,indsL(j):indsR(j))),minp)*2*rad_c*I);
            end            
            fj_tags = [];
            rad_c = eemax;
        elseif isequal(fj_basis,'PC') % particle distance
            minp = min(pdist(1,pdist(1,:)>0));
            fj_tags = [];
            ee = pdist(2,:);
            vv = cumsum(pdist(1,:));
            vv = vv/vv(end);
            vvmax = vv(find(vv>=tol,1));
            eemax = ee(find(vv>=tol,1));
            if isempty(vvmax)
                vvmax=max(vv);
            end
            U = linspace(vv(1),vvmax,J+1+(rad_cfac>0));
            inds = U*0;
            inds(1) = 1;
            for j=2:length(U)
                inds(j) = find(vv>=U(j),1);
                if inds(j)<=inds(j-1)
                    inds(j) = inds(j-1)+1;
                end
            end
            if rad_cfac ==0                
                phi = @(r,a,b) and(r>=a,r<b);
            elseif rad_cfac ==1
                phi = @(r,a,b) max(1-abs((r-mean([a b]))/(b-a)),0);
            else
                phi = @(r,a,b) max(1-((r-mean([a b]))/(b-a)).^2,0).^(rad_cfac-1);
            end
            for j=1:J
                if rad_cfac>0
                    fj_fcn_cell{j} = @(r) phi(r,ee(inds(j)),ee(inds(j+2)))/sqrt((vv(inds(j+2))-vv(inds(j))));
                else
                    fj_fcn_cell{j} = @(r) phi(r,ee(inds(j)),ee(inds(j+1)))/sqrt((vv(inds(j+1))-vv(inds(j))));
                end
            end
            rad_c = eemax;
        elseif isequal(fj_basis,'DG') % particle distance
            ee = pdist(2,:);
            vv = cumsum(pdist(1,:));
            vv = vv/vv(end);
            eemax = ee(find(vv>=tol,1));
            if isempty(eemax)
                eemax = max(ee);
            end
            c = linspace(0,eemax,J+1);
            rad_c = eemax/(J/p);
            phi0 = @(r) (abs(r)<=rad_c/2)/sqrt(rad_c);
            phil = @(r) sqrt(3/rad_c^3)*r.*(abs(r)<=rad_c/2);
            phi2 = @(r) sqrt(5/2/rad_c)*(3/2*(r/rad_c).^2-1).*(abs(r)<=rad_c/2);
            inds = c*0;
            for j=1:J+1
                inds(j) = find(ee>=c(j),1);
            end
            for j=1:2:J-1
                fj_fcn_cell{j} = @(r) phic(r-(c(j+1)+c(j))/2)/sqrt(pdist(1,inds(j))*rad_c);
                fj_fcn_cell{j+1} = @(r) phil(r-(c(j+1)+c(j))/2)/sqrt(pdist(1,inds(j))*rad_c)/sqrt((rad_c/2)^3*2/3);
            end
            if mod(J,2)~=0
                fj_fcn_cell{J} = @(r) r*0;
            end
            fj_tags = [];
        elseif isequal(fj_basis,'dlin') % positive damping basis
            rad_c = tol*max(pdist(2,:));
            fj_tags = linspace(0,rad_c,J);
            fj_fcn_cell{1} = @(x) x*0+1;
            for j=2:J
                fj_fcn_cell{j} = @(x) min(max(1/fj_tags(j)*x+1,0),2);%/sqrt(trapz(pdist(2,:),pdist(1,:).*fj(pdist(2,:)).^2));
            end
        elseif isequal(fj_basis,'BPd') % bernstein polynomial
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:(J-1); fj_tags = [fj_tags' -fj_tags']+[0 J-1];
            for j=1:J
                pows = fj_tags(j,:);
                fj_fcn_cell{j} = @(x) nchoosek(J-1,pows(1))*(((x+rad_c)/2/rad_c).^pows(1)).*(1-((x+rad_c)/2/rad_c)).^pows(2)/(1/2)^(J-1)/nchoosek(J-1,floor((J-1)/2)).*(abs(x)<=rad_c);
            end
        elseif isequal(fj_basis,'mon')
            rad_c = tol;%*max(pdist(2,:));
            fj_tags = 0:(J-1);
            for j=1:J
                fj_fcn_cell{j} = @(x) (x/rad_c/2).^fj_tags(j);%.*(x<=rad_c);
            end
            for j=1:J
                fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/sqrt(trapz(pdist(2,:),pdist(1,:).*fj_fcn_cell{j}(pdist(2,:)).^2));
            end
        elseif isequal(fj_basis,'bess1') % bernstein polynomial
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:(J-1);
            for j=1:J
                fj_fcn_cell{j} = @(x) (x/rad_c).^fj_tags(j);
            end
        elseif isequal(fj_basis,'cos') % shifted cosine series
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:(J-1);
            fj_fcn_cell{1} = @(x) 0*x+1;
            for j=2:J
                fj_fcn_cell{j} = @(x) cos((j-1)*x)+1;
            end
            for j=1:J
                fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/sqrt(3*pi);
            end
        elseif isequal(fj_basis,'cos2') % shifted cosine series
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:(J-1);
            fj_fcn_cell{1} = @(x) 0*x+1;
            for j=2:J
                fj_fcn_cell{j} = @(x) cos((j-1)*x);
            end
            for j=1:J
                fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x);
            end
        elseif isequal(fj_basis,'trig') % shifted cosine series
            rad_c = tol*max(pdist(2,:));
            fj_tags = 0:(J-1);
            fj_fcn_cell{1} = @(x) 0*x+1;
            for j=1:(J-1)/2
                fj_fcn_cell{2*j} = @(x) cos(j*x)+1;
                fj_fcn_cell{2*j+1} = @(x) -cos(j*x)+1;
            end
            for j=1:J
                fj_fcn_cell{j} = @(x) fj_fcn_cell{j}(x)/sqrt(3*pi);
            end
        elseif isequal(fj_basis,'laguerre') % shifted cosine series
            rad_c = tol*rad_cfac;
%             syms a;
            for j=1:J
%                 fj_fcn_cell{j} = matlabFunction(36/rad_c*laguerreL(j-1,a/rad_c*36)*exp(-a/rad_c/2*36));
                c = 36/rad_c*LaguerrePoly(j-1);
                fj_fcn_cell{j} = @(x) horn(c,x/rad_c*36).*exp(-x/rad_c/2*36).*(x<=tol);
                fj_tags{j}=c;
            end
        end
    else
        fj_fcn_cell = {};
        fj_tags = {};
        rad_c = 0;
    end
end
