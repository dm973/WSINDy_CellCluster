function [KLvmat,KLv,eh,h1]=crossvel(X,Y,i1,i2,tpts,numbins,Rcap)
    N2=length(i2);
    X=X(:,:,tpts);
    Y=Y(:,:,tpts);
    R = min(Rcap,min(max(max(vecnorm(X(i1,:,:),2,2))),max(max(vecnorm(Y(i2,:,:),2,2)))));
    e=linspace(0,R,numbins+1);
    eh=(e(1:end-1)+e(2:end))/2;

    h1 = histcounts(vecnorm(X(i1,:,:),2,2),e);
    KLv=zeros(N2,1);
    KLvmat=zeros(numbins,N2);

    for i=1:N2
        h2 = histcounts(vecnorm(Y(i2(i),:,:),2,2),e);
        supph12=and(h1,h2);
        KLv(i)=-sum(h1(supph12).*log(h2(supph12)./h1(supph12)));
        KLvmat(:,i)=h2(:);
    end
end
