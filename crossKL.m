function [KLxxmat,KLxx,eh,h1]=crossKL(X,Y,i1,i2,tpts,numbins,Rcap)
    N2=length(i2);
    X=X(:,:,tpts);
    Y=Y(:,:,tpts);
    R = min(min([sqrt(max(range(X(:,1,:),1))^2+max(range(X(:,2,:),1))^2);...
        sqrt(max(range(Y(:,1,:),1))^2+max(range(Y(:,2,:),1))^2)]),Rcap);
    e=linspace(0,R,numbins+1);
    eh=(e(1:end-1)+e(2:end))/2;

    h1 = crossparticle_dist(X(i1,:,:),X,e);
    KLxx=zeros(N2,1);
    KLxxmat=zeros(numbins,N2);

    for i=1:N2
        h2 = crossparticle_dist(Y(i2(i),:,:),Y,e);
        supph12=and(h1,h2);
        KLxx(i)=-sum(h1(supph12).*log(h2(supph12)./h1(supph12)));    
        KLxxmat(:,i)=h2(:);
    end
end