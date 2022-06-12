
function [mt,pt,sig_est,corners,ufft] = findcorners_sde(X,t,tau,tauhat,max_dt)

T = length(t);
l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*log(tau);
[corners,sig_est,ufft] = findcornerpts(X,t);
k = corners(2);
mstar1 = sqrt(3)/pi*T/2/k*tauhat;
mstar2 = 1/pi*tauhat*(T/2)/k*sqrt(log(exp(1)^3/tau^8));
mt = fzero(@(m)l(m,k,T), [mstar1 mstar2]);
if mt>T/2-1
    mt = T/2-1;
end
mt = min(floor((T-1)/2),ceil(mt));
pt = max(max_dt,ceil(log(tau)/log(1-(1-1/mt)^2)));

end

function [corners,sig_est,halfufft] = findcornerpts(X,t)
    dims = size(X);
    t = t';
    L = length(t);
    wn = ((0:L-1)-floor(L/2))'*(2*pi)/range(t);
    tt = wn(1:ceil(end/2));
    NN = length(tt);
    Ufftfull = zeros(dims(2)*dims(1),dims(3));
    for i=1:dims(2)
        Ufftfull((1:dims(1))+(i-1)*dims(1),:) = abs(fftshift(fft(reshape(X(:,i,:),dims(1),dims(3)),[],2)));
    end
    Ufft = mean(Ufftfull(:,1:ceil((L+1)/2)));
    corners=findchangepts(fliplr(log(Ufft)));
    corners=[-tt(corners) corners];
    sig_est = sqrt(mean(reshape(Ufftfull(:,[1:NN-corners(2)-1 NN+corners(2)+1:end]),[],1).^2)/L);
    halfufft=fliplr(log(Ufft));

%     Ufft = cumsum(Ufft);            
%     errs = zeros(NN-2,1);
%     for k=2:NN-1
%        subinds1 = 1:k;
%        subinds2 = k:NN;
%        Ufft_av1 = Ufft(subinds1);
%        Ufft_av2 = Ufft(subinds2);
%        m1 = range(Ufft_av1)/range(tt(subinds1));
%        m2 = range(Ufft_av2)/range(tt(subinds2));
%        L1 = min(Ufft_av1)+m1*(tt(subinds1)-tt(1));
%        L2 = max(Ufft_av2)+m2*(tt(subinds2)-tt(end));
%        errs(k-1) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2)); % relative l2 
%     end
%     [~,tstarind] = min(errs);
%     tstar = -tt(tstarind);
%     corners = [tstar max(NN-tstarind-3,1)]
%     Ufft = abs(fftshift(fft(reshape(X(:,1,:),dims(1),dims(3))')));
%     sig_est = sqrt(mean(reshape(Ufft([1:tstarind-3 2*NN-tstarind+2:end],:),[],1).^2)/L);
end
