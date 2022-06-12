function oi=oppinds(inds,n)
    oi=1:n;
    oi=find(~ismember(oi,inds));
end