function fn = findfile(d,str,defaultnum)
    if isempty(defaultnum)
        defaultnum=1;
    end
    if isempty(d)
        d = dir;
    else
        d=dir(d);
    end
    
    s=cell(length(d),1);
    [s{:}]=d.name;
    inds=find(cellfun(@(x) contains(x,str),s));
    
    if ~isempty(inds)
        fn =s{inds(defaultnum)};
    else
        fn = '';
        disp('no file found')
    end
end
