function fn = findfilestrloc(d,str,defaultind)
    if isempty(defaultind)
        defaultind=1;
    end
    if isempty(d)
        d = dir;
    else
        d=dir(d);
    end
    
    s=cell(length(d),1);
    [s{:}]=d.name;
    inds=find(cellfun(@(x) isequal(x(min(length(x),defaultind:defaultind-1+length(str))),str),s));
    
    if ~isempty(inds)
        if length(inds)==1
            fn =s{inds};
        else
            fn =s(inds);
            disp('multiple files found');
        end
    else
        fn = '';
        disp('no file found');
    end
end