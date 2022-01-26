function idxs = get_mrio_index(r,n_sec,n_y,n_v,format,entity)
%get_mrio_index
    
    if strcmp(format,'sut') && strcmp(entity,'industry')
        idxs = (r-1)*n_sec*2 + 1 : (r-1)*n_sec*2 + n_sec;
    elseif strcmp(format,'sut') && strcmp(entity,'commodity')
        idxs = (r-1)*n_sec*2 + n_sec + 1 : (r-1)*n_sec*2 + n_sec*2;
    elseif strcmp(format,'iiot') && strcmp(entity,'industry')
        idxs = (r-1)*n_sec + 1 : (r-1)*n_sec + n_sec;
    elseif strcmp(entity,'fd')
         idxs = (r-1)*n_y + 1 : (r-1)*n_y + n_y;
    elseif strcmp(entity,'va')
         idxs = (r-1)*n_v + 1 : (r-1)*n_v + n_v;
    else
        error('Unknown context');
    end

end

