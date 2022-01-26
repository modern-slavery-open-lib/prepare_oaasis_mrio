function [T_iiot,v_iiot, y_iiot] = sut_to_iiot_conversion(T,v,y,n_y,n_v,n_reg,n_sec)
%sut_to_iiot_conversion

    disp('. Performing SUT to IIOT conversion');
    
    % Preallocate
    T_iiot = zeros(n_reg*n_sec,n_reg*n_sec);
    y_iiot = zeros(n_reg*n_sec,n_reg*n_y);
    v_iiot = zeros(n_reg*n_v,n_reg*n_sec);

    % Rewrite
    for r = 1:n_reg
        
        % Origin indices
        idx_i_o_sut = get_mrio_index(r,n_sec,n_y,n_v,'sut','industry'); 
        idx_c_o_sut = get_mrio_index(r,n_sec,n_y,n_v,'sut','commodity');
        idx_v_o = get_mrio_index(r,n_sec,n_y,n_v,'sut','va');
        idx_i_o_iiot = get_mrio_index(r,n_sec,n_v,n_y,'iiot','industry'); 
        
        assert(isempty(intersect(idx_i_o_sut,idx_c_o_sut)));

        for s = 1:n_reg
            
            % Dest indices
            idx_i_d_sut = get_mrio_index(s,n_sec,n_y,n_v,'sut','industry'); 
            idx_c_d_sut = get_mrio_index(s,n_sec,n_y,n_v,'sut','commodity');
            idx_f_d = get_mrio_index(s,n_sec,n_y,n_v,'sut','fd');
            idx_i_d_iiot = get_mrio_index(s,n_sec,n_y,n_v,'iiot','industry'); 
            
            assert(isempty(intersect(idx_i_d_sut,idx_c_d_sut)));
            
            if n_dim_collapse(T(idx_c_o_sut,idx_i_d_sut)) > 0
            
                % Supply matrix (or importing country's supply block)
                V = T(idx_i_d_sut,idx_c_d_sut);

                % Patch zeros in supply
                if any(sum(V,2) == 0)
                    mask = eye(size(V,1));
                    nnz_rows = find(sum(V,2) > 0);
                    mask(nnz_rows,:) = 0; mask(mask > 0) = 0.1*min(V(V>0));
                    V = V + mask;
                    assert(sum(mask(:)) > 0 && ~any(sum(V,2) == 0));
                end

                % Transactions matrix
                if r ~= s % Trade blocks

                    % Use matrix
                    U = T(idx_c_o_sut,idx_i_d_sut);

                    % Make industry-by-indutry transactions matrix
                    g = sum(V,2);
                    Tr = pinv(diag(g))*V; % transformation matrix
                    S = U*Tr; 
                    S = S * sum(U(:))/sum(S(:));

                else % Domestic block - convert SUT to CIOT

                    % Use matrix (the origin and dest labels are the same for domestic case)
                    U = T(idx_c_o_sut,idx_i_o_sut);

                    assert(sum(U(:)) > 0.001);

                    % Make industry-by-indutry transactions matrix
                    g = sum(V,2); 
                    Tr = pinv(diag(g))*V; % transformation matrix
                    S = U*Tr; 
                    S = S * sum(U(:))/sum(S(:));
                end

                assert(all_finite(S) && all_positive(S) && all(size(S) == n_sec));
                assert(is_within_tolerance(sum(S(:)),sum(U(:)),0.001));
                if n_dim_collapse(T(union(idx_i_o_sut,idx_c_o_sut),union(idx_i_d_sut,idx_c_d_sut))) > 0
                    assert(sum(S(:)) > 0);
                end
                assert(n_dim_collapse(T_iiot(idx_i_o_iiot,idx_i_d_iiot)) == 0);

                T_iiot(idx_i_o_iiot,idx_i_d_iiot) = S;
                
            end
            
            % Final demand
            assert(n_dim_collapse(y_iiot(idx_i_o_iiot,idx_f_d)) == 0);

            y_iiot(idx_i_o_iiot,idx_f_d) = y(idx_c_o_sut,idx_f_d);
            
            % VA
            assert(n_dim_collapse(v_iiot(idx_v_o,idx_i_d_iiot)) == 0);
            v_iiot(idx_v_o,idx_i_d_iiot) = v(idx_v_o,idx_i_d_sut) + v(idx_v_o,idx_c_d_sut);
                
        end
    end
    
    % Tests
    assert(all_finite(T_iiot) && all_finite(v_iiot) && all_finite(y_iiot));
    assert(is_within_tolerance(sum(y(:)),sum(y_iiot(:)),0.0001));
    assert(is_within_tolerance(sum(v(:)),sum(v_iiot(:)),0.0001));
    
    r = randi([1 n_reg]);
    idxs_sut = get_mrio_index(r,n_sec,n_y,n_v,'sut','commodity');
    x_sut = sum(T(idxs_sut,:),2) + sum(y(idxs_sut,:),2); 
    
    idxs_iiot = get_mrio_index(r,n_sec,n_y,n_v,'iiot','industry');
    x_iiot = sum(T_iiot(idxs_iiot,:),2) + sum(y_iiot(idxs_iiot,:),2); 
    
    assert(is_within_tolerance(sum(x_sut),sum(x_iiot),0.001));

end

