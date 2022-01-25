
disp('Making OAASIS MRIO');

%% Options

squash_mrio = 0;
mirror_vy = 1;
iiot_conversion = 1;

options.source.phase = '002';
options.source.loop = '054';
options.base.phase = '333';
options.base.loop = '103';
options.overwrite_savedir = 1;
options.timeseries = 2009:2018;

%% Initialise

% Env 
current_dir = [fileparts(which(mfilename)) '/']; 

assert(exist([current_dir 'env.mat'],'file') > 0,'Could not discover env structure');
load([current_dir 'env.mat'],'env');

addpath(genpath(env.open_iel_tools));
addpath(genpath([env.rearrange_gloria_path 'lib/']));
addpath(genpath(current_dir));

local_store = [current_dir 'local_store/'];
create_dir_or_rebirth(local_store,0);

%% Constants

n_v = 6;
n_y = 6; 

options.root_to_base_secagg = [env.ielab_path '/Roots/GlobalIELab/ConcordanceLibrary/' 'Sector Aggregators/HSCPC_Eora25_secagg.csv']; 
options.root_to_base_regagg =  [env.ielab_path '/Roots/GlobalIELab/ConcordanceLibrary/' 'Region Aggregators/UNEP_IRP_164RegAgg.csv'];

%% Squash GLORIA MRIO

if squash_mrio
    
    disp('Squashing GLORIA MRIO');

    addpath(env.rearrange_gloria_path);
    rearrange_global_mrio(options);
    
end

%% Aggregators

%[reg_agg_src,sec_agg_src] = get_basetable_aggregators(env.basetable_path,conc_dir,phase_str_global,loop_str_global);

reg_agg = csvread(options.root_to_base_regagg);
n_reg = size(reg_agg,1);
assert(n_reg < size(reg_agg,2));

sec_agg = csvread(options.root_to_base_secagg);
n_sec = size(sec_agg,1);
assert(n_sec < size(sec_agg,2));

%% Transformations

% Table read settings
tbl_read_opts.enforce_exist = 1; 
tbl_read_opts.use_rebalanced_results  = 0;
margin = 1;
basetable_path = env.basetable_path;
phase_str = options.base.phase; 
loop_str = options.base.loop;

% Read and transform MRIO for each year
for mrio_year = options.timeseries
    
    disp(['Processing: ' num2str(mrio_year)]);
    
    % Read the MRIO
    [T,y,v,~,~] = get_basetable_blocks_ii(basetable_path,num2str(mrio_year),phase_str,loop_str,margin,tbl_read_opts);
    %[T,y,v] = get_tables(env.basetable_path,mrio_year,options.base.phase,options.base.loop,margin,tbl_read_opts);
    assert(~isempty(T) && ~isempty(v) && ~isempty(y));
    
    % Mirror global table
    if mirror_vy && (any(v(:) < 0) || any(y(:) < 0))
        
        [v_mirr,y_mirr] = apply_mirroring(v,y,n_y,n_v,n_reg,'Global');
              
        v = v_mirr;
        y = y_mirr;
        
        assert(all_positive(y) && all_positive(v));
        
        clear v_mirr y_mirr
         
        if any(T(:) < 0)
            T(T < 0) = 0;
        end
    
    end
    
    % IIOT conversion
    if iiot_conversion
        sut_to_iiot_conversion(T,v,y,n_y,n_v,n_reg,n_sec)
    end
    
    % Total output
    x_in = sum(T,1) + sum(v,1); x_out = sum(T,2) + sum(y,2);
    x_mean = mean([x_in',x_out],2);
    
    assert(all_finite(x_mean) && all_positive(x_mean));
    
    x = x_mean; x(x==0) = 10^-8;

    % Check table balance
    
    % Coefficients
    A = T*inv(diag(x)); A(A > 1) = 1;
    
    assert(all_finite(A) && min(A(:)) >= 0)
    
    % Inverse
    L = inv(eye(size(A,1)) - A);
    assert(all_finite(L));
    
    % Aggregate y
    y_agg = zeros(size(y,1),n_reg);
    for r = 1:n_reg
        idx_start = (r-1)*n_y + 1;
        idx_fin = idx_start + n_y - 1;
        y_agg(:,r) = sum(y(:,idx_start:idx_fin),2);
    end
    
    assert(is_within_tolerance(sum(y(:)),sum(y_agg(:)),0.001));
    
    % Save
    phase_str = options.base.phase;
    loop_str = options.base.loop;
    
    fname = [local_store 'A_' phase_str '_' loop_str '_' num2str(mrio_year) '.mat']; save(fname,'A','-v7.3');
    fname = [local_store 'A_' phase_str '_' loop_str '_' num2str(mrio_year) '.csv']; csvwrite(fname,A);
    
    fname = [local_store 'L_' phase_str '_' loop_str '_' num2str(mrio_year) '.mat']; save(fname,'L','-v7.3');
    fname = [local_store 'L_' phase_str '_' loop_str '_' num2str(mrio_year) '.csv']; csvwrite(fname,L);
    
    fname = [local_store 'x_' phase_str '_' loop_str '_' num2str(mrio_year) '.mat']; save(fname,'x','-v7.3');
    fname = [local_store 'x_' phase_str '_' loop_str '_' num2str(mrio_year) '.csv']; csvwrite(fname,x);
        
    fname = [local_store 'y_' phase_str '_' loop_str '_' num2str(mrio_year) '.mat']; save(fname,'y_agg','-v7.3');
    fname = [local_store 'y_' phase_str '_' loop_str '_' num2str(mrio_year) '.csv']; csvwrite(fname,y_agg);
    
end

disp('.'); disp('Finished');
