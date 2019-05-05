% Checking parameters set for problem

% model
assert(exist('model','var')==1,...
    'ERROR: ''model'' variable must be defined')
assert(strcmp(model,'Euler'),...
    'ERROR: ''model'' must be set to ''Euler''')

% gas_const
assert(exist('gas_const','var')==1,...
    'ERROR: ''gas_const'' variable must be defined')
assert((isnumeric(gas_const) & gas_const > 0.0),...
    'ERROR: ''gas_const'' must be a positive number')

% gas_gamma
assert(exist('gas_gamma','var')==1,...
    'ERROR: ''gas_gamma'' variable must be defined')
assert((isnumeric(gas_gamma) & gas_gamma > 1.0),...
    'ERROR: ''gas_gamma'' must be a positive number')

% test_name
assert(exist('test_name','var')==1,...
    'ERROR: ''test_name'' variable must be defined')

% IC
assert(exist('InitialCond','var')==1,...
    'ERROR: ''InitialCond'' variable must be defined')

% Boundary Flags
assert(exist('BC_flags','var')==1,...
    'ERROR: ''BC_flags'' variable must be defined to set boundary conditions')
num_pfaces = length(find(not(BC_flags(:,2)-BC_ENUM.Periodic)));
assert(mod(num_pfaces,2)==0,...
    'ERROR: the number of periodic face tags should be a multiple of 2')
if(num_pfaces > 0)
    UseMeshPerData = true;
else
    UseMeshPerData = false;
end

% FinalTime
assert(exist('FinalTime','var')==1,...
    'ERROR: ''FinalTime'' variable must be defined')
assert((isnumeric(FinalTime) & FinalTime >= 0.0),...
    'ERROR: ''FinalTime'' must be a non-negative number')

% CFL and fixed_dt
assert((exist('CFL','var')==1 & exist('fixed_dt','var')==0) | ...
    (exist('CFL','var')==0 & exist('fixed_dt','var')==1),...
    'ERROR: Either ''CFL'' or ''fixed_dt'' must be defined')

if(exist('CFL','var')==1)
    assert((isnumeric(CFL) & CFL > 0.0),...
        'ERROR: ''CFL'' must be a positive number')
    fixed_dt = [];
else
    assert((isnumeric(fixed_dt) & fixed_dt > 0.0),...
        'ERROR: ''fixed_dt'' must be a positive number')
    CFL = [];
end

% tstamps
assert(exist('tstamps','var')==1,...
    'ERROR: ''tstamps'' variable must be defined')
assert((floor(tstamps)==tstamps & tstamps > 0),...
    'ERROR: ''tstamps'' must be a positive integer')

% N
assert(exist('N','var')==1,...
    'ERROR: ''N'' variable must be defined')
assert((floor(N)==N & N >= 0),...
    'ERROR: ''N'' must be a non-negative integer')

% Indicator and Limiter
assert(exist('Limiter','var')==1,...
    'ERROR: ''Limiter'' variable must be defined')
if(strcmp(Limiter,'NONE'))
    
elseif(strcmp(Limiter,'BJES') || strcmp(Limiter,'VENK'))
    
    assert(exist('Indicator','var')==1,...
        'ERROR: ''Indicator'' variable must be defined')
    switch Indicator
        
        case 'NONE'
            
        case 'ALL'
            
        case 'TVB'
            assert(exist('TVBM','var')==1,...
                'ERROR: ''TVBM'' variable must be defined since Indicator = TVB')
            assert((isnumeric(TVBM) & TVBM >= 0.0),...
                'ERROR: ''TVBM'' must be a non-negative number')
            
            assert(exist('TVBnu','var')==1,...
                'ERROR: ''TVBnu'' variable must be defined since Indicator = TVB')
            assert((isnumeric(TVBnu) & TVBnu > 1.0),...
                'ERROR: ''TVBnu'' must be a number greater than 1')
            
        case 'NN'
            
            assert(exist('nn_model','var')==1,...
                'ERROR: ''nn_model'' variable must be defined since Indicator = NN')
            
            
        otherwise
            error('Unknown indicator type %s',Indicator)
    end
    
    assert(exist('ind_var','var')==1,...
        'ERROR: ''ind_var'' variable must be defined')
    
    assert((strcmp(ind_var,'density')  | ...
        strcmp(ind_var,'velocity') | ...
        strcmp(ind_var,'pressure') | ...
        strcmp(ind_var,'prim')     | ...
        strcmp(ind_var,'con')),...
        ['ERROR: Unknown ''ind_var'' type ', ind_var])
    
    assert(exist('lim_var','var')==1,...
        'ERROR: ''lim_var'' variable must be defined')
    
    assert((strcmp(lim_var,'prim')     | ...
        strcmp(lim_var,'con')),...
        ['ERROR: Unknown ''lim_var'' type ', lim_var])
    
    if(~strcmp(Indicator,'NONE') &&  ~strcmp(Indicator,'ALL'))
        assert(exist('Filter_const','var')==1,...
            'ERROR: ''Filter_const'' variable must be defined')
        assert(islogical(Filter_const),...
            'ERROR: ''Filter_const'' must be a logical variable')
        
%         assert(exist('Remove_iso','var')==1,...
%             'ERROR: ''Remove_iso'' variable must be defined')
%         assert(islogical(Remove_iso),...
%             'ERROR: ''Remove_iso'' must be a logical variable')
        
    end
    
else
    error('Unknown indicator type %s',Indicator)
end

% Mesh file details
assert(exist('msh_file','var')==1,...
    'ERROR: ''msh_file'' variable must be defined')



% Plot and solution saving
assert(exist('plot_iter','var')==1,...
    'ERROR: ''plot_iter'' variable must be defined')
assert((floor(plot_iter)==plot_iter & plot_iter > 0),...
    'ERROR: ''plot_iter'' must be a positive integer')

assert(exist('show_plot','var')==1,...
    'ERROR: ''show_plot'' variable must be defined')
assert(islogical(show_plot),...
    'ERROR: ''show_plot'' must be a logical variable')

if(show_plot)
    
    assert(exist('plot_var','var')==1,...
        'ERROR: ''plot_var'' variable must be defined to create plots')
    assert(iscellstr(plot_var),...
        'ERROR: ''plot_var'' variable must be a cell of strings')
    for i = 1:length(plot_var)
        assert((strcmp(plot_var{i},'density')  | ...
            strcmp(plot_var{i},'velx')     | ...
            strcmp(plot_var{i},'vely')     | ...
            strcmp(plot_var{i},'pressure') | ...
            strcmp(plot_var{i},'energy')),...
            ['ERROR: Unknown variable ', plot_var{i}, ' in ''plot_var'''])
    end
    
    assert(exist('xran','var')==1,...
        'ERROR: ''xran'' variable must be defined to save plots')
    assert((isnumeric(xran) & length(xran) == 2),...
        'ERROR: ''xran'' must be a numeric array of length 2')
    
    assert(exist('yran','var')==1,...
        'ERROR: ''yran'' variable must be defined to save plots')
    assert((isnumeric(yran) & length(yran) == 2),...
        'ERROR: ''yran'' must be a numeric array of length 2')
    
    assert(exist('clines','var')==1,...
        'ERROR: ''clines'' variable must be defined to save plots')
    assert((iscell(clines) & length(clines) == length(plot_var)),...
        'ERROR: ''clines'' must be a numeric cell of size equal to that of plot_var')
    
    for i = 1:length(clines)
        assert((isnumeric(clines{i}) & length(clines{i}) >= 2),...
            'ERROR: each element of ''clines'' must be an real array of size >= 2')
    end
    
end



assert(exist('save_soln','var')==1,...
    'ERROR: ''save_soln'' variable must be defined')
assert(islogical(save_soln),...
    'ERROR: ''save_soln'' must be a logical variable')




%% Assigning data to structure object



Problem.model        = model;
Problem.gas_gamma    = gas_gamma;
Problem.gas_const    = gas_const;
Problem.test_name    = test_name;
Problem.InitialCond  = InitialCond;
Problem.FinalTime    = FinalTime;
Problem.CFL          = CFL;
Problem.fixed_dt     = fixed_dt;
Problem.tstamps      = tstamps;

Mesh.N         = N;
Mesh.msh_file  = msh_file;
Mesh.BC_ENUM   = BC_ENUM;
Mesh.BC_flags  = BC_flags;
Mesh.UseMeshPerData = UseMeshPerData;

Limit.Limiter    = Limiter;
Limit.ind_var    = [];
Limit.lim_var    = [];
if(~strcmp(Limiter,'NONE'))
    Limit.Indicator  = Indicator;
    Limit.lim_var    = lim_var;
    if(strcmp(Indicator,'TVB'))
        Limit.TVBM  = TVBM;
        Limit.TVBnu = TVBnu;
    elseif(strcmp(Indicator,'NN'))
        Limit.nn_model   = nn_model;
    end
end

if(~strcmp(Indicator,'NONE') && ~strcmp(Indicator,'ALL'))
    Limit.ind_var = ind_var;
end

Limit.Filter_const = Filter_const;
% Limit.Remove_iso   = Remove_iso;

Output.plot_iter  = plot_iter;
Output.show_plot  = show_plot;
if(show_plot)
    Output.plot_var   = plot_var;
    Output.xran       = xran;
    Output.yran       = yran;
    Output.clines     = clines;
end
Output.save_soln  = save_soln;