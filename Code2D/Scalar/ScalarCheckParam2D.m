% Checking parameters set for problem

% model
assert(exist('model','var')==1,...
    'ERROR: ''model'' variable must be defined')
assert((strcmp(model,'Advection') | strcmp(model,'Burgers')),...
    'ERROR: ''model'' must be set to ''Advection'' or ''Burgers''')

% Advection Velocity
if(strcmp(model,'Advection'))
    assert(exist('AdvectionVelocity','var')==1,...
        'ERROR: ''AdvectionVelocity'' variable must be defined')
    assert((isnumeric(AdvectionVelocity) & length(AdvectionVelocity) ==2),...
        'ERROR: ''AdvectionVelocity'' must be numeric array of length 2')
end

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
            
            %
            %         case 'NN_modal_patch_Pwise'
            %             assert(exist('nn_model','var')==1,...
            %                 'ERROR: ''nn_model'' variable must be defined since Indicator = NN_modal_patch_Pwise')
            
        otherwise
            error('Unknown indicator type %s',Indicator)
    end
    
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
assert(exist(msh_file,'file')==2,...
    'ERROR: Specified ''msh_file'' not found')


% Plot and solution saving
assert(exist('plot_iter','var')==1,...
    'ERROR: ''plot_iter'' variable must be defined')
assert((floor(plot_iter)==plot_iter & plot_iter > 0),...
    'ERROR: ''plot_iter'' must be a positive integer')

assert(exist('save_soln','var')==1,...
    'ERROR: ''save_soln'' variable must be defined')
assert(islogical(save_soln),...
    'ERROR: ''save_soln'' must be a logical variable')

assert(exist('show_plot','var')==1,...
    'ERROR: ''show_plot'' variable must be defined')
assert(islogical(show_plot),...
    'ERROR: ''show_plot'' must be a logical variable')

if(show_plot)
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
    assert((isnumeric(clines) & length(clines) >=2),...
        'ERROR: ''clines'' must be a numeric array of size >= 2')
end


%% Assigning data to structure object



Problem.model     = model;
if(strcmp(model,'Advection'))
    Problem.AdvectionVelocity = AdvectionVelocity;
end
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
if(~strcmp(Limiter,'NONE'))
    Limit.Indicator  = Indicator;
    if(strcmp(Indicator,'TVB'))
        Limit.TVBM  = TVBM;
        Limit.TVBnu = TVBnu;
    elseif(strcmp(Indicator,'NN'))
        Limit.nn_model   = nn_model;
    end
end
Limit.ind_var      = [];
Limit.lim_var      = [];
Limit.Filter_const = Filter_const;
%Limit.Remove_iso   = Remove_iso;



Output.plot_iter  = plot_iter;
Output.save_soln  = save_soln;
Output.show_plot  = show_plot;
if(show_plot)
    Output.xran       = xran;
    Output.yran       = yran;
    Output.clines     = clines;
end

