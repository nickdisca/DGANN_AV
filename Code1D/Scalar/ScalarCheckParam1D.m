% Checking parameters set for problem

% model
assert(exist('model','var')==1,...
    'ERROR: ''model'' variable must be defined')
assert((strcmp(model,'Advection') | strcmp(model,'Burgers') | strcmp(model,'BuckLev')),...
    'ERROR: ''model'' must be set to ''Advection'' or ''Burgers''')

% test_name
assert(exist('test_name','var')==1,...
    'ERROR: ''test_name'' variable must be defined')

% N
assert(exist('N','var')==1,...
    'ERROR: ''N'' variable must be defined')
assert((floor(N)==N & N >= 0),...
    'ERROR: ''N'' must be a non-negative integer')

% IC
assert(exist('u_IC','var')==1,...
    'ERROR: ''u_IC'' variable must be defined')

% mesh
assert(exist('bnd_l','var')==1,...
    'ERROR: ''bnd_l'' variable must be defined')
assert((isnumeric(bnd_l)),...
    'ERROR: ''bnd_l'' must be a real number')

assert(exist('bnd_r','var')==1,...
    'ERROR: ''bnd_r'' variable must be defined')
assert((isnumeric(bnd_r)),...
    'ERROR: ''bnd_r'' must be a real number')

assert(exist('mesh_pert','var')==1,...
    'ERROR: ''mesh_pert'' variable must be defined')
assert((isnumeric(mesh_pert) && mesh_pert >=0.0 && mesh_pert < 1),...
    'ERROR: ''mesh_pert'' must be a number in [0,1)')

assert(exist('K','var')==1,...
    'ERROR: ''K'' variable must be defined')
assert((floor(K)==K & K > 0),...
    'ERROR: ''K'' must be a  positive negative integer')

% Boundary Flags
assert(exist('bc_cond','var')==1,...
    'ERROR: ''bc_cond'' variable must be defined to set boundary conditions')
Check_BC1D(bc_cond,1);

% FinalTime
assert(exist('FinalTime','var')==1,...
    'ERROR: ''FinalTime'' variable must be defined')
assert((isnumeric(FinalTime) & FinalTime >= 0.0),...
    'ERROR: ''FinalTime'' must be a non-negative number')

% CFL
assert(exist('CFL','var')==1,...
    'ERROR: ''CFL'' variable must be defined')
assert((isnumeric(CFL) & CFL > 0.0),...
    'ERROR: ''CFL'' must be a positive number')

% Runge-Kutta scheme
assert(exist('RK','var')==1,...
    'ERROR: ''RK'' variable must be defined')


% Indicator and Limiter
assert(exist('Limiter','var')==1,...
    'ERROR: ''Limiter'' variable must be defined')
if(strcmp(Limiter,'NONE'))
    
elseif(strcmp(Limiter,'MINMOD'))
    
    assert(exist('Indicator','var')==1,...
        'ERROR: ''Indicator'' variable must be defined')
    switch Indicator
        
        case 'NONE'
            
        case 'ALL'
            
        case 'MINMOD'
            
        case 'TVB'
            assert(exist('TVBM','var')==1,...
                'ERROR: ''TVBM'' variable must be defined since Indicator = TVB')
            assert((isnumeric(TVBM) & TVBM >= 0.0),...
                'ERROR: ''TVBM'' must be a non-negative number')
            
        case 'NN'
            
            assert(exist('nn_model','var')==1,...
                'ERROR: ''nn_model'' variable must be defined since Indicator = NN')
            
            
        otherwise
            error('Unknown indicator type %s',Indicator)
    end
    
else
    error('Unknown limiter type %s',Limiter)
end


%Viscosity
assert(exist('Visc_model','var')==1,...
    'ERROR: ''Visc_model'' variable must be defined')
switch Visc_model
    case 'NONE'
        
    case 'MDH'
        assert(exist('c_A','var')==1 && exist('c_k','var')==1 && exist('c_max','var')==1, ...
            'ERROR: ''c_A,c_k,c_max'' variables must be defined since Visc_model = MDH')
        
    case 'MDA'
        assert(exist('c_max','var')==1, ...
            'ERROR: ''c_max'' variable must be defined since Visc_model = MDA')
        
    case 'EV'
        assert(exist('c_E','var')==1 && exist('c_max','var')==1, ...
            'ERROR: ''c_E,c_max'' variables must be defined since Visc_model = EV')
        
    case 'NN'
        assert(exist('nn_visc_model','var')==1,...
            'ERROR: ''nn_visc_model'' variable must be defined since Visc_model = NN')
        
    otherwise
            error('Unknown viscosity_model %s',Visc_model)
end


% Output flags
assert(exist('plot_iter','var')==1,...
    'ERROR: ''plot_iter'' variable must be defined')
assert((floor(plot_iter)==plot_iter & plot_iter > 0),...
    'ERROR: ''plot_iter'' must be a positive integer')

assert(exist('save_soln','var')==1,...
    'ERROR: ''save_soln'' variable must be defined')
assert(islogical(save_soln),...
    'ERROR: ''save_soln'' must be a logical variable')

assert(exist('save_ind','var')==1,...
    'ERROR: ''save_ind'' variable must be defined')
assert(islogical(save_soln),...
    'ERROR: ''save_ind'' must be a logical variable')

assert(exist('save_plot','var')==1,...
    'ERROR: ''save_plot'' variable must be defined')
assert(islogical(save_plot),...
    'ERROR: ''save_plot'' must be a logical variable')

if(save_plot)
    
    assert(exist('ref_avail','var')==1,...
        'ERROR: ''ref_avail'' variable must be defined')
    assert(islogical(ref_avail),...
        'ERROR: ''ref_avail'' must be a logical variable')
    
    assert(exist('ref_fname','var')==1,...
        'ERROR: ''ref_fname'' variable must be defined')
    
    assert(exist('var_ran','var')==1,...
        'ERROR: ''var_ran'' variable must be defined to save plots')
    Check_ran1D(var_ran,1);
    
end


%% Assigning data to structure object

Problem.model     = model;
Problem.test_name = test_name;
Problem.u_IC      = u_IC;
Problem.bc_cond   = bc_cond;
Problem.FinalTime = FinalTime;
Problem.CFL       = CFL;
Problem.RK        = RK;


Mesh.N         = N;
Mesh.bnd_l     = bnd_l;
Mesh.bnd_r     = bnd_r;
Mesh.mesh_pert = mesh_pert;
Mesh.K         = K;

Limit.Limiter    = Limiter;
Limit.Indicator  = Indicator;
if(strcmp(Indicator,'TVB'))
    Limit.TVBM = TVBM;
elseif(strcmp(Indicator,'NN'))
    Limit.nn_model   = nn_model;
end

Viscosity.model = Visc_model;
switch Visc_model     
    case 'MDH'
        Viscosity.c_A=c_A; 
        Viscosity.c_k=c_k; 
        Viscosity.c_max=c_max;
        
    case 'MDA'
        Viscosity.c_max=c_max;
        
    case 'EV'
        Viscosity.c_E=c_E; 
        Viscosity.c_max=c_max;
        
    case 'NN'
        Viscosity.nn_visc_model = nn_visc_model;
end


Output.plot_iter  = plot_iter;
Output.save_soln  = save_soln;
Output.save_ind   = save_ind;
Output.save_visc   = save_visc;
Output.save_plot  = save_plot;
if(Output.save_plot)
    Output.ref_avail  = ref_avail;
    Output.ref_fname  = ref_fname;
    Output.var_ran    = var_ran;
end
