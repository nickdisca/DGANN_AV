% Checking parameters set for problem

% model
assert(exist('model','var')==1,...
    'ERROR: ''model'' variable must be defined')
assert((strcmp(model,'Advection') | strcmp(model,'Burgers')),...
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
assert(exist('InitialCond','var')==1,...
    'ERROR: ''InitialCond'' variable must be defined')

% Periodic Flag
assert(exist('UseMeshPerData','var')==1,...
    'ERROR: ''UseMeshPerData'' variable must be defined')
assert(islogical(UseMeshPerData),...
    'ERROR: ''UseMeshPerData'' must be a logical variable')

% FinalTime
assert(exist('FinalTime','var')==1,...
    'ERROR: ''FinalTime'' variable must be defined')
assert((isnumeric(FinalTime) & FinalTime >= 0.0),...
    'ERROR: ''FinalTime'' must be a non-negative number')

% Advection Velocity
if(strcmp(model,'Advection'))
    assert(exist('AdvectionVelocity','var')==1,...
        'ERROR: ''AdvectionVelocity'' variable must be defined')
    assert((isnumeric(AdvectionVelocity) & length(AdvectionVelocity) ==2),...
        'ERROR: ''AdvectionVelocity'' must be numeric array of length 2')
end

% Indicator and Limiter
assert(exist('Limiter','var')==1,...
    'ERROR: ''Limiter'' variable must be defined')
if(strcmp(Limiter,'NONE'))
    
elseif(strcmp(Limiter,'BJES') || strcmp(Limiter,'VENK'))
    
    assert(exist('Indicator','var')==1,...
        'ERROR: ''Indicator'' variable must be defined')
    switch Indicator
        
        case 'NONE'
        
        case 'TVB'
            assert(exist('TVBM','var')==1,...
                'ERROR: ''TVBM'' variable must be defined since Indicator = TVB')
            assert((isnumeric(TVBM) & TVBM >= 0.0),...
                'ERROR: ''TVBM'' must be a non-negative number')
            
            assert(exist('TVBnu','var')==1,...
                'ERROR: ''TVBnu'' variable must be defined since Indicator = TVB')
            assert((isnumeric(TVBnu) & TVBnu > 1.0),...
                'ERROR: ''TVBnu'' must be a number greater than 1')
            
        case 'TVB2'
            assert(exist('TVBM','var')==1,...
                'ERROR: ''TVBM'' variable must be defined since Indicator = TVB2')
            assert((isnumeric(TVBM) & TVBM >= 0.0),...
                'ERROR: ''TVBM'' must be a non-negative number')
            
            assert(exist('TVBnu','var')==1,...
                'ERROR: ''TVBnu'' variable must be defined since Indicator = TVB2')
            assert((isnumeric(TVBnu) & TVBnu > 1.0),...
                'ERROR: ''TVBnu'' must be a number greater than 1')    
            
        case 'NN'

            assert(exist('nn_model','var')==1,...
                'ERROR: ''nn_model'' variable must be defined since Indicator = NN')
            
        case 'NN_Pwise'
            assert(exist('nn_model','var')==1,...
                'ERROR: ''nn_model'' variable must be defined since Indicator = NN_Pwise')    
        
        case 'NN_modal_Pwise'
            assert(exist('nn_model','var')==1,...
                'ERROR: ''nn_model'' variable must be defined since Indicator = NN_modal_Pwise') 
            
        case 'NN_modal_patch_Pwise'
            assert(exist('nn_model','var')==1,...
                'ERROR: ''nn_model'' variable must be defined since Indicator = NN_modal_patch_Pwise')     
            
        otherwise
            error('Unknown indicator type %s',Indicator)
    end
    
else
    error('Unknown indicator type %s',Indicator)
end

% Mesh file details
assert(exist('mtail','var')==1,...
    'ERROR: ''mtail'' variable must be defined')
assert(exist('msh_file','var')==1,...
    'ERROR: ''msh_file'' variable must be defined')

% Plot saving
assert(exist('Save_plots','var')==1,...
    'ERROR: ''Save_plots'' variable must be defined')
assert(islogical(Save_plots),...
    'ERROR: ''Save_plots'' must be a logical variable')

if(Save_plots)

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