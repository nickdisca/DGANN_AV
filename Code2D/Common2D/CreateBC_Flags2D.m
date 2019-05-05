[m,n] = size(BC_cond);

BC_flags = zeros(m,n);

In = 1; Out = 2; Slip = 3; Far = 4; Dirichlet = 5; Sym = 6; Periodic=7;

BC_ENUM.In        = 1;
BC_ENUM.Out       = 2;
BC_ENUM.Slip      = 3;
BC_ENUM.Far       = 4;
BC_ENUM.Dirichlet = 5;
BC_ENUM.Sym       = 6;
BC_ENUM.Periodic  = 7;

for i=1:m
    BC_flags(i,1) = BC_cond{i,1};
    
    switch BC_cond{i,2}
        
        case 'I'
            BC_flags(i,2) = BC_ENUM.In;
            
        case 'O'
            BC_flags(i,2) = BC_ENUM.Out;
            
        case 'S'
            BC_flags(i,2) = BC_ENUM.Slip;    
            
        case 'F'
            BC_flags(i,2) = BC_ENUM.Far;
            
        case 'D'
            BC_flags(i,2) = BC_ENUM.Dirichlet;
            
        case 'Sym'
            BC_flags(i,2) = BC_ENUM.Sym;
            
        case 'P'
            BC_flags(i,2) = BC_ENUM.Periodic;    
            
        otherwise
            
            error('Unknown bc type %s',BC_cond{i,2})
            
    end
end
            