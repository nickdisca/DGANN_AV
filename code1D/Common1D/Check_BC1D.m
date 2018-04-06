function Check_BC1D(bc_cond,dim)

[m,n] = size(bc_cond);

if(m~=dim || n ~=4)
    error('Assigned bc_cond is not valid!!')
end

vbc_type = ['P','D','N'];

for i=1:dim
    if(~ismember(bc_cond{i,1},vbc_type) || ~ismember(bc_cond{i,3},vbc_type))
        error('bc_types in bc_cond can only be ''P'', ''D'' or ''N''!!');
    end
    if(~isfloat(bc_cond{i,2}) || ~isfloat(bc_cond{i,4}))
        error('bc_val in bc_cond can only be floats!!');
    end
end

return
