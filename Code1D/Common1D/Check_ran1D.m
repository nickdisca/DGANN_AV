function Check_ran1D(var_ran,dim)

[m,n] = size(var_ran);

if(m~=dim || n ~=2)
    error('Assigned var_ran does not have the right dimension!!')
end

return
