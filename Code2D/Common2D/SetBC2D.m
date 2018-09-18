function BCind = SetBC2D(face_list,EToV)

BCind  = false(size(EToV));
Nfaces = 3; 
for j = 1:length(face_list(:,1))
    
    f1                 = face_list(j,:);
    [elem_listA,dummy] = find(not(EToV - f1(1)));
    [elem_listB,dummy] = find(not(EToV - f1(2)));
    elem               = intersect(elem_listA,elem_listB);
            
    f_ind              = find(EToV(elem,:)==f1(1));
    if(~(EToV(elem,1+mod(f_ind,Nfaces))==f1(2)))
        f_ind  = 1+mod(f_ind+1,Nfaces);
    end
    
    BCind(elem,f_ind)  = true;
    
end




return