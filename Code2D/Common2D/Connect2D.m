function [EToE, EToF, PShift, BCTag] = Connect2D(EToV,BFaces,PerBToB_map,PerBFToF_map,BC_flags,...
                                  UseMeshPerData,VX,VY,Pflag)

% function [EToE, EToF] = Connect2D(EToV,BFaces,PerBToB_map,PerBFToF_map,UseMeshPerData)
% Purpose  : Build global connectivity arrays for grid based on
%            standard EToV input array from grid generator
%            Periodic conditions also considered (if available)

Nfaces = 3;

% Find number of elements and vertices
K = size(EToV,1); Nv = max(max(EToV));

% Create face to node connectivity matrix
TotalFaces = Nfaces*K;

% List of local face to local vertex connections
vn = [[1,2];[2,3];[1,3]];

% Build global face to node sparse array
% tic
% SpFToV = spalloc(TotalFaces, Nv, 2*TotalFaces);
% sk = 1;
% for k=1:K
%   for face=1:Nfaces
%     SpFToV( sk, EToV(k, vn(face,:))) = 1;
%     sk = sk+1;
%   end
% end
% toc

i       = [(1:K*Nfaces);(1:K*Nfaces)];
i       = i(:);
j       = [EToV(:,1)';EToV(:,2)';EToV(:,2)';EToV(:,3)';EToV(:,1)';EToV(:,3)'];
j       = j(:);
vals    = ones(2*TotalFaces,1);
SpFToV  = sparse(i,j,vals);


% Build global face to global face sparse array
SpFToF = SpFToV*SpFToV' - 2*speye(TotalFaces);

% Find complete face to face connections
[faces1, faces2] = find(SpFToF==2);

% Convert face global number to element and face numbers
element1 = floor( (faces1-1)/Nfaces )  + 1; face1    =   mod( (faces1-1), Nfaces ) + 1;
element2 = floor( (faces2-1)/Nfaces )  + 1; face2    =   mod( (faces2-1), Nfaces ) + 1;

% Rearrange into Nelements x Nfaces sized arrays
ind = sub2ind([K, Nfaces], element1, face1);

EToE = (1:K)'*ones(1,Nfaces); EToF = ones(K,1)*(1:Nfaces);
EToE(ind) = element2; EToF(ind) = face2;

% Create periodic cell data if available AND requested

% PShift is a map whose key is the index of a triangle with a periodic face
% and the value is an array of corresponding triangle pair and the
% (x,y) distance from the former. Note that some triangles may contain two
% periodic pairs. However, this usually happens at corners.
PShift = containers.Map('KeyType','uint32','ValueType','any');
if(UseMeshPerData)
    
    % Loop of periodic faces
    masterface = keys(PerBToB_map);
    for i=1:PerBToB_map.Count
        bf1 = masterface{i};
        bf2 = PerBToB_map(bf1);
        
        face_list1    = BFaces(bf1);
        face_list2    = BFaces(bf2);
        pflink        = PerBFToF_map(bf1);
        for j = 1:length(face_list1)
            
            f1                 = face_list1(j,:);
            [elem_listA,dummy] = find(not(EToV - f1(1)));
            [elem_listB,dummy] = find(not(EToV - f1(2)));
            elem1      = intersect(elem_listA,elem_listB);
            
            lind1      = find(EToV(elem1,:)==f1(1));
            if(EToV(elem1,1+mod(lind1,Nfaces))==f1(2))
                lfind1  = lind1;
            else
                lfind1  = 1+mod(lind1+1,Nfaces);
            end
            
            f2         = face_list2(abs(pflink(j)),:);
            [elem_listA,dummy] = find(not(EToV - f2(1)));
            [elem_listB,dummy] = find(not(EToV - f2(2)));
            elem2      = intersect(elem_listA,elem_listB);
            
            lind2      = find(EToV(elem2,:)==f2(1));
            if(EToV(elem2,1+mod(lind2,Nfaces))==f2(2))
                lfind2  = lind2;
            else
                lfind2  = 1+mod(lind2+1,Nfaces);
            end
            
            
            EToE(elem1,lfind1) = elem2;
            EToE(elem2,lfind2) = elem1;
            EToF(elem1,lfind1) = lfind2;   
            EToF(elem2,lfind2) = lfind1;

            % Calculating perpendicular distance between periodic faces
            v1 = f1(1);
            if(pflink(j)>0)
                v2 = f2(1);
            else
                v2 = f2(2);
            end
            
            x1 = VX(v1); x2 = VX(v2); y1 = VY(v1); y2 = VY(v2);
            
            if(isKey(PShift,elem1))
                PShift(elem1) = [PShift(elem1);elem2,x2-x1,y2-y1];
            else
                PShift(elem1) = [elem2,x2-x1,y2-y1];
            end
            
            if(isKey(PShift,elem2))
                PShift(elem2) = [PShift(elem2);elem1,x1-x2,y1-y2];
            else
                PShift(elem2) = [elem1,x1-x2,y1-y2];
            end
                                                                                 
        end
        
    end    
end


% Setting boundary conditions for remaining boundary triangles/faces.
% Periodic boundary face tag is kept as 0
BFace_keys = keys(BFaces);
BCTag      = zeros(size(EToV));
for k = 1:length(BFace_keys)
    ckey   = BFace_keys{k};
    f_ind  = find(BC_flags(:,1)==ckey);
    if(~(BC_flags(f_ind,2) == Pflag))
       c_ind  = SetBC2D(BFaces(ckey),EToV); 
       BCTag  = BCTag + BC_flags(f_ind,1)*c_ind;
    end

end    

return;
