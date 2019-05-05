function [VX,VY,K,Nv,EToV,BFaces, PerBToB_map, PerBFToF_map] = read_gmsh_file(mshfile)

fprintf('... reading mesh from %s\n',mshfile)
fid = fopen(mshfile);

% REMOVE DUMMY LINES
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);

line = fgetl(fid);
assert(strcmp(line,'$Nodes'));

% Exctract number of vertices
Nv = str2double(fgetl(fid));
VX = zeros(1,Nv);
VY = zeros(1,Nv);

for i=1:Nv
    line    = fgetl(fid);
    C       = strsplit(line,' ');
    VX(1,i) = str2double(C{1,2});
    VY(1,i) = str2double(C{1,3});
end

line = fgetl(fid);
assert(strcmp(line,'$EndNodes'));

line = fgetl(fid);
assert(strcmp(line,'$Elements'));

% Get number of lines, faces etc
Nelem = str2double(fgetl(fid));

line = fgetl(fid);
C    = strsplit(line,' ');

nbface  = 0;

% BFaces is a map for each boundary face in the mesh file, where the key
% corresponds to the boundary tag, while the map value is an array
% containing the pair of vertices for faces with the same tag. Note that
% BFaces key uses the user defined tags.
% g2utag is a mapping between tag assigned by gmsh and user defined tag 
BFaces  = containers.Map('KeyType','uint32','ValueType','any');
g2utag  = containers.Map('KeyType','uint32','ValueType','uint32');
while(strcmp(C{1,2},'1'))
    nbface  = nbface+1;
    gftag    = str2double(C{1,end-2});
    uftag    = str2double(C{1,end-3});
    if(~isKey(g2utag,gftag))
        g2utag(gftag) = uftag;
    end
    if(isKey(BFaces,uftag))
        BFaces(uftag) = [BFaces(uftag);str2double(C{1,end-1}),str2double(C{1,end})];
    else
        BFaces(uftag) = [str2double(C{1,end-1}),str2double(C{1,end})];
    end
    
    line = fgetl(fid);
    C    = strsplit(line,' ');
end

K    = Nelem - nbface; 
EToV = zeros(K,3);
for i=1:K-1
    EToV(i,:) = [str2double(C{1,end-2}),str2double(C{1,end-1}),str2double(C{1,end})];
    line      = fgetl(fid);
    C         = strsplit(line,' ');
end
EToV(K,:) = [str2double(C{1,end-2}),str2double(C{1,end-1}),str2double(C{1,end})];

line = fgetl(fid);
assert(strcmp(line,'$EndElements'));

line = fgetl(fid);

% PerBToB_map is a map where the key is the tag as the master domain boundary,
% and the value is the tag of the slave domain boundary.
% PerBFToF_map is a map where the key is the tag as the master domain boundary,
% and the value is an array containing the indices of the perioidic face
% in the slave domain boundary. Using the index and the tag of the slave 
% boundary form PerBToB_map, the vertices of the periodic face can be
% extracted from BFaces. Furthermore, if an index is negative, this means
% that the orientation of the periodic slave face needs to be flipped.
PerBToB_map  = containers.Map('KeyType','uint32','ValueType','uint32');
PerBFToF_map = containers.Map('KeyType','uint32','ValueType','any');

Per_avail     = false;
while ischar(line)
    if(strcmp(line,'$Periodic'))
        n_per_elem = str2double(fgetl(fid));
        Per_avail  = true;
        for i = 1:n_per_elem
            line  = fgetl(fid);
            C     = strsplit(line,' ');
            n_ent = str2double(fgetl(fid));
            if(strcmp(C{1,1},'1'))
                bf1   = str2double(C{1,2}); bf1 = g2utag(bf1);
                bf2   = str2double(C{1,3}); bf2 = g2utag(bf2);
                nface = length(BFaces(bf1));
                assert(nface == length(BFaces(bf2)))
                PerBToB_map(bf1) = bf2;
                f2f       = zeros(nface,nface);
                bf1_vlist = BFaces(bf1);
                bf2_vlist = BFaces(bf2);
                p2p = containers.Map('KeyType','uint32','ValueType','uint32');
                for j = 1:n_ent
                    line  = fgetl(fid);
                    nodes = strsplit(line,' ');
                    p1    = str2double(nodes{1,1});
                    p2    = str2double(nodes{1,2});
                    p2p(p1) = p2;
                    ind1 = find((bf1_vlist(:,1) == p1)|(bf1_vlist(:,2) == p1));
                    ind2 = find((bf2_vlist(:,1) == p2)|(bf2_vlist(:,2) == p2));
                    for i1=ind1
                        for i2=ind2
                            f2f(i1,i2) = f2f(i1,i2) + 1;
                        end
                    end
                end
                f2flist = zeros(nface,1);
                for f1 = 1:nface
                    f2 = find(f2f(f1,:)==2);
                    if(p2p(bf1_vlist(f1,1)) ~= bf1_vlist(f2,1)) 
                        f2 = -1*f2; % Negative integer means periodic face has opposite orienatation
                    end
                    f2flist(f1) = f2;
                end
                PerBFToF_map(bf1) = f2flist;
            else
                for j = 1:n_ent
                    line  = fgetl(fid);
                end
            end
                    
        end 
    end
    line = fgetl(fid);
end

fclose(fid);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

% Deisplay mesh data
fprintf('    -> Number of cells    = %d\n',K)
fprintf('    -> Number of vertices = %d\n',Nv)
if(Per_avail)
    fprintf('    -> Periodic faces available\n')
else
    fprintf('    -> No periodic faces available\n')
end

return
