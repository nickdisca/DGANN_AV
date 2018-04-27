% This function is used to find th relative position of a script with
% respect to the main DGANN folder containing mypath.m

REL_PATH  = '';
path_file = sprintf('%smypath.m',REL_PATH);
found     = exist(path_file,'file');

while(found == 0)
    REL_PATH  = sprintf('../%s',REL_PATH);
    path_file = sprintf('%smypath.m',REL_PATH);
    found     = exist(path_file,'file');
end