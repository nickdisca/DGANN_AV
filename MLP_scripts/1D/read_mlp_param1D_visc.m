function Net = read_mlp_param1D_visc(nn_model,REL_PATH,N)

fprintf('... loading viscosity network weights\n')

Net.avail    = true;
Net.nn_model = nn_model;
Net.NN_Dir   = horzcat(REL_PATH,'Trained_networks/1D/',nn_model);

% Add NN_dir to path. This will be removed at the end of the solver run
addpath(Net.NN_Dir);

param_file = horzcat(Net.NN_Dir,'/model_parameters.dat');

fid = fopen(param_file);
PAR = textscan(fid,'%s %f','delimiter',' ');
fclose(fid);

assert(strcmp(PAR{1}(1),'NHL'));
assert(strcmp(PAR{1}(2),'LEAKY_P'));

Net.n_input        = N+1;
Net.n_output       = N+1;
Net.n_hidden_layer = round(PAR{2}(1));
Net.leaky_alpha    = PAR{2}(2);

Net.WEIGHTS = cell(Net.n_hidden_layer + 1,1);
Net.BIASES  = cell(Net.n_hidden_layer + 1,1);

for i=1:Net.n_hidden_layer
    Net.WEIGHTS{i,1} = load(horzcat(Net.NN_Dir,'/w_h',num2str(i),'_m',num2str(N),'.dat'))';
    Net.BIASES{i,1}  = load(horzcat(Net.NN_Dir,'/b_h',num2str(i),'_m',num2str(N),'.dat'));
end

Net.WEIGHTS{end,1} = load(horzcat(Net.NN_Dir,'/w_out_m',num2str(N),'.dat'))';
Net.BIASES{end,1}  = load(horzcat(Net.NN_Dir,'/b_out_m',num2str(N),'.dat'));

end