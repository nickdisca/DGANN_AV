function [n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES] = ...
    read_mlp_param1D(nn_model)

Model_dir = horzcat('../Trained_networks/1D/',nn_model);

param_file = horzcat(Model_dir,'/model_parameters.dat');

fid = fopen(param_file);
PAR = textscan(fid,'%s %f','delimiter',' ');
fclose(fid);

assert(strcmp(PAR{1}{1},'IDIM'));
assert(strcmp(PAR{1}{2},'ODIM'));
assert(strcmp(PAR{1}{3},'NHL'));
assert(strcmp(PAR{1}{4},'LEAKY_P'));

n_input        = round(PAR{2}(1));
n_output       = round(PAR{2}(2));
n_hidden_layer = round(PAR{2}(3));
leaky_alpha    = PAR{2}(4);

WEIGHTS = cell(n_hidden_layer + 1,1);
BIASES  = cell(n_hidden_layer + 1,1);

for i=1:n_hidden_layer
    WEIGHTS{i,1} = load(horzcat(Model_dir,'/w_h',num2str(i),'.dat'))';
    BIASES{i,1}  = load(horzcat(Model_dir,'/b_h',num2str(i),'.dat'));
end

WEIGHTS{end,1} = load(horzcat(Model_dir,'/w_out.dat'))';
BIASES{end,1}  = load(horzcat(Model_dir,'/b_out.dat'));

end