clc
clear all
close all


current_model = 'OP1';	
data_set      = 'I';
data_subset   = '1';

x = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15]

[n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES] = ...
    read_mlp_param(current_model,data_set,data_subset);

ind = ind_MLP(x,n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES)

