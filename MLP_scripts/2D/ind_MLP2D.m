function badc_prob = ind_MLP2D(x, Net)
   
   % CHECKS THAT SIZES ARE COMPATIBLE
   [m,n] = size(x);
   assert(m == Net.n_input);
   assert(length(Net.WEIGHTS) == Net.n_hidden_layer+1);
   assert(length(Net.BIASES) == Net.n_hidden_layer+1);
   assert(m == size(Net.WEIGHTS{1},2));
   
   for i=1:Net.n_hidden_layer
       assert(size(Net.WEIGHTS{i},1) == size(Net.BIASES{i},1));
       assert(size(Net.WEIGHTS{i},1) == size(Net.WEIGHTS{i+1},2));
   end
   
   assert(size(Net.WEIGHTS{end},1) == Net.n_output);
   assert(size(Net.BIASES{end},1) == Net.n_output);
   
   % SCALE DATA BEFORE USING IN NETWORK
   x_in = Scaling(x);
   
   % HIDDEN LAYERS
   for i=1:Net.n_hidden_layer

       y = leaky_ReLU(Net.WEIGHTS{i}*x_in + repmat(Net.BIASES{i},1,n),Net.leaky_alpha);
       clear x_in
       x_in = y;
       clear y
   end          
   
   %OUTPUL LAYER AND SOFTMAX
   y = softmax(Net.WEIGHTS{end}*x_in + repmat(Net.BIASES{end},1,n));
   %y = WEIGHTS{end}*x_in + BIASES{end};
   badc_prob = y(1,:);
                 
end                 