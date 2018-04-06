function badc_prob = ind_MLP1D(x, n_input,n_output,n_hidden_layer,leaky_alpha,...
                     WEIGHTS,BIASES)
   
   % CHECKS THAT SIZES ARE COMPATIBLE
   [m,n] = size(x);
   assert(m == n_input);
   assert(length(WEIGHTS) == n_hidden_layer+1);
   assert(length(BIASES) == n_hidden_layer+1);
   assert(m == size(WEIGHTS{1},2));
   
   for i=1:n_hidden_layer
       assert(size(WEIGHTS{i},1) == size(BIASES{i},1));
       assert(size(WEIGHTS{i},1) == size(WEIGHTS{i+1},2));
   end
   
   assert(size(WEIGHTS{end},1) == n_output);
   assert(size(BIASES{end},1) == n_output);
   
   x_in = x;
   
   %Scaling Data
   fact = repmat(max(max(abs(x_in),[],1),1),m,1);
   x_in = x_in./fact;
   % HIDDEN LAYERS
   for i=1:n_hidden_layer

       y = leaky_ReLU(WEIGHTS{i}*x_in + repmat(BIASES{i},1,n),leaky_alpha);
       clear x_in
       x_in = y;
       clear y
   end          
   
   %OUTPUL LAYER AND SOFTMAX
   y = softmax(WEIGHTS{end}*x_in + repmat(BIASES{end},1,n));
   %y = WEIGHTS{end}*x_in + BIASES{end};
   badc_prob = y(1,:);
                 
end                 