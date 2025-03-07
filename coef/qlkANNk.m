%Takes as input the physical R/LTi, q, s, Ti/Te temperature ratio, and the neural net weights and biases

%Input can be vectors, so may points can be calculated at the same time
%Output is the ion heat flux vector

%%% NETWORK VALIDITY RANGE: Ati [2 12]; ti/te [0.3 3]; q [1 5]; s [0.1 3]

% rlti input R/LTi , q input q-profile, shear input magnetic shear, max_qlk/min_qlk the output normalizations
% net the neural network, prepos sets output normalization procedure depending on the preprocessing of output

function NNout = NNoutput(rlti,q,tite,shear,max_qlk,min_qlk,net,prepros)

%initialize input matrix in order expected by the trained network. 
invec=[q,shear,tite,rlti]';

%find size of input vector (if rlti, q, tite, shear are not same length, the function will crash)
npoints = numel(rlti);

%Normalize input matrix to match NN input normalization
scalefac =[2 ; 1.45 ; 1.35 ; 5] ; 
dispfac  =[3 ; 1.55 ; 1.65 ; 7] ;

%normalize input vector to match input normalization in NN
unorm = bsxfun(@times,bsxfun(@minus,invec,dispfac),1./scalefac);

%calculate regression neural network piece by piece, for all radii
g=bsxfun(@plus,net.IW*unorm,net.b1);
sg = sigmo(g);
f=bsxfun(@plus ,net.L1W*sg,net.b2);
sf = sigmo(f);

output=(net.L2W*sf + net.b3*ones(1,length(rlti)))';

% if outputs processing between 0 and 1%
if prepros==0
 NNout= output.*(ones(length(rlti),1)*max_qlk'-ones(length(rlti),1)*min_qlk')+ones(length(rlti),1)*min_qlk';
end 

% if outputs processing between -1 a +1
if prepros==-1
 NNout= (output+ones(size(output))).*(ones(length(rlti),1)*max_qlk'-ones(length(rlti),1)*min_qlk')/2+ones(length(rlti),1)*min_qlk';
end

%nonlinear neuron function. 
function s = sigmo(x)
s= 2./(1+exp(-2.*x))-1; %define sigmoid function (nonlinear transfer functions in NN)
return

