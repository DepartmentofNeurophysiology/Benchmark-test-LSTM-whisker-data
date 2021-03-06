function [ connections ] = thalamus_cone_connection( N_x, N_y, N_thx, N_thy, K_out, sigmaffwd )
% This function makes connections for each thalamus neurons to K_out
% reservoir neurons within a range of sigmaffwd. This is done by projecting
% the 1D thalamic array onto a 2D rectangle of size N_th x N_thy. For each
% neuron there is then added a distance generated by a gaussian
% distribution. The closest neuron in the 2D rectangle of reservoir neurons
% is then taken to make a connection with. This is done K_out times and
% thus there each thalamic neuron makes K_out connections with K_out
% reservoir neurons. 

% Define number of neurons along each dimension for thalamic feedforward
% and for reservoir network.

% N_x = 6;
% N_y = 3;
% N_thx = 2;
% N_thy = 1;

N_th = N_thy*N_thx;
N = N_x*N_y;

% Define number of connections each thalamus neurons makes with the
% reservoir neurons.
%
% K_out = 2;
%
% Feedforward connection width
% sigmaffwd = 0.1;


% Define array in which to store the connections each thalamic neuron makes
% with the reservoir. Each thalamic neuron makes K_out connections so this
% is stored in a N_th x K_out array. 

connections = zeros(N_th, K_out);

for n = 1:N_th
    
    % Transform 1D index to 2D index
    
    j = ceil(n / N_thx);
    i = n - (j -1)* N_thx;
    
    % For each neuron generate outgoing connections
    
    for k = 1:K_out
        
        % Generate pair of numbers from Gausian distribution to represent the
        % distance in each direction from the presynaptic neuron to the
        % postsynaptic neuron.
        
        z = normrnd(0, sigmaffwd, [1,2]);
        
        % Generate coordinates of presynaptic neurons in square [0,1]
        
        y = [i / N_thx , j / N_thy];
        
        % Generate target coordinates in square [0,1]
        
        x = mod(z + y, 1);
        
        % Find index closest to target coordinates
        
        j1 = ceil(N_x * x(1));
        j2 = ceil(N_y * x(2));
        
        % Transform 2D index to 1D index
        connections(n, k)  = (j1-1)*N_y + j2;
    end
end

end





