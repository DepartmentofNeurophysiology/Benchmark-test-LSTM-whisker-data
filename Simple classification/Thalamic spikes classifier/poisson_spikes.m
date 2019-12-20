
% Make poisson spikes for the 200 thalamic neurons

dt = 0.001;                             % 1 msec
rate = 50;                              % 50 spikes per second, on average
T = 1.000;                             % 1 sec simulation
T_vec = 0:dt:T;                         % a vector with each time step	


for n = 1:200
    vt = rand(size(T_vec));
    spikes = (rate*dt) > vt;
    thalamus_poisson(n).spike_times = find( spikes == 1);
end

figure(1)
for i = 1:200
    neuron = i*ones(length(thalamus_poisson(i).spike_times) ,1);
    plot(thalamus_poisson(i).spike_times, neuron,'k.')
    hold on  
end
hold off
%%
N = 50;                                    % Number of neurons in the reservoir.
Win = 0.1;                                   % Scale input weights.
Ein = Win*(rand(N, 200) - 0.0);            % These are the input weights connecting the thalamic neurons to the reservoir.
T = 1000;

[ neuron_input] = make_poisson_spikes_weighted(N, Ein, T, thalamus_poisson);

























