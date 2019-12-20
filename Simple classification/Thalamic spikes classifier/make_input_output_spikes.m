function [ z_all, neuron_input, curve_trace ] = make_input_output_spikes( N, n_t, n_r, SpikeTrainStruct,trials,dat, Ein, pole, rate, WhiskerTrace, random)
% This function takes a number of trials as input and creates an input
% signal of n trials long where the inputs of 200 thalamus neurons are randomly arranged.
% The corresponding target output is also returned as is the target label.
% rate is the average firing rate of the poisson spikes inbetween the
% trials.
for n = 1:n_t
    spike_array{n}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,n}));
end
for n = 1:n_t
    for i = 1:200
        spike_array{n}.trial(i,SpikeTrainStruct{1, 1}.SpikeTimes{i, n}) = ones(length(SpikeTrainStruct{1, 1}.SpikeTimes{i, n}) ,1);
    end
end
for n = 1:n_t
    spike_array{1,n}.neuron_input = zeros(N,length(spike_array{1,n}.trial));
    for t=1:length(spike_array{n}.trial)
        index_input = find(spike_array{1,n}.trial(:,t) == 1);  % Find thalamic input neurons that have spiked
        if length(index_input) ~= 0
            spike_array{1,n}.neuron_input(:,t) = sum(Ein(:,index_input),2);
        end
    end
end
% spike_array{1}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,1}));
% spike_array{2}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,2})); 
% spike_array{3}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,3})); 
% spike_array{4}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,4})); 
% for t = 1:4
%     for i = 1:200
%         spike_array{t}.trial(i,SpikeTrainStruct{1, 1}.SpikeTimes{i, t}) = ones(length(SpikeTrainStruct{1, 1}.SpikeTimes{i, t}) ,1);
%     end
% end
% for n = 1:4
%     spike_array{1,n}.neuron_input = zeros(N,length(spike_array{1,n}.trial));
%     for t=1:length(spike_array{n}.trial)
%         index_input = find(spike_array{1,n}.trial(:,t) == 1);  % Find thalamic input neurons that have spiked
%         if length(index_input) ~= 0
%             spike_array{1,n}.neuron_input(:,t) = sum(Ein(:,index_input),2);
%         end
%     end
% end

% Some parameters for the pulse
reset = 800; % Time inbetween the trials
pulse_length = 400; % Length of pulse
amp = 3; % Amplitude of pulse
decay = 80; % decay of exponential pulse
constant = 200; % How long the pulse is kept constant
neuron_input = [];
z_all = [];
curve_trace = [];
start_early = 500; % How long to start before the end of the input
%target_label = zeros(1000,1)';

% Poisson spikes for inbetween trials 
dt = 0.001;                                     % 1 msec                          
T_vec = 0:dt:reset*dt;                         % a vector with each time step	

if random == 1
if pole == 1
    for t = 1:n_r
        trial = randi([1,length(trials)],1,1);
        l_t = length(spike_array{1,trial}.neuron_input(1,:)) + reset;
        z_t = zeros(l_t,1);
        z_t(length(spike_array{1,trial}.neuron_input(1,:))-start_early:1:length(spike_array{1,trial}.neuron_input(1,:))+constant) = -dat(trials(trial)).pole*amp*ones(start_early+constant+1,1);
        z_t(length(spike_array{1,trial}.neuron_input(1,:))+constant:1:length(spike_array{1,trial}.neuron_input(1,:))+pulse_length+constant) =  -dat(trials(trial)).pole*amp*exp(-(0:1:pulse_length)./decay);
        %target_label = [target_label dat(trial).pole*ones(length(X(trial).cnz) +pulse_length+constant,1)' zeros(reset -pulse_length-constant,1)' ];
        z_all = [z_all z_t'];
        
        % make the poisson input 
        for n = 1:200
            vt = rand(size(T_vec) - [0 1]);
            spikes = (rate*dt) > vt;
            thalamus_poisson(n).spike_times = find( spikes == 1);
        end
        [ poisson_input] = make_poisson_spikes_weighted(N, Ein, reset, thalamus_poisson);

        neuron_input = [neuron_input spike_array{1,trial}.neuron_input poisson_input];  
        curve_trace = [curve_trace WhiskerTrace.Recording{2,trial} zeros(reset,1)'];
    end
else
    for t = 1:n_t
        trial = randi([1,length(trials)],1,1);
        l_t = length(spike_array{1,trial}.neuron_input(1,:)) + reset;
        z_t = zeros(l_t,1);
        z_t(dat(trials(trial)).pole_times(2)-start_early:1:dat(trials(trial)).pole_times(2)+constant) = -dat(trials(trial)).pole*amp*ones(start_early+constant+1,1);
        z_t(dat(trials(trial)).pole_times(2)+constant:1:dat(trials(trial)).pole_times(2)+pulse_length+constant) =  -dat(trials(trial)).pole*amp*exp(-(0:1:pulse_length)./decay);
        %target_label = [target_label dat(trial).pole*ones(length(X(trial).cnz) +pulse_length+constant,1)' zeros(reset -pulse_length-constant,1)' ];
        z_all = [z_all z_t'];
        neuron_input = [neuron_input spike_array{1,trial}.neuron_input zeros(reset,N)'];  
    end
end

else
    for t = 1:n_t
        trial = t;
        l_t = length(spike_array{1,trial}.neuron_input(1,:)) + reset;
        z_t = zeros(l_t,1);
        z_t(length(spike_array{1,trial}.neuron_input(1,:))-start_early:1:length(spike_array{1,trial}.neuron_input(1,:))+constant) = -dat(trials(trial)).pole*amp*ones(start_early+constant+1,1);
        z_t(length(spike_array{1,trial}.neuron_input(1,:))+constant:1:length(spike_array{1,trial}.neuron_input(1,:))+pulse_length+constant) =  -dat(trials(trial)).pole*amp*exp(-(0:1:pulse_length)./decay);
        %target_label = [target_label dat(trial).pole*ones(length(X(trial).cnz) +pulse_length+constant,1)' zeros(reset -pulse_length-constant,1)' ];
        z_all = [z_all z_t'];
        
        % make the poisson input 
        for n = 1:200
            vt = rand(size(T_vec) - [0 1]);
            spikes = (rate*dt) > vt;
            thalamus_poisson(n).spike_times = find( spikes == 1);
        end
        [ poisson_input] = make_poisson_spikes_weighted(N, Ein, reset, thalamus_poisson);

        neuron_input = [neuron_input spike_array{1,trial}.neuron_input poisson_input];  
        curve_trace = [curve_trace WhiskerTrace.Recording{2,trial} zeros(reset,1)'];
    end
    
end
end


