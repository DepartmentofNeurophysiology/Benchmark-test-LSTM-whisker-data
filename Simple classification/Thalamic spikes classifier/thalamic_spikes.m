addpath(genpath('Cortical-representation-of-touch-in-silico'))
addpath(genpath('Data'))
addpath(genpath('Helper functions'))
load('whiskingstruct_an197522_2013_03_08')
load('Kernels')
load('spikes_5_8_15_18')

%%
plot_trial(18,dat)
%%
trials = [5 8 15 18];
nt = 4;

% figure(1)
% subplot(3,1,1)
% plot(WhiskerTrace.Recording{1,nt})
% title('angle')
% %ylim([-3 3])
% subplot(3,1,2)
% plot(WhiskerTrace.Recording{2,nt})
% %ylim([-8 6])
% title('curvature')
% subplot(3,1,3)
for i = 1:200
    neuron = i*ones(SpikeTrainStruct{1, 1}.SpikeCount{i,nt} ,1);
    plot(SpikeTrainStruct{1,1}.SpikeTimes{i,nt}, neuron,'k.')
    hold on  
end
hold off
%%
% figure(2)
% nt = 4;
% subplot(3,1,1)
% plot(WhiskerTrace.Recording{1,nt})
% title('angle')
% %ylim([-3 3])
% subplot(3,1,2)
% plot(WhiskerTrace.Recording{2,nt})
% %ylim([-8 6])
% title('curvature')
% subplot(3,1,3)
for i = 1:200
    neuron = i*ones(SpikeTrainStruct{1, 1}.SpikeCount{i,nt} ,1);
    plot(SpikeTrainStruct{1,1}.SpikeTimes{i,nt}, neuron,'k.')
    hold on  
end
hold off

%%

spike_array{1}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,1}));
spike_array{2}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,2})); 
spike_array{3}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,3})); 
spike_array{4}.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,4})); 
for t = 1:4
    for i = 1:200
        spike_array{t}.trial(i,SpikeTrainStruct{1, 1}.SpikeTimes{i, t}) = ones(length(SpikeTrainStruct{1, 1}.SpikeTimes{i, t}) ,1);
    end
end
%% 
N = 10;
Ein = rand(N, 200); % these are the input weights connecting the thalamic neurons to the reservoir.

for n = 1:4
    spike_array{1,n}.neuron_input = zeros(N,length(spike_array{1,n}.trial));
    for t=1:length(spike_array{n}.trial)
        index_input = find(spike_array{1,n}.trial(:,t) == 1);  % Find thalamic input neurons that have spiked
        if length(index_input) ~= 0
            spike_array{1,n}.neuron_input(:,t) = sum(Ein(:,index_input),2);
        end
    end
end





