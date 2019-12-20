addpath(genpath('../Cortical-representation-of-touch-in-silico'))
addpath(genpath('../Data'))
addpath(genpath('../Helper functions'))
load('whiskingstruct_an197522_2013_03_08')
load('Kernels')
%load('Spike_train_5_8_15_18')
%load('WhiskerTrace')
%%
Nbx = 1;                                                % # of barrels 'x-direction'
Nby = 1;                                                % # of barrels 'y-direction'
barrelstruct = cell(Nbx, Nby);
for nbx = 1:Nbx
    for nby = 1:Nby
        barrelstruct{nbx, nby}.xpos         = (nbx-2)*300;      % barrel position
        barrelstruct{nbx, nby}.ypos         = (nby-1)*300;
        barrelstruct{nbx, nby}.Nthalamic    = 200;              % # filter neurons for this barrel
        barrelstruct{nbx, nby}.mainbarrel   = 3;
    end
end
    % Choose main and secondary barrels
% barrelstruct{2,1}.mainbarrel    = 1; % main
% barrelstruct{1,1}.mainbarrel    = 2; % secondary
% barrelstruct{3,1}.mainbarrel    = 2; % tertiary
barrelstruct{1,1}.mainbarrel    = 1; % main
%%
Nbarrel = Nbx*Nby;
SpikeTrainStruct = cell(Nbx,Nby);
SpikeGenStruct = cell(Nbx,Nby);
nb = 0;
pole = 1; % Choose whether to make whisker trace 
% easy trials training set = [5 8 15 18 24 76 28 29 33 34 82 93]
% hard trials training set = [5 8 15 18 24 76 28 29 33 34 82 93 30 50 71 83]

% easy trials testing set = [2 7 19 25]
% hard trials testing set = [64 73 84 87 91 100]
trials = [64 73 84 87 91 100];
for i = 1:length(trials)
    [WhiskerTrace.Recording{2,i} ,WhiskerTrace.Recording{1,i}]  = make_whisker_trace( trials(i), dat, pole);
end

save('WhiskerTrace_hard_test','WhiskerTrace')
%%
WhiskerTrace.binsize = 1;
Barrelstruct = barrelstruct;
for nbx = 1:Nbx
    for nby = 1:Nby
        nb = nb+1;
        disp(['Making Thalamic spike trains for barrel ' num2str(nb) '/' num2str(Nbarrel)])
        SpikeGenStruct{nbx,nby}.refra             = 3; % refractory period (ms)
        SpikeGenStruct{nbx,nby}.Ntrial_pertrace   = 1; % # trials for each Deflection trace
        SpikeGenStruct{nbx,nby}.binsize           = 1; % binsize spike trains (ms)
        if Barrelstruct{nbx,nby}.mainbarrel == 1
            
            SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
            SpikeGenStruct{nbx,nby}.scaling           = 1; % Scaling of PSTH
        elseif Barrelstruct{nbx,nby}.mainbarrel == 2
            
            SpikeGenStruct{nbx,nby}.delay             = 2.5; % (ms)
            SpikeGenStruct{nbx,nby}.scaling           = .3;
        elseif Barrelstruct{nbx,nby}.mainbarrel == 3
            
            SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
            SpikeGenStruct{nbx,nby}.scaling           = 0;
        end
        
        plotyn = 1;
        SpikeTrainStruct{nbx,nby} = kernel_recording_to_spiketrain(WhiskerTrace, KernelStruct{nbx,nby}, SpikeGenStruct{nbx,nby}, [], plotyn);
    end
end


%%
save('spikes_hard_test','SpikeTrainStruct')


