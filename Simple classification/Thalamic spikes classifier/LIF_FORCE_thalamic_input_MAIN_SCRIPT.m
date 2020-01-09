%%  Leaky integrate and fire FORCE neural network training with thalamic spike trains as input.
% This script takes a series of whisker traces (in this case trial 5, 8,
% 15, 18) as input. Trial 5 and 18 have the same positive output and 
% trial 8 and 15 the same negative output. These corrsepond to two different
% pole locations. The thalamic spike trains of these trials are taken as
% input to the network reservoir. There are 200 thalamic neurons and N
% neurons in the reservoir. The 200 input thalamic neurons are weighted by
% a 200xN input matrix. This weights matrix is created randomly. It can be
% chosen what kind of connections from thalamus to reservoir to create.
% First an all to all connection matrix was generated with random weights
% for each connection. This however seemed rather unrealistic. So then 
% connections from thalamus to reservoir were made based on a 'projection
% cone'. This basically means that all of the thalamus neurons and
% reservoir neurons are projected onto a circular disk. For each thalamus
% neuron K_out reserovir neurons were connected with according to a
% gaussian distribution. These connections were then weighted with a random
% or not random weight.
%
% To make this more computationally efficient the
% spikes of the input are stored per ms and are multiplied with the input
% weights. This create a 2000x1 input array for each ms of the trial.
%
% This input is then fed into the network. This is done by filtering te spikes through
% a double exponential filter. 
%
% The network is then trained through RLS to give the correct output for a
% specific input trial, and thus is trained to distinguish different
% thalamic spike train inputs. 
%
% The trained network is then tested by turning the RLS off and checking to
% see if the network can classify the input correctly on its own. 

%% First load all the whisker traces and the ( earlier made) thalamic spike trains.
addpath(genpath('..\Cortical-representation-of-touch-in-silico'))
addpath(genpath('..\Data'))
addpath(genpath('..\Helper functions'))

load('whiskingstruct_an197522_2013_03_08')                                      % Structure containing the whisker traces for each trial. 
WhiskerTrace_pole = load('WhiskerTrace_5_8_15_18_24_76_28_29_33_34_82_93');     % Structure containing the whisker traces for pole in reach times of trial 5,8,15, and 18.
SpikeTrainStruct_pole = load('spikes_5_8_15_18_24_76_28_29_33_34_82_93');       % Structure containing the thalamic spike trains for trial 5, 8, 15, and 18.
SpikeTrainStruct_whole = load('spikes_5_8_15_18_whole_trial');                  % Structure containing thalamic spikes for the whole trial.
WhiskerTrace_whole = load('WhiskerTrace_5_8_15_18_whole_trial');                % Structure containing whisker traces for whole trial.
SpikeTrain_hard_train = load('spikes_hard');                                    % Also contains 4 'hard' trials
WhiskerTrace_hard_train = load('WhiskerTrace_hard');                            % Also contains 4 'hard' trials

WhiskerTrace_test = load('WhiskerTrace_2_7_19_25_pole_trial');                  % Structure containing the whisker traces to test on.
SpikeTrainStruct_test = load('spikes_2_7_19_25_pole_trial');                    % Structure containing the thalamic spike trains to test on.
SpikeTrain_hard_test = load('spikes_hard_test');                                % 'Hard' trials to test on
WhiskerTrace_hard_test = load('WhiskerTrace_hard_test');                        % 'Hard' trials to test on

%% Plot the curvature and angle trace to compare to the thalamic spike trains.
trials_train = [5 8 15 18 24 76 28 29 33 34 82 93];
trials_test = [2 7 19 25];
trials = trials_train;
% 
% trials_train_hard = [5 8 15 18 24 76 28 29 33 34 82 93 30 50 71 83];
% trials_test_hard = [64 73 84 87 91 100];
% trials = trials_train_hard; 
% WhiskerTrace = WhiskerTrace_hard_train.WhiskerTrace;                                  % This is the training set
WhiskerTrace = WhiskerTrace_pole.WhiskerTrace;      % testing set
% SpikeTrainStruct = SpikeTrain_hard_train.SpikeTrainStruct;
SpikeTrainStruct = SpikeTrainStruct_pole.SpikeTrainStruct;

nt = 5;      % Select which trial to plot
figure(8)
subplot(3,1,1)
plot(WhiskerTrace.Recording{1,nt})
ylim([-1 0.0])
title('Angle')
%ylim([-3 3])
subplot(3,1,2)
plot(WhiskerTrace.Recording{2,nt})
ylim([-0.01 0.01])
title('Curvature')
subplot(3,1,3)
for i = 1:200
    neuron = i*ones(SpikeTrainStruct{1, 1}.SpikeCount{i,nt} ,1);
    plot(SpikeTrainStruct{1,1}.SpikeTimes{i,nt}, neuron,'k.')
    hold on  
end
xlabel('Time /ms')
title('Thalamic spikes')
hold off
%% Average firing rate per trial 
% The firing rate for trial 8 and 15 are quite low. This is useful to
% determine the poisson firing rate inbetween trials.
for nt = 1:length(trials)
    rate = (sum(cell2mat(SpikeTrainStruct{1,1}.SpikeCount(:,nt))))/(200*length(WhiskerTrace.Recording{1,nt})/1000);
    fprintf('rate_%d = %f \n', trials(nt), rate)
end
%% 1. Create random 'all to all' input weights matrix
N = 2000;                                   % Number of neurons in the reservoir.
Win = 10;                                   % Scale input weights.
Ein = Win*(rand(N, 200) - 0.0);             % These are the input weights connecting the thalamic neurons to the reservoir.

%% 2. Or create cone connections with corresponding input weights matrix
N_th = 200;                                 % This is the number of thalamus neurons
N = 1600;                                   % Number of neurons in the reservoir.
sigmaffwd = 0.05;                           % Radius of the gaussian projection.
K_out = 50;                                 % Number of connections each thalamus neuron makes.
[ indices ] = disk_connections( N_th, N, sigmaffwd, K_out); % This function returns the reservoir neurons each thalamus neuron makes a connection with.
%%  Make the input weights matrix
Win = 30;                                    % Scale the input weights.
Ein = zeros(N,N_th);                        

for i = 1:N_th
    %Ein(indices(i,:), i) = Win * ones(K_out,1); %  Give all the neurons to connect with a weight.
    Ein(indices(i,:), i) = Win * rand(K_out,1); %  Give all neurons to connect with a random weight.
end

figure(2)
imagesc(Ein')
colorbar
%% This makes poisson input, but is not neccessary to be made for the simulation.
% Make poisson spikes for the 200 thalamic neurons
dt = 0.001;                             % 1 msec
rate = 20;                              % 50 spikes per second, on average
T = 1.000;                              % 1 sec simulation
T_vec = 0:dt:T;                         % a vector with each time step	
for n = 1:200
    vt = rand(size(T_vec)-[0 1]);
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
% Take the poisson spikes and multiply with the input weights matrix to
% make a N x T array which can be used as input for the neural network
N = 1600;
Win = 2;
Ein = Win*(rand(N, 200) - 0.0);    
T = 1000;
[ poisson_input] = make_poisson_spikes_weighted(N, Ein, T, thalamus_poisson);
%% Create the array of neural input to the reservoir and the target function the output has to match.
% The neural input is made by multiplying the spikes with the input
% weights. 
n_t = length(trials_train);                                   % Number of trials that are repeated after each other. So the input becomes 'nt' long.
pole = 1;                                  % This tells the funtion whether trials are just 'pole in reach' or the 'whole' trial.
rate = 1.5;                                 % This is the average firing rate of the poisson input spikes 
random = 0;                                 % Arrange trials randomly or not, if not random then n_t needs to equal the number of input trials.
n_r = 12;
[ z_all, neuron_input, curve_trace ] = make_input_output_spikes( N,n_t, n_r, SpikeTrainStruct,trials,dat, Ein, pole, rate, WhiskerTrace, random);

%% Plot to check
T = length(neuron_input);
subplot(3,1,1)
plot(z_all)
xlim( [0 T])
title('Target output')
subplot(3,1,2)
plot(curve_trace,'LineWidth',1.5)
title('Curvature')
xlim( [0 T])
subplot(3,1,3)
imagesc(neuron_input)
xlim( [0 T])
%ylim([1 5])
colorbar('southoutside')
ylabel('Neuron input')
xlabel('Time /ms')
title('Weighted spikes')

%% Setup network and give it the thalamic spikes as input but no RLS. Then look at the reservoir spikes.
% Set last three terms to 0 so that the random weights are taken. Check if
% you want Dale's law or not. BIAS controls the firing rate of the neurons
% (-40 is the threshold). 

% G scales the static weights and Q scales the feeback weights.
 G = 1; Q = 1; alpha = 0.25; dt = 0.05; td = 10; dale = 0; INPUT = 1;BIAS = -41; Bphi = 0; Omega = 0; Efb = 0; wij = 0;
tic
[ E, OMEGA, BPhi, tspike_train,current ,ns,v_rec] = LIF_network_spikes( neuron_input, z_all, N, G, Q, alpha, dt, td, dale,INPUT,BIAS,Bphi,Omega,Efb,wij );
toc

%% Plot the spikes of the reservoir

T = length(neuron_input(1,:));
tspike = tspike_train;
x_lim= [0, T];
y_lim= [0, N];

figure(1)
subplot(4,1,1)
imagesc(neuron_input)
ylim([0,200])
xlim(x_lim)
ylabel('Neuron')

subplot(4,1,2)
plot(tspike(:,2),tspike(:,1),'k.')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Network spikes')
subplot(4,1,3)
plot(v_rec)
xlim(x_lim)
ylim([-70 -35])
ylabel('Voltage')
title('Ten random neurons')
 subplot(4,1,4)
 plot(z_all)
% % for i = 1:200
% %     neuron = i*ones(length(thalamus_poisson(i).spike_times) ,1);
% %     plot(thalamus_poisson(i).spike_times, neuron,'k.')
% %     hold on  
% % end
% % hold off
% xlim(x_lim)
% % title('Thalamic spikes')
 title('Target')
 subplot(5,1,5)
[A , edges] = histcounts(tspike(:,2),1:1:length(neuron_input(1,:)));
A(1) = 0;

plot(A)
xlim(x_lim)
xlabel('Time /ms')
ylabel('Av. pop. act.')
suptitle('Spikes with input but untrained')
%% Setup network and FORCE train so the output matches the target
% Save the weights E, OMEGA, and BPhi so we can test the classification
% performance of the network in the next section.

%zx = z_all; G = 10; Q = 30; alpha = 0.25; dt = 0.05; td = 20; 
G = 1; Q = 1; alpha = 0.25; dale = 0; INPUT = 1;BIAS = -41; step = 20; dt = 0.05; td = 10; rep = 1;
Bphi = 0; Omega = 0; Efb = 0; % To train again with the already trained weights or not
tic
[ E, OMEGA, BPhi, tspike_train,current, wij, BPhi_rec, delta_wij_3,r_rec ] = ...
LIF_network_spikes_train( neuron_input, z_all, N, G, Q, alpha, dt, td, dale,INPUT,BIAS,step,rep,Bphi,Omega,Efb);
toc
%%

%%
figure(2)
T = length(z_all);
i = T/dt;
% plot(dt*(1:1:i),z_plot(1:1:i),'g--','LineWidth',2), hold on
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r.','LineWidth',0.8), hold off
ylim([ -4 4])
xlim([ 0 T])
title('Training')
%%
plot(BPhi_rec)
xlim([0 T])
%ylim([-65 -35])
ylabel('BPhi')
title('Ten random output weights')
%%
tspike = tspike_train;
x_lim= [0,T];
y_lim= [0,1000];

figure(3)
subplot(4,1,1)
imagesc(neuron_input)
ylim([0,200])
xlim(x_lim)
ylabel('Neuron')
title('Input')
subplot(4,1,2)
plot(tspike(:,2),tspike(:,1),'k.')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Reservoir spikes')
subplot(4,1,3)
plot(z_all)
xlim(x_lim)
title('Target')
subplot(4,1,4)
[A , edges] = histcounts(tspike(:,2),1:1:length(neuron_input(1,:)));
A(1) = 0;
plot(A)
xlim(x_lim)
xlabel('Time /ms')
ylabel('pop. activity')
title('Average population activity')

%% Test trained network with learnt weights. Input the learned weights into the function so the network can classify the thalamic spike input.
% select test data, these should be 'easy to classify' trials.
n_t = 4;                                   % Number of trials.
pole = 1;                                  % This tells the funtion whether trials are just 'pole in reach' or the 'whole' trial.
rate = 1.5;                                 % This is the average firing rate of the poisson input spikes 
n_r = 10;                                    % Number of trials that are repeated after each other. So the input becomes 'nt' long.

SpikeTrainStruct = SpikeTrainStruct_test.SpikeTrainStruct;
WhiskerTrace = WhiskerTrace_test.WhiskerTrace;
% SpikeTrainStruct = SpikeTrain_hard_test.SpikeTrainStruct;
% WhiskerTrace = WhiskerTrace_hard_test.WhiskerTrace;
random = 0;
trials = [2 7 19 25];
% trials = trials_test_hard; 
[ z_all, neuron_input, curve_trace ] = make_input_output_spikes( N,n_t, n_r, SpikeTrainStruct,trials,dat, Ein, pole, rate, WhiskerTrace, random);

%
% Plot to check
subplot(3,1,1)
plot(z_all)
title('Target output')
subplot(3,1,2)
plot(curve_trace)
title('Curvature')
subplot(3,1,3)
imagesc(neuron_input)
%ylim([1 5])

colorbar('southoutside')
ylabel('Neuron input')
xlabel('Time /ms')
title('Weighted spikes')


%%
Bphi = BPhi; Omega = OMEGA; Efb = E; BIAS = -41; 
% The wij term is to make sure that the feedback weights interact with the
% reservoir in a way that respect dales law.
tic
[ E, OMEGA, BPhi, tspike_train,current,ns ] = LIF_network_spikes( neuron_input, z_all, N, G, Q, alpha, dt, td, dale,INPUT,BIAS,Bphi,Omega,Efb,wij );
toc
%% Plot to check result

figure(4)

T = length(z_all);
i = T/dt;
subplot(2,1,1)
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
subplot(2,1,2)
plot(curve_trace,'k', 'linewidth',1.5)
xlim([ 0 T])
ylabel('Curvature')
xlabel('Time /ms')
title('Whisker trace')

%% Plot spikes

tspike = tspike_train;
x_lim= [2500,6000];
y_lim= [1200,1600];

figure(5)
subplot(4,1,1)
% imagesc(neuron_input)
% ylabel('Neuron')
plot(curve_trace)
title('Curvature')
xlim(x_lim)

subplot(4,1,2)
plot(tspike(:,2),tspike(:,1),'k.')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Reservoir spikes')

subplot(4,1,3)
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
xlim(x_lim)

subplot(4,1,4)
[A , edges] = histcounts(tspike(:,2),1:1:length(neuron_input(1,:)));
A(1) = 0;
plot(A)
xlim(x_lim)
xlabel('Time /ms')
ylabel('Pop. activity')
title('Average population activity')
%% Calculate average firing rate of network
% maximum firing rate in barrel cortex during whisker deflection ~ 60 Hz.
R = ns/(N*T/1000);
%% Plot weights matrix
% OMEGA is the G*w_0 and E*BPhi is the Q*eta*phi, added togethes this gives
% weights matrix.subplot(2,2,1)
%wij = E*BPhi';
x_lim = [1 1600];
y_lim = [1 1600];
figure(6)
subplot(2,2,1)
imagesc(OMEGA)
axis 'square';
colorbar
title('Static weights')
xlim(x_lim)
ylim(y_lim)
subplot(2,2,2)
if dale == 1
    imagesc(wij)
else
    imagesc(E*BPhi')
end
xlim(x_lim)
ylim(y_lim)
colorbar
axis 'square';
title('Learned weights')
subplot(2,2,[3,4])
if dale == 1 
    imagesc(OMEGA + wij)
else
    imagesc(OMEGA + E*BPhi')
end
axis 'square';
xlim(x_lim)
ylim(y_lim)
colorbar
title('Static weights + learned weights')

%% Calculate input rate of thalamic neurons
n_spikes = 0;
for i =1:200
    n_spikes = n_spikes + length(SpikeTrainStruct{1, 1}.SpikeTimes{i,4});
end

R = n_spikes/(200*1.175);

%% ISI
ISI = [];
for n=1:N
    index = find(tspike(:,1)==n);
    spike_neuron(n).time = tspike(index,2);
    spike_neuron(n).ISI = diff(spike_neuron(n).time);
    ISI = [ISI diff(spike_neuron(n).time)'];
end
%%
figure(7)
subplot(2,1,1)
histogram(ISI,10000,'DisplayStyle','stairs')
xlim([0 200])
title('ISI')
% Coefficient of variation
Cv = zeros(N,1);
for n = 1:N
    Cv(n,1) = std(spike_neuron(n).time)/mean(spike_neuron(n).time);
end
subplot(2,1,2)
histogram(Cv,100,'DisplayStyle','stairs')
title('Coefficient of variation')



%% Colour code the spikes according to its BPhi value
c = zeros(ns,3);            

% This gives each BPhi value a RGB 'jet' value
C = jet(64); 
L = size(C,1);
Gs = round(interp1(linspace(min(BPhi(:)),max(BPhi(:)),L),1:L,BPhi));
H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
H_s = squeeze(H);

%% Assign each colour to neuron
tspike = tspike_train;
for neur = 1:N
   which = find(tspike(:,1) == neur);
   c(which,:) = H_s(neur,:) .* ones(length(which),3);
end

%% Remove zeros from spike array
    % Remove zeros
tspike_nz = [];
for t=1:length(tspike)
    if tspike(t,:) ~= [0 0]
        tspike_nz = [tspike_nz, tspike(t,:)];
    end
end
tspike_nz2 = zeros(length(tspike_nz)/2 , 2);

for i=1:length(tspike_nz)/2
    tspike_nz2(i,:) = [tspike_nz( 2*i -1) tspike_nz( 2*i)];
end
%% Plot coloured spikes with each neuron showing its BPhi value
x_lim= [0,T];
y_lim= [0,1600];
i = T/dt;
figure(1)
subplot(4,1,1)
% imagesc(neuron_input)
% ylabel('Neuron')
plot(curve_trace,'LineWidth',1.5)
title('Curvature')
xlim(x_lim)

subplot(4,1,2)
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1.5), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
xlim(x_lim)

subplot(4,1,3)
scatter(tspike_nz2(:,2),tspike_nz2(:,1),1,c,'filled')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
xlabel(' Time /ms')
title('Reservoir spikes')

subplot(4,1,4 )
BPhi_order = sort(BPhi);
imagesc(BPhi)
colorbar
set(gca,'YDir','normal')
%ylim(y_lim)
colormap('jet')
axis 'square';
title('BPhi')

suptitle('Spikes coloured with BPhi')

%% Sort neurons based on BPhi value

[~,indx] = sort(BPhi);


for i = 1:N
    %spike_time = zeros(sum(tspike_nz2(:,1) == i));
    ind = find(tspike_nz2(:,1) == i);
    spike_time = tspike_nz2(ind,2);
    spike{i}.time = spike_time;
    
end
        
%%
subplot(2,1,1)
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1.5), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
xlim(x_lim)

subplot(2,1,2)
for in = 1:N
    neuron = in*ones(length(spike{indx(in)}.time) ,1);
    plot(spike{indx(in)}.time, neuron,'k.')
    hold on  
end
xlabel('Time /ms')
title('spikes sorted on BPhi value')
hold off


%%
histogram(BPhi,100,'DisplayStyle','stairs')

%% Plot spikes with each neuron showing its sum of input weights 

Ein_sum = sum(Ein,2);

histogram(Ein_sum,100)

c_e = zeros(ns,3);            

% This gives each Ein_sum value a RGB 'jet' value
C = jet(64); 
L = size(C,1);
Gs = round(interp1(linspace(min(Ein_sum(:)),max(Ein_sum(:)),L),1:L,Ein_sum));
H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
H_s = squeeze(H);

%% Assign each colour to neuron
tspike = tspike_train;
for neur = 1:N
   which = find(tspike(:,1) == neur);
   c_e(which,:) = H_s(neur,:) .* ones(length(which),3);
end

%% Plot Input
x_lim= [0,T];
y_lim= [0,1600];
i = T/dt;
figure(1)
subplot(4,1,1)
% imagesc(neuron_input)
% ylabel('Neuron')
plot(curve_trace,'LineWidth',1.5)
title('Curvature')
xlim(x_lim)

subplot(4,1,2)
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1.5), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
xlim(x_lim)

subplot(4,1,3)
scatter(tspike_nz2(:,2),tspike_nz2(:,1),1,c_e,'filled')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Reservoir spikes')
xlabel(' Time /ms')

subplot(4,1,4 )
imagesc(Ein_sum)
colorbar
set(gca,'YDir','normal')
%ylim(y_lim)
colormap('jet')
axis 'square';
title('Ein sum')

suptitle('Spikes coloured with input')

%% Plot spikes with each neuron showing its feedback weights 

histogram(Efb,100)

c_fb = zeros(ns,3);            

% This gives each Ein_sum value a RGB 'jet' value
C = jet(64); 
L = size(C,1);
Gs = round(interp1(linspace(min(Efb(:)),max(Efb(:)),L),1:L,Efb));
H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
H_s = squeeze(H);

%% Assign each colour to neuron
tspike = tspike_train;
for neur = 1:N
   which = find(tspike(:,1) == neur);
   c_fb(which,:) = H_s(neur,:) .* ones(length(which),3);
end

%% Plot
x_lim= [0,T];
y_lim= [0,1600];
i = T/dt;
figure(1)
subplot(4,1,1)
% imagesc(neuron_input)
% ylabel('Neuron')
plot(curve_trace,'LineWidth',1.5)
title('Curvature')
xlim(x_lim)

subplot(4,1,2)
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1.5), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
xlim(x_lim)

subplot(4,1,3)
scatter(tspike_nz2(:,2),tspike_nz2(:,1),1,c_fb, 'filled')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Reservoir spikes')
xlabel(' Time /ms')

subplot(4,1,4 )
imagesc(Efb)
colorbar
set(gca,'YDir','normal')
%ylim(y_lim)
colormap('jet')
axis 'square';
title('E feedback')

suptitle('Spikes coloured with feedback')


%% Plot Bphi spikes and feedback spikes

x_lim= [0,T];
y_lim= [0,1600];
circle_size = 1;

subplot(4,2,[1 2])
plot((1:1:T),z_all(1:1:T),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r','LineWidth',1.5), hold off
legend('Target', 'Output')
title('Testing')
ylabel('Ouput')
xlim(x_lim)

subplot(4,2,[3 4])
scatter(tspike_nz2(:,2),tspike_nz2(:,1),circle_size,c,'filled')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Output coloured spikes')

subplot(4,2,[5 6])
scatter(tspike_nz2(:,2),tspike_nz2(:,1),circle_size,c_fb, 'filled')
ylim(y_lim)
xlim(x_lim)
ylabel('Neuron')
title('Feedback coloured spikes')
xlabel(' Time /ms')

subplot(4,2,7 )
BPhi_order = sort(BPhi);
imagesc(BPhi)
colorbar
set(gca,'YDir','normal')
%ylim(y_lim)
colormap('jet')
axis 'square';
title('BPhi')

subplot(4,2,8 )
imagesc(Efb)
colorbar
set(gca,'YDir','normal')
%ylim(y_lim)
colormap('jet')
axis 'square';
title('E feedback')


























