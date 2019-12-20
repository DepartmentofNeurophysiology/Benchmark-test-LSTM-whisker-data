function [ X_add_curve, X_add_angle, z_all, target_label] = make_output_target( trials, n, dat )
% This function takes a number of trials as input and creates an input
% signal of n trials long where the input trials of angle traces and curvature traces
% are randomly arranged. The corresponding target output is also returned
% as is the target label.

% This loop makes sure there are no incorrect trials
for i = 1:length(trials)
    if dat(trials(i)).correct == 0
        show = [' Trial ',num2str(trials(i)),' is incorrect.'];
        disp(show)
        break
    end
end
% Take trials and standardise and remove zeros and store in structure
for i = 1:length(trials)
    [X(i).cnz , X(i).anz] = standard_nozero( trials(i),dat);
end

% Some parameters for the pulse
reset = 800; % Time inbetween the trials
pulse_length = 400; % Length of pulse
amp = 3; % Amplitude of pulse
decay = 80; % decay of exponential pulse
constant = 200; % How long the pulse is kept constant
X_curve = -0.5*ones(1000,1)';
X_angle = -0.5*ones(1000,1)';
z_all = zeros(1000,1)';
start_early = 500; % How long to start before the end of the input
target_label = zeros(1000,1)';

for t=1:n
    trial = randi([1,length(trials)],1,1);
    l_t = length(X(trial).cnz) + reset;
    z_t = zeros(l_t,1);
    z_t(length(X(trial).cnz)-start_early:1:length(X(trial).cnz)+constant) = dat(trial).pole*amp*ones(start_early+constant+1,1);
    z_t(length(X(trial).cnz)+constant:1:length(X(trial).cnz)+pulse_length+constant) =  dat(trial).pole*amp*exp(-(0:1:pulse_length)./decay);
    target_label = [target_label dat(trial).pole*ones(length(X(trial).cnz) +pulse_length+constant,1)' zeros(reset -pulse_length-constant,1)' ];
    z_all = [z_all z_t'];
    X_curve = [X_curve X(trial).cnz zeros(reset,1)'];
    X_angle = [X_angle X(trial).anz zeros(reset,1)'];
end

X_add_curve = X_curve + 0.5;
X_add_angle = X_angle + 0.5;
end
