%% IZ FORCE classification with whisker input
clear all
clc 
%%
load('an197522_2013_03_08_session');
load('Kernels');
load('whiskingstruct');
%load('weights');

%% Load data of whisking
for i=1:length(s.trialIds)
    [dat(i).thetaVec, dat(i).kappaVec ,dat(i).timeVec, dat(i).pole,dat(i).touch_times,...
        dat(i).pole_vec,dat(i).pole_times,dat(i).lick_times ]...
     = angle_curve_position( s, i); % this creates structure for all trials with relevant info.
    if dat(i).pole_vec(1)==1 || dat(i).pole_vec(2)==1 % Label correct trials
        dat(i).correct = 1 ;        
    else
        dat(i).correct = 0 ;
    end
end
    
%% Plot to check

plot_trial(15,dat);


%% Create supervisor and input

for i=2:length(s.trialIds)
    if dat(i).lick_times ~= 0
        dat(i).z = dat(i).timeVec;
        %dat(i).z(1)=0;
        for t=1:length(dat(i).timeVec)
            %if dat(i).pole_times(2) < dat(i).timeVec(t) && dat(i).timeVec(t) < dat(i).lick_times(1) % Time between pole not more in reach and first lick
            if dat(i).pole_times(2) < dat(i).timeVec(t) % This puts target signal to +/- 1 for the rest of the trial
                if dat(i).pole == 1
                    dat(i).z(t) = 1;            % Label correct output for pole
                elseif dat(i).pole == -1
                    dat(i).z(t) = -1;           % Label correct output for pole
                end
            else
                dat(i).z(t) = 0;                % Rest state?
            end
        end
    else
        dat(i).z = 0; 
    end
end
    
save('whiskingstruct', 'dat')   
%% Create input signals and target signals
% Create one input array with all trials and create a supervisor that 
% correctly labels the desired output for each time point of the input
% array.
X1 = [];
X2 = [];
zx = [];
count = 1;
for i=1:length(s.trialIds)
    if dat(i).correct == 1 
        if dat(i).lick_times ~= 0   
        X1 = [X1,zeros(1,5000), dat(i).thetaVec];          % takes theta values per correct trial
        X2 = [X2,zeros(1,5000), dat(i).kappaVec];          % takes kappa values per correct trial
        zx = [zx,zeros(1,5000), dat(i).z];
        count = count + 1;
        end
    end  
end
X = [X1 ; X2];
%% Convolve with kernels
[conv_whisk, conv_curve] = convolve_kernel_whisk_curve( KernelStruct, length(s.trialIds), dat);
%% Create signal for only two trials
% Take trial 5 and 8, pole is -1 and +1 respectively
% Only take curvature for now, maybe learn to ignore angle later..?
% Select trials 't1' and 't2'
t1 = '5'; % touch 
t2 = '8'; % no touch
noNAN = dat('t1').kappaVec;
noNAN(isnan(noNAN)) = 0;
trial_1 = (noNAN-mean(noNAN))/std(noNAN);

noNAN = dat('t2').kappaVec;
noNAN(isnan(noNAN)) = 0;
trial_2 = (noNAN-mean(noNAN))/std(noNAN);

pole_1 = dat('t1').pole_times;
pole_2 = dat('t2').pole_times;
time_1 = dat('t1').timeVec;
time_2 = dat('t2').timeVec;
%%
X_1 = zeros(time_1(end),1);
X_2 = zeros(time_2(end),1);
i=1;
for t=1:length(X_1)
    X_1(t) = trial_1(i);
    if mod(t,2) == 0
        i = i + 1;
    end
end
%%
i = 1;
for t=1:length(X_2)
    X_2(t) = trial_2(i);
    if mod(t,2) == 0
        i = i + 1;
    end
end
%% Only input between poles
for t=1:length(X_1)
    if t < pole_1(1) || t > pole_1(2)
        X_1(t) = 0;
    end
end
for t=1:length(X_2)
    if t < pole_2(1) || t > pole_2(2)
        X_2(t) = 0;
    end
end
%% Remove zeros
Xnz_1 = [];
Xnz_2 = [];
for t=1:length(X_1)
    if X_1(t) ~= 0
        Xnz_1 = [Xnz_1, X_1(t)];
    end
end
for t=1:length(X_2)
    if X_2(t) ~= 0
        Xnz_2 = [Xnz_2, X_2(t)];
    end
end
    
%% Target pulse
reset = 800; % Time inbetween the trials
pulse_length = 400; % Length of the sine pulse for target function
amp = 3; % Amplitude of pulse
decay = 80;
y=amp*exp(-(0:1:pulse_length)./decay);
plot(y)
%%
constant = 200;
n=10; % Number of times a trial is repeated
X_all = -0.5*ones(1000,1)';
z_all = zeros(1000,1)';
target_label = zeros(1000,1)';
for t=1:n
    trial = randi([0,1],1,1);
    if trial == 0 % take trial 5
        l_t = length(Xnz_5) + reset;
        z_t = zeros(l_t,1);
        %z_t(pole_5(2):1:pole_5(2)+pulse_length) = amp*sin((-pi/pulse_length)*(0:1:pulse_length));
        %z_t(pole_5(2):1:pole_5(2)+pulse_length) =  -1.5*fpdf(0.05*(0:1:pulse_length),10,5);
        %z_t(length(Xnz_5):1:length(Xnz_5)+pulse_length) =  -3*fpdf(0.05*(0:1:pulse_length),10,5);
        z_t(length(Xnz_5):1:length(Xnz_5)+constant) = -amp*ones(constant+1,1);
        z_t(length(Xnz_5)+constant:1:length(Xnz_5)+pulse_length+constant) =  -amp*exp(-(0:1:pulse_length)./decay);
        target_label = [target_label -ones(length(Xnz_5) +pulse_length+constant,1)' zeros(reset -pulse_length-constant,1)' ];
        z_all = [z_all z_t'];
        X_all = [X_all Xnz_5 zeros(reset,1)'];
    else % take trial 8
        l_t = length(Xnz_8) + reset;
        z_t = zeros(l_t,1);
        %z_t(pole_8(2):1:pole_8(2)+pulse_length) = amp*sin((pi/pulse_length)*(0:1:pulse_length));
        %z_t(pole_8(2):1:pole_8(2)+pulse_length) =  1.5*fpdf(0.05*(0:1:pulse_length),10,5);
        %z_t(length(Xnz_8):1:length(Xnz_8)+pulse_length) =  3*fpdf(0.05*(0:1:pulse_length),10,5);
        z_t(length(Xnz_8):1:length(Xnz_8)+constant) = amp*ones(constant+1,1);
        z_t(length(Xnz_8)+constant:1:length(Xnz_8)+pulse_length+constant) =  amp*exp(-(0:1:pulse_length)./decay);
        target_label = [target_label ones(length(Xnz_8) +pulse_length+constant,1)' zeros(reset-pulse_length-constant,1)' ];
        z_all = [z_all z_t'];
        X_all = [X_all Xnz_8 zeros(reset,1)'];
    end
end
%%
X_add = X_all + 0.5;
plot(z_all)
hold on
plot(X_add)
hold on
plot(target_label)
hold off

%% Test other target
%Load the data set for classification.  P contains class, z contains the points.  
inputfreq = 4; %Present a data point; 
j =0;
Xin = zeros(nt,2); 
z2 = zeros(nt,1); 
nx = round(1000/(dt*inputfreq)); 
z2 = 0.5*abs(sin(2*pi*(1:1:nt)*dt*inputfreq/2000)); 
%% Create supervisor and inputs.  
j = 1; 
k2 = 0;
for i =1:1:nt 
 z2(i) = z2(i)*mod(k2,2);        % This loop makes the points nx=6250 long... and then 6250 points of rest.
if mod(i,nx)==1 
    k2 = k2 + 1; 
end
end
%%
X=X';
%% INITIALISE SIMULATION PARAMETERS
T = length(X_all); %Total time 
dt = 0.005; %integration time step in ms
nt = round(T/dt);
N =  20;  %number of neurons, 500 works
% Izhikevich Parameters
C = 250;
vr = -60; 
b = 0; 
k = 2.5; 
vpeak = 30; 
vreset = -65;
vt = vr+40-(b/k); %threshold 
Er = 0; %Reversal Potential 
u = zeros(N,1); 
a = 0.01;
d = 200; 
tr = 2; 
td = 20; 
p = 0.1; %sparsity 
G =6*10^3; %scale weight matrix  
BIAS = 1000; %Bias current
OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(p*sqrt(N));
% Initialize currents, FORCE method, other parameters
IPSC = zeros(N,1); %post synaptic current 
h = zeros(N,1);
r = zeros(N,1);
hr = zeros(N,1);
JD = zeros(N,1);

%-----Initialization---------------------------------------------
v = vr+(vpeak-vr)*rand(N,1); %initial distribution 
v_ = v; %These are just used for Euler integration, previous time step storage
IPSC = zeros(N,1);
Q = 5*10^3; %Scale feedback term, Q in paper
E = (2*rand(N,1)-1)*Q; %scale feedback term
z = 0; 
tspike = zeros(nt,2);
ns = 0;
% RLS parameters.
current = zeros(nt,1);
RECB = zeros(nt,5);
REC = zeros(nt,20);
i=1;%1;%5000/dt;
ilast = i ;
I_list = zeros(nt,1);
%X(isnan(X)) = 0;
%z1=1*zx;
in=1;
z_plot=zeros(nt,1);
I_c = zeros(nt,5);
w_dot = zeros(nt,1);

Q = 1*10^0; %Scale feedback term, Q in paper
E = (2*rand(N,1)-1)*Q; %scale feedback term
G =1*10^3; %scale weight matrix, this induces chaos 
BIAS = 1000; %Bias current
OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(p*sqrt(N));
imin =  1000;%round(500/dt);
icrit = nt;
WE2 = 5*10^3; %scale input weights
Psi = 2*pi*rand(N,1); 
%Ein = [cos(Psi),sin(Psi)]*WE2;
Ein = randi([-WE2 WE2],N,1);
BPhi =  0.1*ones(N,1); % Ones or zeros?
Pinv = eye(N)*10;
step = 5;
f_r = zeros(nt,1);
%% SIMULATION
z_t = z_all; % Set target function
X = X_add';
tic
ilast = 1;
for i = ilast:1:nt
    if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented, as dt = 0.05 ms
        in=in+1;
    end
    % EULER INTEGRATE
    I = IPSC + E*z+ BIAS + Ein*X(in,:)';
    v = v + dt*(( k.*(v-vr).*(v-vt) - u + I))/C ; % v(t) = v(t-1)+dt*v'(t-1)
    u = u + dt*(a*(b*(v_-vr)-u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
    
    %
    index = find(v>=vpeak);
    if length(index)>0
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        %tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i];
        ns = ns + length(index);
    end
    if tr == 0
        IPSC = IPSC*exp(-dt/td)+   JD*(length(index)>0)/(td);
        r = r *exp(-dt/td) + (v>=vpeak)/td;
    else
        IPSC = IPSC*exp(-dt/tr) + h*dt;
        h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td);  %Integrate the current
        
        r = r*exp(-dt/tr) + hr*dt;
        hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
    end
    %Apply RLS
    z = BPhi'*r;
    err = z - z_t(in);
    if mod(i,step)==1
        if i > imin
            if i < icrit
                if z_t(in) ~= 0 %Only RLS when z is not zero%|| X(in,:) == 0 % No RLS when there is input
                    %if sum(X(in,:))==0 || z_t(in) == 1 || z_t(in) == -1 % No RLS when there is input before pole withdrawn
                    cd = Pinv*r;
                    BPhi = BPhi - (cd*err');
                    Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
                end
            end
        end
    end
    z_plot(i,1) = z_t(in);
    % COMPUTE S, APPLY RESETS
    u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component.
    v = v+(vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    v_ = v;  % sets v(t-1) = v for the next itteration of loop
    REC(i,1:10) = [v(1:5)',u(1:5)'];
    current(i,:) = z;
    %I_list(i,:) = mean(I);
    %I_c(i,:) = Ein*X(in,:)'(1:5);
    I_c(i,:) = mean(Ein*X(in,:)');
    %f_r(i,:) = mean(r);
    w_dot(i,:) = sum(abs(cd*err'));
        
    if mod(i,round(100/dt))==1
        drawnow
        figure(2)
        subplot(2,1,1)
        plot(dt*(1:1:i),current(1:1:i,1),'r.','Linewidth',0.8), hold on
        plot(dt*(1:1:i),z_plot(1:1:i),'g--','Linewidth',2), hold off
        ylim([-4,4])
        xlim([dt*i-1000,dt*i])
        xlabel('Time')
        ylabel('Network Response')
        %legend('Network Output','Target Signal')
        subplot(2,1,2)
        plot(dt*(1:1:i),I_c(1:1:i,1),'b','Linewidth',0.8)
        xlim([dt*i-1000,dt*i])
        ylim([-100,100])
        xlabel('Time')
        ylabel('Mean input')
%         subplot(3,1,3)
%         plot(dt*(1:1:i),f_r(1:1:i,1),'r','Linewidth',0.8)
%         xlim([dt*i-1000,dt*i])
%         %ylim([-100,100])
%         xlabel('Time')
%         ylabel('Mean firing rate')
        figure(3)
        for k=1:5
        plot(dt*(1:1:i),REC(1:1:i,k))
        hold on
        end
        hold off
%         hold on 
%         plot(dt*(1:1:i),REC(1:1:i,2))
%         hold on
%         plot(dt*(1:1:i),REC(1:1:i,3))
%         hold on
%         plot(dt*(1:1:i),REC(1:1:i,4))
%         hold on
%         plot(dt*(1:1:i),REC(1:1:i,5))
%         hold off
        xlim([dt*i-1000,dt*i])
        xlabel('Time')
        ylabel('Membrane potential')
        title('five neurons')
        
        figure(4)
        plot(dt*(1:1:i),w_dot(1:1:i,1))
        %xlim([dt*i-1000,dt*i])
        xlabel('Time')
        ylabel('Magnitude of change in weights')
        
        pause(0.1)
    end
    if mod(i,1000) == 0
        i
        i/nt
    end
end
toc
%%
E=zeros(n,1);
output_5 = zeros(length(Xnz_5) + pulse_length + constant,1);
output_8 = zeros(length(Xnz_8) + pulse_length + constant,1);

ti=1;
for t=1:nt
    if target_label(t) == 1
        output_5(ti) = current(t);
        ti = ti + 1;
    elseif target_label(t) == -1
        output_8(ti) = current(t);
        ti = ti + 1;
    else
        ti=1;
    end
end
        
%%

index_p1 = find([0 diff(target_label)]==1 );
index_m1 = find([0 diff(target_label)]==-1);

index = [index_p1 index_m1];

index = sort(index);
index=index./dt;

%%
mean_out = zeros(n,1);

for trial =1:n
    mean_out(trial) = mean(current(index(2*trial-1):index(2*trial)));
    
end

   



