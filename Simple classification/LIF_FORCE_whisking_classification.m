clear all
clc
%% Load data
load('whiskingstruct')
%load('all_variables_40')

%% Create input and target
[ X_curve, X_angle, z_all, target_label] = make_output_target( [ 5 8 ], 10, dat );
plot(z_all)
hold on
plot(X_curve)
hold on
%plot(target_label)
plot(X_angle)
hold off
xlabel( 'time /ms')
ylabel( 'activity')
legend('target output', 'curvature', 'angle','location','southwest')
%% Network parameters
N = 50;  %Number of neurons 
dt = 0.005; %time step
tref = 2; %Refractory time constant in milliseconds 
tm = 10; %Membrane time constant 
vreset = -65; %Voltage reset 
vpeak = -40; %Voltage peak. 
rng(1);
td = 20; tr = 2;
T = length(z_all); 

imin = round(1000/dt); nt = round(T/dt); icrit = nt; step = 20; Q = 150; G = 100; %% IMPORTANT PARAMETERS

alpha = dt*10^3 ; %Sets the rate of weight change, too fast is unstable, too slow is bad as well.  
Pinv = eye(N)*alpha; %initialize the correlation weight matrix for RLMS
p = 0.1; %Set the network sparsity 

k = min(size(z_all));
IPSC = zeros(N,1); %post synaptic current storage variable 
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
JD = 0*IPSC; %storage variable required for each spike time 
tspike = zeros(4*nt,2); %Storage variable for spike times 
ns = 0; %Number of spikes, counts during simulation  
z = zeros(k,1);  %Initialize the approximant 
z_plot=zeros(nt,1); % This is to store the target values for a plot
w_dot = zeros(nt,1); % Storage of the change in weights
 
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
v_ = v;  %v_ is the voltage at previous time steps  
RECB = zeros(nt,10);  %Storage matrix for the synaptic weights (a subset of them) 
OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p); %The initial weight matrix with fixed random weights  
BPhi = zeros(N,k); %The initial matrix that will be learned by FORCE method
%set the row average weight to be zero, explicitly.
for i = 1:1:N 
    QS = find(abs(OMEGA(i,:))>0);
    OMEGA(i,QS) = OMEGA(i,QS) - sum(OMEGA(i,QS))/length(QS);
end

%% Input weights
WE2 = 5*10^2; %scale input weights
%Ein = -WE2 + (WE2+WE2)*rand(N,1); % If only one input

Psi = 2*pi*rand(N,1); 
Ein = [cos(Psi),sin(Psi)]*WE2; % If two inputs

%%
E = (2*rand(N,k)-1)*Q;  %n
REC2 = zeros(nt,20);
REC = zeros(nt,10);
current = zeros(nt,k);  %storage variable for output current/approximant 
i = 1; 
tlast = zeros(N,1); %This vector is used to set  the refractory times 
BIAS = vpeak; %Set the BIAS current, can help decrease/increase firing rates.  0 is fine. 
I_c = zeros(nt,1);
%% Test trained weights
load('weights')
imin = nt;
 
%%
zx = z_all;
X = [X_curve ; X_angle];
ilast = i; 
in = 1;
%icrit = ilast;
for i = ilast:1:nt
    if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented
        in=in+1;
    end
    I = IPSC + E*z + BIAS + Ein*X(:,in); %Neuronal Current 
    dv = (dt*i>tlast + tref).*(-v+I)/tm; %Voltage equation with refractory period
    v = v + dt*(dv);   
    index = find(v>=vpeak);  %Find the neurons that have spiked   
    %Store spike times, and get the weight matrix column sum of spikers
    if length(index)>0
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i];
        ns = ns + length(index);  % total number of spikes so far
    end
    
    tlast = tlast + (dt*i -tlast).*(v>=vpeak);  %Used to set the refractory period of LIF neurons
    
    % Code if the rise time is 0, and if the rise time is positive
    if tr == 0
        IPSC = IPSC*exp(-dt/td)+   JD*(length(index)>0)/(td);
        r = r *exp(-dt/td) + (v>=vpeak)/td;
    else
        IPSC = IPSC*exp(-dt/tr) + h*dt;
        h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td);  %Integrate the current
        
        r = r*exp(-dt/tr) + hr*dt;
        hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
    end  
    %Implement RLMS with the FORCE method
    z = BPhi'*r; %approximant
    err = z - zx(:,in); %error
    %RLMS
    if mod(i,step)==1
        if i > imin
            if i < icrit
                if zx(in) ~= 0
                    cd = Pinv*r;
                    BPhi = BPhi - (cd*err');
                    Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
                end
            end
        end
    end    
    v = v + (30 - v).*(v>=vpeak);
%     REC(i,:) = v(1:10); %Record a random voltage
    v = v + (vreset - v).*(v>=vpeak); %reset with spike time interpolant implemented.
    current(i,:) = z;
%     RECB(i,:) = BPhi(1:10);
%     REC2(i,:) = r(1:20);
    z_plot(i,1) = zx(1,in);
    I_c(i,:) = mean(Ein*X(:,in)); % Record mean input
    w_dot(i,:) = sum(abs(cd*err')); % Record change in weights
    % Plot
    if mod(i,round(100/dt))==1
        drawnow
        
        figure(1)
        plot(tspike(1:1:ns,2),tspike(1:1:ns,1),'k.')
        xlim([dt*i-1000,dt*i])
        ylim([0,N])
        
        figure(2)
        
        subplot(2,1,1)
        plot(dt*(1:1:i),z_plot(1:1:i),'g--','LineWidth',2), hold on
        plot(dt*(1:1:i),current(1:1:i,:),'r.','LineWidth',0.8), hold off    
        xlim([dt*i-1000,dt*i])
        ylim([-4,4])
        
        subplot(2,1,2)
        plot(dt*(1:1:i),I_c(1:1:i,1),'b','Linewidth',0.8)
        xlim([dt*i-1000,dt*i])
        ylim([-100,100])
        xlabel('Time')
        ylabel('Mean input')
        %xlim([dt*i,dt*i])
%         figure(3)
%         plot(dt*(1:1:i),REC(1:1:i,1:5),'.')
%         xlim([dt*i-100,dt*i])
%         pause(0.1)
    end
end

%%
save('weights_trial_5_8_15_18_curve_angle', 'Ein', 'E', 'OMEGA', 'BPhi')

%% Find times were the target is
index_p1 = find([0 diff(target_label)]==1 );
index_m1 = find([0 diff(target_label)]==-1);
index = [index_p1 index_m1];
index = sort(index);
index=index./dt;

%% Calculate mean output for trials
mean_out = zeros(n,1);
for trial = 1:n
    mean_out(trial) = mean(current(index(2*trial-1):index(2*trial)));
end

%%
correct_out = zeros(n,1);
%%
values = diff(target_label);
values = values(values ~= 0);
%%
for trial = 1:n
    correct_out(trial) = values( 2*trial-1);
end
%% Plot correct and mean output
plot(correct_out, 'go','Linewidth',2)
hold on
plot(mean_out, 'r*', 'Linewidth',2)
hold off
xlabel('trial number')
ylabel('response')
legend('Correct response', 'Mean output')
ylim([-3, 3])
%% Plot results
subplot(3,1,1)
%plot(tspike(1:1:ns,2),tspike(1:1:ns,1),'k.')
plot(tspike(:,2),tspike(:,1),'k.')
ylim([0,N])
xlim([0,length(z_all)])
title('Network spikes')
ylabel('Neuron')
subplot(3,1,2)
plot(dt*(1:1:i),current(:,:),'r')
title('Network output (no RLS)')
ylabel('Output')
xlim([0,length(z_all)])
ylim([-3.5 3.5])
subplot(3,1,3)
plot(z_all,'g')
hold on
plot(X_add)
hold off
xlim([0,length(z_all)])
title('Curvature trace and correct output')
xlabel('Time /ms')
ylabel('Output')
legend('Correct ouput', 'curvature','Location', 'southeast')    
suptitle('Trained network output for ten trials')









