function [wij, Ein, E, OMEGA, BPhi, tspike_train, tspike_test,z_plot,current ] = train_LIF_network( Xin, Xtrain, zx, ztrain, N, G, Q, alpha, Win, dt, step, td, dale,train,BIAS)
% This function takes as input 
% Xin: the current input to be classified (e.g. whisker trace) This should
% be of the form Xin = X_curve for  one inputo or Xin = [X_curve ; X_angle] for
% two inputs. 
% zx: the correct target function for the output to match.
% N: the number of neurons.
% G: the synaptic scaling parameter. This controls the chaos in the
% reservoir. 
% Q: the scaling of the feedback weights. This can 'tame' the chaos.
% alpha: this is the learning parameter ('sets the rate of weights change.
% Win: the scaling of the input weights. 
% dt: is the integration time step.
% step: number of iteration step RLS is performed.
% td: decay time constant, usually 10, 20, or 30 ms. 
% dale: set to 1 if you want weights to respect Dale's law.
% train: set to 1 if you want to train the network. Otherwise it will
% just let the reservoir run so the activity can be checked. 
% BIAS = -65; %Set the BIAS current, can help decrease/increase firing rates.  0 is fine. 
%
% This function gives as output
% Ein: the Nxin input weights
% E: the Nxk feedback weights
% OMEGA: the NxN recurrent weights of the reservoir
% BPhi: the Nxk learned output weights.
%
% The inputs defined above are the most interesting to change to see how
% the network functions. All the other parameters can be adjusted below. 

%% Network parameters
tref = 2; %Refractory time constant in milliseconds 
tm = 10; %Membrane time constant 
vreset = -65; %Voltage reset 
vpeak = -40; %Voltage peak. 
rng(1);
tr = 2;
T = length(zx); 
nt = round(T/dt); imin = round(500/dt); % start RLS
icrit = nt;% stop RLS
%% Storage parameters and some others
Pinv = eye(N)*alpha; %initialize the correlation weight matrix for RLMS
p = 0.1; %Set the network sparsity 
k = min(size(zx));
IPSC = zeros(N,1); %post synaptic current storage variable 
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
JD = 0*IPSC; %storage variable required for each spike time 
tspike_train = zeros(4*nt,2); %Storage variable for spike times 
ns = 0; %Number of spikes, counts during simulation  
z = zeros(k,1);  %Initialize the approximant 
z_plot=zeros(nt,1); % This is to store the target values for a plot
w_dot = zeros(nt,1); % Storage of the change in weights
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
v_ = v;  %v_ is the voltage at previous time steps  
RECB = zeros(nt,10);  %Storage matrix for the synaptic weights (a subset of them) 
BPhi = zeros(N,k); %The initial matrix that will be learned by FORCE method set the row average weight to be zero, explicitly.
E = (2*rand(N,k)-1)*Q;  
REC2 = zeros(nt,20);
REC = zeros(nt,10);
current = zeros(nt,k);  %storage variable for output current/approximant 
i = 1; 
tlast = zeros(N,1); %This vector is used to set  the refractory times 
I_c = zeros(nt,1);

if dale == 1
    % This creates a static weights matrix that respects Dales law.
    Ne = round(0.8*N);
    Ni = round(0.2*N);
    OMEGA_E =  abs(G*(randn(Ne,N)).*(rand(Ne,N)<p)/(sqrt(N)*p));
    OMEGA_I = -abs(G*(randn(Ni,N)).*(rand(Ni,N)<p)/(sqrt(N)*p));
    OMEGA = [ OMEGA_E ; OMEGA_I]';
    % The next loop sets the sample mean to zero by adjusting the inhibitory
    % weights
    for i = 1:1:N
        QS = find(abs(OMEGA(i,:))>0);
        QS_ind_in = find(QS > Ne); % These are the inhibitory weights
        QS_in = QS(QS_ind_in);
        OMEGA(i,QS_in) = OMEGA(i,QS_in) - sum(OMEGA(i,QS))/length(QS_in); % Sets sum of input weights to each neuron to zero.
    end
else
    OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p);
    for i = 1:1:N
        QS = find(abs(OMEGA(i,:))>0);
        OMEGA(i,QS) = OMEGA(i,QS) - sum(OMEGA(i,QS))/length(QS);
    end
end

%% Create Input weights
Psi = 2*pi*rand(N,1); 
if min(size(Xin)) == 1
    Ein = -Win + (Win+Win)*rand(N,1); % If only one input
elseif min(size(Xin)) == 2
    Ein = [cos(Psi),sin(Psi)]*Win; % If two inputs
end

%% Network Training
if train == 1
    ilast = 1;
    in = 1;
    for i = ilast:1:nt
        
        if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented
            in=in+1;
        end
        if in > T
            break
        end
        % Make sure the RLS trained weight matrix also respects Dales law.
        if dale == 1
            wij = E*BPhi';
            % This makes sure that the approximant interacts with the reservoir in a way which obeys Dales law.
            for ni = 1:Ne
                ind = find(wij(:,ni) < 0);
                wij(ind,ni) = 0;
            end
            for ni = Ne+1:N
                ind = find(wij(:,ni) > 0);
                wij(ind,ni) = 0;
            end
            I = IPSC + wij*r + BIAS + Ein*Xin(:,in); %Neuronal Current
        else
            I = IPSC + E*z + BIAS + Ein*Xin(:,in); %Neuronal Current
            wij = 0;
        end
        
        dv = (dt*i>tlast + tref).*(-v+I)/tm; %Voltage equation with refractory period
        v = v + dt*(dv);
        index = find(v>=vpeak);  %Find the neurons that have spiked
        %Store spike times, and get the weight matrix column sum of spikers
        if length(index)>0
            JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
            tspike_train(ns+1:ns+length(index),:) = [index,0*index+dt*i];
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
                    if zx(in) ~= 0 % Makes sure there is only RLS when the target is non-zero
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
        %    I_c(i,:) = mean(Ein*X(:,in)); % Record mean input
        %   w_dot(i,:) = sum(abs(cd*err')); % Record change in weights
    end
    
    %% Plot target and output for training
    figure(1)
    plot(dt*(1:1:i),z_plot(1:1:i),'g--','LineWidth',2), hold on
    plot(dt*(1:1:i),current(1:1:i,:),'r.','LineWidth',0.8), hold off
    ylim([ -4 4])
    xlim([ 1000 T])
    title(' Training')
end

%% Create new input and target

T = length(ztrain); 
nt = round(T/dt); 

%% Reset storage variables

IPSC = zeros(N,1); %post synaptic current storage variable 
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
JD = 0*IPSC; %storage variable required for each spike time 
tspike_test = zeros(4*nt,2); %Storage variable for spike times 
ns = 0; %Number of spikes, counts during simulation  
z = zeros(k,1);  %Initialize the approximant 
z_plot=zeros(nt,1); % This is to store the target values for a plot
w_dot = zeros(nt,1); % Storage of the change in weights
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
v_ = v;  %v_ is the voltage at previous time steps  
RECB = zeros(nt,10);  %Storage matrix for the synaptic weights (a subset of them)  
REC2 = zeros(nt,20);
REC = zeros(nt,10);
current = zeros(nt,k);  %storage variable for output current/approximant 
i = 1; 
tlast = zeros(N,1); %This vector is used to set  the refractory times 

%% Network testing
ilast = 1; 
in = 1;
if train ~=1
    Ein = 0;
end
for i = ilast:1:nt
    if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented
        in=in+1;
    end
    if in > T
        break
    end
    I = IPSC + E*z + BIAS + Ein*Xtrain(:,in); %Neuronal Current 
    dv = (dt*i>tlast + tref).*(-v+I)/tm; %Voltage equation with refractory period
    v = v + dt*(dv);   
    index = find(v>=vpeak);  %Find the neurons that have spiked   
    %Store spike times, and get the weight matrix column sum of spikers
    if length(index)>0
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        tspike_test(ns+1:ns+length(index),:) = [index,0*index+dt*i];
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
    z = BPhi'*r; %approximant
    v = v + (30 - v).*(v>=vpeak);
    v = v + (vreset - v).*(v>=vpeak); %reset with spike time interpolant implemented.
    current(i,:) = z;
    z_plot(i,1) = ztrain(1,in);
end
%% Plot for testing
figure(2)
plot(dt*(1:1:i),z_plot(1:1:i),'g--','LineWidth',2), hold on
plot(dt*(1:1:i),current(1:1:i,:),'r.','LineWidth',0.8), hold off
ylim([ -4 4])
xlim([ 1000 T])
title(' Testing')

%% Average firing rate
f_r = ns / (T /1000) ;
show = ['Average firing rate is ',num2str(f_r)];
disp(show)
end





