function [ E, OMEGA, BPhi, tspike_train,current,ns, v_rec ] = LIF_network_spikes( neuron_input, zx, N, G, Q, alpha, dt, td, dale,INPUT,BIAS,Bphi,Omega,Efb,wij )
% This function is very similar to train_LIF_network.m, except that here the
% input is in the form of 'weighted' thalamic spikes. This is in the form
% of a NxT array with N the number of nerons in the reservoir and T the
% total time of the simulation.
% 
% neuron_input: thalamic spike trains to be classified. 
% zx: the correct target function for the output to match.
% N: the number of neurons.
% G: the synaptic scaling parameter. This controls the chaos in the
% reservoir. 
% Q: the scaling of the feedback weights. This can 'tame' the chaos.
% alpha: this is the learning parameter, sets the rate of weights change. 
% dt: is the integration time step.
% step: number of iteration step RLS is performed.
% td: decay time constant, usually 10, 20, or 30 ms. 
% dale: set to 1 if you want weights to respect Dale's law.
% train: set to 1 if you want to train the network. Otherwise it will
% just let the reservoir run so the activity can be checked. 
% BIAS = -65; %Set the BIAS current, can help decrease/increase firing rates.  0 is fine. 
% INPUT: set ot either 1 or 0 to give thalamic input or not.
%
% This function gives as output
% Ein: the Nxin input weights
% E: the Nxk feedback weights
% OMEGA: the NxN recurrent weights of the reservoir
% BPhi: the Nxk learned output weights.
% tpike_train: the spikes of the network.
% z_plot: the target function the network has to match.
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
T = length(neuron_input(1,:)); 
nt = round(T/dt); imin = round(500/dt); % start RLS
icrit = nt;% stop RLS
%% Storage parameters and some others
Pinv = eye(N)*alpha; %initialize the correlation weight matrix for RLMS
p = 0.1; %Set the network sparsity 
k = min(size(zx));
IPSC = zeros(N,1); %post synaptic current storage variable hmm
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
JD = 0*IPSC; %storage variable required for each spike time 
tspike_train = zeros(4*nt,2); %Storage variable for spike times 
ns = 0; %Number of spikes, counts during simulation  
z = zeros(k,1);  %Initialize the approximant 
z_plot=zeros(nt,1); % This is to store the target values for a plot
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
v_ = v;  %v_ is the voltage at previous time steps  
v_rec = zeros(T ,10); % record voltage of 5 neurons
current = zeros(nt,k);  %storage variable for output current/approximant 
i = 1; 
tlast = zeros(N,1); %This vector is used to set  the refractory times 
I_c = zeros(nt,1);

if Bphi == 0 % This statement checks if new weights have to be made or not.
    BPhi = zeros(N,k); %The initial matrix that will be learned by FORCE method set the row average weight to be zero, explicitly.
    E = (2*rand(N,k)-1)*Q; % Feedback weights.
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
else
    BPhi = Bphi;
    E = Efb;
    OMEGA = Omega;
end
%% Make neuron input on the time scale of the network intergation
% As the neural input is every 1 ms and the integration step is 0.05 ms
% this loop makes the neural input T/dt long but only presents the spikes
% every 1 ms step. 

% neuron_dt = zeros(N,T/dt);
% 
% for t = 1:nt
%     if mod(t,1/dt) == 0
%         neuron_dt(:,t) = neuron_input(:,t*dt);
%     end  
% end
%% Run network with no RLS just input spikes

ilast = 1;
in = 1; 

for i = ilast:1:nt
    
    if dale == 1
        I = IPSC + wij*r + BIAS;
    else
        I = IPSC + E*z + BIAS ; %Neuronal Current
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
    
    IPSC = IPSC*exp(-dt/tr) + h*dt;
    if mod(i,1/dt) == 0
        h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td) + INPUT*neuron_input(:,i*dt)/(tr*td);  % THE LAST TERM ARE THE THALAMIC SPIKES
    else
        h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td);
    end
    
    r = r*exp(-dt/tr) + hr*dt;
    hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
 
    %Implement RLMS with the FORCE method
    z = BPhi'*r; %approximant
    v = v + (30 - v).*(v>=vpeak);
    v = v + (vreset - v).*(v>=vpeak); %reset with spike time interpolant implemented.
    current(i,:) = z;
    
    % record 5 random neuron's voltages. 
    if mod(i,1/dt) == 0
        v_rec(in,:) = v(1:length(v_rec(1,:)));
        in = in + 1;
    end
        
end

end



