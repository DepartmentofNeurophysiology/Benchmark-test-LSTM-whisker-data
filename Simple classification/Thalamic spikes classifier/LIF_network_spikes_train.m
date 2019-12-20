function [ E, OMEGA, BPhi, tspike_train,current, wij, BPhi_rec,delta_wij,r_rec ] = LIF_network_spikes_train( neuron_input, zx, N, G, Q, alpha, dt, td, dale,INPUT,BIAS,step,rep,BPhi,Omega,Efb)
% This function is very similar to LIF_network_spikes.m except that here the
% network is FORCE trained. This means that the output weights are updated
% through the RLS method so that the output matches the correct target
% function. Most parameters are already described in LIF_network_spikes.m

% rep: the number of times the training of the trial is repeated with the
% same neural input.
%% Network parameters
tref = 2; %Refractory time constant in milliseconds 
tm = 10; %Membrane time constant 
vreset = -65; %Voltage reset 
vpeak = -40; %Voltage peak. 
rng(1);
tr = 2;
T = length(zx); 
nt = round(T/dt); 
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
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
v_ = v;  %v_ is the voltage at previous time steps  

current = zeros(nt,k);  %storage variable for output current/approximant 
i = 1; 
tlast = zeros(N,1); %This vector is used to set the refractory times 
BPhi_rec = zeros(T,10);
delta_wij = zeros(T,1);
r_rec = zeros(T,10);

if BPhi == 0% This statement checks if new weights have to be made or not.
    BPhi = zeros(N,k); %The initial matrix that will be learned by FORCE method set the row average weight to be zero, explicitly.
    E = (2*rand(N,k)-1)*Q;  
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
    BPhi = BPhi;
    E = Efb;
    OMEGA = Omega;
end
%% Make neuron input on the time scale of the network intergation
% As the neural input is every 1 ms and the integration step is 0.05 ms
% this loop makes the neural input T/dt long but only presents the spikes
% every 1 ms step. 

neuron_dt = zeros(N,T/dt);
wij = 0;
for t = 1:nt
    if mod(t,1/dt) == 0
        neuron_dt(:,t) = neuron_input(:,t*dt);
    end  
end
%% Run network with RLS to update output weights. 
for repeat = 1:rep
ilast = 1;
in = 1;
ir = 1;

if dale == 1
for i = ilast:1:nt
    if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented
        in=in+1;
    end
    if in > T
        break
    end
    
    I = IPSC + wij*r + BIAS ; %Neuronal Current
    
    dv = (dt*i>tlast + tref).*(-v+I)/tm; %Voltage equation with refractory period
    v = v + dt*(dv);
    index = find(v>=vpeak);  %Find the neurons that have spiked
    %Store spike times, and get the weight matrix column sum of spikers
    if length(index)>0
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        tspike_train(ns+1:ns+length(index),:) = [index,0*index+dt*i];
        ns = ns + length(index);  % total number of spikes so far
    end
    
    %Implement RLMS with the FORCE method
    z = BPhi'*r; %approximant
    err = z - zx(:,in); %error
    
    % RLS
    if mod(i,step) == 1
        if i > imin
            if zx(in) ~= 0 % Makes sure there is only RLS when the target is non-zero
                cd = Pinv*r;
                BPhi = BPhi - (cd*err');
                Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
                
                
                    wij = E*BPhi';
                    wij(:,1:Ne) = abs(wij(:,1:Ne));
                    wij(:,Ne+1:N) = -abs(wij(:,Ne+1:N));
                
            end
        end
    end
    
    tlast = tlast + (dt*i -tlast).*(v>=vpeak);  %Used to set the refractory period of LIF neurons
    
    IPSC = IPSC*exp(-dt/tr) + h*dt;
    h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td) + INPUT*neuron_dt(:,i)/(tr*td);  % THE LAST TERM ARE THE THALAMIC SPIKES
    
    r = r*exp(-dt/tr) + hr*dt;
    hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);

    %Implement RLMS with the FORCE method
    z = BPhi'*r; %approximant
    v = v + (30 - v).*(v>=vpeak);
    v = v + (vreset - v).*(v>=vpeak); %reset with spike time interpolant implemented.
    current(i,:) = z;
    
    % record 5 random neuron's voltages. 
    if mod(i,1/dt) == 0
        BPhi_rec(ir,:) = BPhi(1:length(BPhi_rec(1,:)));
        ir = ir + 1;
    end
end

else 
for i = ilast:1:nt
    if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented
        in=in+1;
    end
    if in > T
        break
    end
    
    I = IPSC + E*z + BIAS ; %Neuronal Current
    
    dv = (dt*i>tlast + tref).*(-v+I)/tm; %Voltage equation with refractory period
    v = v + dt*(dv);
    index = find(v>=vpeak);  %Find the neurons that have spiked
    %Store spike times, and get the weight matrix column sum of spikers
    if length(index)>0
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        tspike_train(ns+1:ns+length(index),:) = [index,0*index+dt*i];
        ns = ns + length(index);  % total number of spikes so far
    end
    
    %Implement RLMS with the FORCE method
    z = BPhi'*r; %approximant
    err = z - zx(:,in); %error
    
    % RLS
    if mod(i,step) == 1   
        if zx(in) ~= 0 % Makes sure there is only RLS when the target is non-zero
            cd = Pinv*r;
            BPhi = BPhi - (cd*err');
            Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
        else
            cd = 0;
        end 
    end
    
    tlast = tlast + (dt*i -tlast).*(v>=vpeak);  %Used to set the refractory period of LIF neurons
    
    IPSC = IPSC*exp(-dt/tr) + h*dt;
    h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td) + INPUT*neuron_dt(:,i)/(tr*td);  % THE LAST TERM ARE THE THALAMIC SPIKES
    
    r = r*exp(-dt/tr) + hr*dt;
    hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);

    v = v + (30 - v).*(v>=vpeak);
    v = v + (vreset - v).*(v>=vpeak); %reset with spike time interpolant implemented.
    current(i,:) = z;
    
    % record 5 random neuron's voltages. 
    if mod(i,1/dt) == 0
        BPhi_rec(ir,:) = BPhi(1:length(BPhi_rec(1,:)));
        delta_wij(ir) = sum(abs(cd*err'));
        r_rec(ir,:) = r(1:10);
        ir = ir + 1;
    end
end
      
end
end

end

