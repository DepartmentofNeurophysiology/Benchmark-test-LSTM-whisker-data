%% IZ FORCE classification with whisker input
clear all
clc 
%%
load('an197522_2013_03_08_session');
load('Kernels');
load('whiskingstruct');

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
plot_trial(5,dat);

%% Create supervisor and input

for i=2:length(s.trialIds)
    if dat(i).lick_times ~= 0
        dat(i).z = dat(i).timeVec;
        %dat(i).z(1)=0;
        for t=1:length(dat(i).timeVec)
            if dat(i).pole_times(2) < dat(i).timeVec(t) && dat(i).timeVec(t) < dat(i).lick_times(1) % Time between pole not more in reach and first lick
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
%%

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
        X1 = [X1, dat(i).thetaVec];          % takes theta values per correct trial
        X2 = [X2, dat(i).kappaVec];          % takes kappa values per correct trial
        zx = [zx, dat(i).z];
        count = count + 1;
        end
    end  
end
X = [X1 ; X2];
[conv_whisk, conv_curve] = convolve_kernel_whisk_curve( KernelStruct, length(s.trialIds), dat);
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
T = 1000000; %Total time 
dt = 0.04; %integration time step 
nt = round(T/dt);
N =  2000;  %number of neurons, 500 works
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
%BPhi = zeros(N,1);
%-----Initialization---------------------------------------------
v = vr+(vpeak-vr)*rand(N,1); %initial distribution 
v_ = v; %These are just used for Euler integration, previous time step storage
IPSC = zeros(N,1);
Q = 5*1E0^3; %Scale feedback term, Q in paper
E = (2*rand(N,1)-1)*Q; %scale feedback term
WE2 = 5*10^2; %scale input weights
Psi = 2*pi*rand(N,1); 
Ein = [cos(Psi),sin(Psi)]*WE2;
z = 0; 
tspike = zeros(nt,2);
ns = 0;
% RLS parameters.
Pinv = eye(N)*30;
step = 10;
imin = round(200/dt);
icrit = length(z1);% round(20000/dt);
current = zeros(nt,1);
RECB = zeros(nt,5);
REC = zeros(nt,10);
i=1;
ilast = i ;
I_list = zeros(nt,1);
X(isnan(X)) = 0;
z1=1*zx;
%% SIMULATION
z_t = z1;
for i = ilast:1:length(z1)
% EULER INTEGRATE
I = IPSC + E*z+ BIAS + Ein*X(i,:)'; 
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
 err = z - z_t(i);
% if mod(i,step)==1
% if i > imin 
%  if i < icrit 
%      %if z_t(i) ~= 0
%          cd = Pinv*r;
%          BPhi = BPhi - (cd*err');
%          Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
%     % end
%  end 
% end 
% end
% COMPUTE S, APPLY RESETS
u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
v = v+(vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
v_ = v;  % sets v(t-1) = v for the next itteration of loop
REC(i,:) = [v(1:5)',u(1:5)']; 
current(i,:) = z; 
I_list(i,:) = mean(I);
RECB(i,:)=BPhi(1:5);
if mod(i,round(100/dt))==1 
drawnow

figure(2)
plot(dt*(1:1:i),current(1:1:i,1),'r.','Linewidth',0.8), hold on 
plot(dt*(1:1:i),z_t(1:1:i),'g--','Linewidth',2), hold off
ylim([-2,2])
xlim([dt*i-1000,dt*i])
xlabel('Time')
ylabel('Network Response')
legend('Network Output','Target Signal')
pause(0.4)
end   
end






