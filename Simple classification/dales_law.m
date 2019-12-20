G=10;
N=8;
p=1;


%%
% This creates a static weights matrix that respects Dales law.
% ??sample mean set to zero??
Ne = round(0.8*N);
Ni = round(0.2*N);
OMEGA_E =  abs(G*(randn(Ne,N)).*(rand(Ne,N)<p)/(sqrt(N)*p));
OMEGA_I = -abs(G*(randn(Ni,N)).*(rand(Ni,N)<p)/(sqrt(N)*p));
OMEGA = [ OMEGA_E ; OMEGA_I]';
%%
for i = 1:1:N
    QS = find(abs(OMEGA(i,:))>0);
    QS_ind_in = find(QS > Ne); % These are the inhibitory weights
    QS_in = QS(QS_ind_in);
    OMEGA(i,QS_in) = OMEGA(i,QS_in) - sum(OMEGA(i,QS))/length(QS_in); % Sets sum of input weights to each neuron to zero. 
end

%%
Q=10;
N =4;
Ne = round(0.8*N);
Ni = round(0.2*N);
BPhi = rand(N,1);
E = (2*rand(N,1)-1)*Q;  
wij = E*BPhi';
%% This makes sure that the approximant interacts with the reservoir in a way which obeys Dales law. 
for i = 1:Ne
    ind = find(wij(:,i) < 0);
    wij(ind,i) = 0;
end
for i = Ne+1:N
    ind = find(wij(:,i) > 0);
    wij(ind,i) = 0;
end



    
    