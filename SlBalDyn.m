J = 2;  % Dimensionality of the input
dt = 0.0001; % Time step
I = eye(J);  % Identity matrix
omega = 0.05;  % Tolerated error
lambda = 10;  % Neuronal time constant 
lambda1 = 2;  % First synaptic time constant (not used in this script)
lambda2 = 1.2;   % First synaptic time constant (not used in this script)


%%%%% Predictable stimulus
TL = 1000000;    % Total duration of the input in time steps 
c0 = [-0.3;0.96];  % Initial value of the input
TR = 1;          % Number of time steps to reach c0    
A = [-0.12 -0.036;1 0];  % Dynamical matrix for the input evolution
stim = dynStim(c0,A,lambda,TL,TR,dt);  % Creating the input


%%%%%%% Explored directions 
[Fact,curr,x,spikes] = ActiveF1(stim.input,omega,dt,lambda,lambda1); % Determines the 
                                %selectivity vectors of the neurons that
                                %will fire in response to the input
F2 = bunch(Fact,0.06);          %adding neurons with feedforward weights in the 
                                %directions adjacent to Fact
F2 = diag(1./sqrt(sum(F2.^2,2)))*F2;  % Normalizing the feedforward weights

%% Initialize the Network
N = size(F2,1);
n = EffNetwork('N',N,'dt',dt,'LambdaFast',lambda,'LambdaSLow',[lambda1,lambda2]);

%%%%%%%%%% Feedforward connections and  Network connectivity and Thresholds
F = F2;
n.F = F2;  %Feedforward connections

D = omega*F'*diag((sum(F.^2,2)).^(-1/2)); %Matrix of selectivity vectors d^hat concatinated together
n.Wfast = -F*D;                           %Fast connectivity matrix  
n.T = omega*sqrt(sum(F.^2,2));            %Thresholds  

Ds = lambda*D;                              % Slow decoder matrix
n.Wslow{1} = -F*Ds;                         % Slow connectivity matrix
n.Wslow{2} = zeros(N,N);

%%%%%%%%%%%% Run the network
nspTot = zeros(N,1);
chat = [];
Ehat = [];

for k=1:TL/10000
    k
    [~,~,nsp,g,~,r] = n.runTheNetwork(n.F*stim.input(:,(k-1)*10000+1:k*10000));  
    if k==1
        V = n.PastV;
    end
    n.removeHist;
    nspTot = nspTot + nsp; %Number of spikes by neuron
    chat = [chat,Ds*g{1}]; %Reconstruction of the input 
    Ehat = [Ehat,D*r];     %Reconstruction of the leaky integral of total (feedforward + slow recurrent) input 
                           % for fig 3b and 3c
end

% Computing the leaky integral of the total input
TotInConv = [];
TotIn = stim.input - chat;
for j = 1:2
    TotInConv(j,:) = conv(TotIn(j,:),n.Inker(dt:dt:dt*10000))*dt;
end
TotInConv = TotInConv(:,1:length(TotIn));
%save('SB','chat','Ehat','V','NspTot','TotIn','TotInConv','F2','lambda','lambda1','nspTot','N','stim')

figure % figure 3a
plot((1:TL)*dt,stim.input')
hold on
plot((1:TL)*dt,chat')
xlim([2,5])
hold off
saveas(1,'SlowBal.fig')

figure % figure 3b
plot((1:TL)*dt,lambda*TotInConv')
hold on
plot((1:TL)*dt,lambda*Ehat')
xlim([2,5])
ylim([-5,5])
hold off

figure % figure 3c
plot((1:TL)*dt,lambda*stim.toRepresent')
hold on
plot((1:TL)*dt,lambda*Ehat'+lambda*(-TotInConv+stim.toRepresent)')
xlim([2,5])
ylim([-5,5])
hold off

%Plot the traces of voltages
figure;
plot((1:10000)*dt,V(randperm(N,5),:)')
xlim([0.6,0.9])

NspTot = sum(nspTot);  %Total number of spikes
