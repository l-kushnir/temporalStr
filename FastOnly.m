J = 2;  % Dimensionality of the input
dt = 0.0001; % Time step
I = eye(J);  % Identity matrix
omega = 0.05; % Tolerated error
lambda = 10;  % Neuronal time constant 
lambda1 = 2;  % First synaptic time constant (not used in this script)
lambda2 = 1.2;   % First synaptic time constant (not used in this script)

%%%%% Predictable stimulus
TL = 1000000;    % Total duration of the input in time steps 
c0 = [-0.3;0.96];  % Initial value of the input
TR = 1;          % Number of time steps to reach c0    
A = [-0.12 -0.036;1 0];  % Dynamical matrix for the input evolution
stim = dynStim(c0,A,lambda,TL,TR,dt);  % Creating the input

%%%%%%% Feedforward connections
[Fact,x,spikes] = ActiveF1Fast(stim.input,omega,dt,lambda); % Determines the 
                                %selectivity vectors of the neurons that
                                %will fire in response to the input
F2 = Fact;
F2 = diag(1./sqrt(sum(F2.^2,2)))*F2;  % Normalizing the feedforward weights
%% Initialize the Network
N = size(F2,1);
n = EffNetwork('N',N,'dt',dt,'LambdaFast',lambda,'LambdaSLow',[lambda1,lambda2]);

%%%%%%%%%% Network connectivity and Thresholds
F = F2;
n.F = F2;
D = omega*F'*diag((sum(F.^2,2)).^(-1/2));  % The decoder matrix (selectivity vectors d^hat)
n.Wfast = -F*D;                            % Fast connectivity matrix 
n.T = omega*sqrt(sum(F.^2,2));             % Thresholds

n.Wslow{1} = zeros(N,N);
n.Wslow{2} = zeros(N,N);

%%%%%%%%%%%% Run the network
nspTot = zeros(N,1);
chat = [];

for k=1:TL/10000
    k
    [~,~,nsp,g,~,r] = n.runTheNetwork(n.F*stim.input(:,(k-1)*10000+1:k*10000));  
    if k==1
        V = n.PastV;
    end
    n.removeHist;
    nspTot = nspTot + nsp;    %Number of spikes by neuron
    chat = [chat,lambda*D*r]; %Reconstruction of the input
end

NspTot = sum(nspTot);

%save('FastOnly','chat','V','NspTot','F2','nspTot','N','lambda','stim')

%Plot the Input
figure
plot((1:TL)*dt,stim.input')

% Plot the Input and the reconstruction
figure
plot((1:TL)*dt,stim.input')
hold on
plot((1:TL)*dt,chat')
xlim([2,5])
hold off
%saveas(1,'FastOnly.fig')

%Plot the leaky integral of the input and the reconstruction
figure;
plot((1:TL)*dt,stim.toRepresent'*lambda)
hold on
plot((1:TL)*dt,chat')
xlim([2,5])
hold off

%Plot the traces of voltages
figure;
plot((1:10000)*dt,V(randperm(N,5),:)')
xlim([0.6,0.9])

NspTot = sum(nspTot)  %Total number of spikes


