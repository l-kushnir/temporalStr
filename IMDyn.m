J = 2;      % Dimensionality of the input
dt = 0.0001;  % Time step
I = eye(J);    % Identity matrix
tau = 0.02*I;  %matrix tau from eq (43)
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
[Fact,curr,x,spikes] = ActiveF2(stim.input,A,omega,tau,dt,lambda,lambda1);% Determines the 
                                %4-dim selectivity vectors of the neurons that
                                %will fire in response to the input
F4 = bunch(Fact,0.06);          %adding neurons with feedforward weights in the 
                                %directions adjacent to Fact
F4 = diag(1./sqrt(sum(F4.^2,2)))*F4; % Normalizing the rows of the matrix F4

%% Initialize the Network
N = size(F4,1);
n = EffNetwork('N',N,'dt',dt,'LambdaFast',lambda,'LambdaSLow',[lambda1,lambda2]);

%%%%%%%%%% Feedforward connections and  Network connectivity and Thresholds
Ftot = F4;         % F concatinated with F^int from eq (43)

Ff = Ftot(:,1:J);  % Feedforward weights 
Fin = Ftot(:,J+1:2*J);  % internal directions
F = [Ff,Fin];
n.F = Ff;

D = omega*F'*diag((sum(F.^2,2)).^(-1/2)); % Matrix of 4-dim selectivity vectors d^hat concatinated together
n.Wfast = -F*D;                           % Fast connections  
n.T = omega*sqrt(sum(F.^2,2));            % Thresholds  

Ds = (lambda*I+A)*D(1:J,:) + (lambda1*I+A)*tau^(-1)*D(J+1:end,:); % Slow decoder
n.Wslow{1} = -Ff*Ds + Fin*tau*Ds;  % Matrix of slow connections
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
end

% Computing the leaky integral of the total input (4-dim)
TotInConv = [];
TotIn = [stim.input(:,1:length(chat)) - chat;tau*chat];

for j = 1:4
    TotInConv(j,:) = conv(TotIn(j,:),n.Inker(dt:dt:dt*10000))*dt;
end
TotInConv = TotInConv(:,1:length(TotIn));
%save('IM','chat','Ehat','V','NspTot','TotIn','TotInConv','F4','lambda','lambda1','N','stim')

figure %figure 4a
plot((1:TL)*dt,stim.input')
hold on
plot((1:TL)*dt,chat')
xlim([2,5])
hold off
%saveas(1,'InputMatch.fig')

NspTot = sum(nspTot); %Total number of spikes


