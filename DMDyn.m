J = 2;      % Dimensionality of the input
dt = 0.0001;  % Time step
I = eye(J);    % Identity matrix
tau1 = 0.02*I;  % Matrix tau_1 from eq (53)
tau2 = [0.025 0;0 0.035];   % Matrix tau_2 from eq (53)
tau1B = [0.035 0;0 0.025];  % Matrix \bar{tau_1}  from eq (53) 
tau2B = [0.035 0;0 0.02];   % Matrix \bar{tau_2} bar from eq (53)    
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
[Fact,curr,x,spikes] = ActiveF3(stim.input,A,omega,tau1,tau2,tau1B,tau2B,dt,lambda,lambda1,lambda2);
                                % Determining the 
                                %6-dim selectivity vectors of the neurons that
                                %will fire in response to the input
F6 = bunch(Fact,0.06);          %adding neurons with feedforward weights in the 
                                %directions adjacent to Fact
F6 = diag(1./sqrt(sum(F6.^2,2)))*F6;  % Normalizing the rows of the matrix F6

%% Initialize the Network
N = size(F6,1);
n = EffNetwork('N',N,'dt',dt,'LambdaFast',lambda,'LambdaSLow',[lambda1,lambda2]);
%%%%%%%%%% Feedforward connections and  Network connectivity and Thresholds

Ftot = F6;          % F concatinated with F^int and with \bar{F^int} from eq (53)
Ff = Ftot(:,1:J);   % Feedforward connections (matrix F)
Fin = Ftot(:,J+1:2*J);  % first 2 internal dimensions (matrix F^int) 
FinB = Ftot(:,2*J+1:3*J);  % second pair of internal dimensions (matrix \bar{F^int}) 
F = Ftot;
n.F = Ff;
D = omega*F'*diag((sum(F.^2,2)).^(-1/2));  % 6-dim selectivity vectors \hat d concatinated
n.Wfast = -F*D;                            % Matrix of fast connections
n.T = omega*sqrt(sum(F.^2,2));             % Thresholds


%%%%%%%%%%%%%% Computing the 2 matrices of slow connections  
B = (tau1^(-1)*tau2 - tau1B^(-1)*tau2B)^(-1)*(tau1^(-1)*D(J+1:2*J,:) - tau1B^(-1)*D(2*J+1:3*J,:));
C = (tau2^(-1)*tau1 - tau2B^(-1)*tau1B)^(-1)*(tau2^(-1)*D(J+1:2*J,:) - tau2B^(-1)*D(2*J+1:3*J,:));
D1 = 1/(lambda2-lambda1)*((lambda2*I+A)*(lambda*I+A)*(D(1:J,:)+B) + (lambda1*I+A)*((lambda+lambda2-lambda1)*I+A)*C);
D2 = 1/(lambda1-lambda2)*((lambda1*I+A)*(lambda*I+A)*(D(1:J,:)+C) + (lambda2*I+A)*((lambda+lambda1-lambda2)*I+A)*B);
n.Wslow{1} = -Ff*D1 + Fin*tau1*D1 + FinB*tau1B*D1;
n.Wslow{2} = -Ff*D2 + Fin*tau2*D2 + FinB*tau2B*D2;

%%%%%%%%%%%% Run the network
nspTot = zeros(N,1);
chat = [];
TotInIn = [];
Ehat = [];

for k=1:TL/10000
    k
    [~,~,nsp,g,~,r] = n.runTheNetwork(n.F*stim.input(:,(k-1)*10000+1:k*10000));  
    if k==1
        V = n.PastV;
    end
    n.removeHist;
    nspTot = nspTot + nsp;   %Number of spikes by neuron
    chat = [chat,D1*g{1}+D2*g{2}];  %Reconstruction of the input 
    TotInIn = [TotInIn,[tau1*D1*g{1}+tau2*D2*g{2};tau1B*D1*g{1}+tau2B*D2*g{2}]];
                                    % Recurrent input in the additional 4-dim internal
                                    % space
    Ehat = [Ehat,D*r];  %Reconstruction of the leaky integral of total (feedforward + slow recurrent) input
end

% Computing the leaky integral of the total input (6-dim)
TotInConv = [];
NspTot = sum(nspTot);
TotIn = [stim.input(:,1:length(chat)) - chat;TotInIn];

for j = 1:6
    TotInConv(j,:) = conv(TotIn(j,:),n.Inker(dt:dt:dt*10000))*dt;
end
TotInConv = TotInConv(:,1:length(TotIn));
%save('DM','chat','Ehat','V','NspTot','TotIn','TotInConv','nspTot','F6','N','stim')

figure %figure 4c
plot((1:TL)*dt,stim.input')
hold on
plot((1:TL)*dt,chat')
xlim([2,5])
hold off
%saveas(1,'DerMatch.fig')


