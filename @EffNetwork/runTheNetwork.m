function [input,Delta,nsp,g,PastV,r,noise] = runTheNetwork(this,input,varargin)

dt = this.dt;
[N,TL] = size(input);

p = inputParser;
addOptional(p,'RecordHis',1);
addOptional(p,'Delta',min(100,TL));
addOptional(p,'RunTillSpike',0);
parse(p,varargin{:});
Delta = p.Results.Delta;

if ~isempty(this.PastV)
    V = this.PastV(:,end);
else
    V = zeros(N,1);
end

if ~isempty(this.PastIn)
    input = [this.PastIn,input];
else
    input = [zeros(N,1),input];
end

PastV = [];
PastSpikes = [];

if ~isempty(this.fastSignal)
    fastsignal = this.fastSignal;
else
    fastsignal = zeros(N,1);
end

for m = 1:length(this.Slowker)
    if ~isempty(this.slowSignal)
        slowsignal{m} = this.slowSignal{m}(:,end);
    else
        slowsignal{m} = zeros(N,1);
    end
end

if ~isempty(this.filtSpikes)
    r = this.filtSpikes;
else
    r = zeros(this.N,1);
end

spikes = zeros(N,1);
if ~isempty(this.PastSpikes)
    if this.PastSpikes(1,end) == this.Tcurr
        spikes(this.PastSpikes(2,end)) = 1;
    end
end
Tcurr = this.Tcurr;
Th = this.T;
THis = [];

nsp = zeros(this.N,1); 
g = cell(size(this.Slowker));
noise = [];

while (size(input,2)>1)&&(~p.Results.RunTillSpike || sum(nsp)==0)
    
    
    %%% Feedforward contribution 
    ffcontr = convn(input(:,1:Delta),this.Inker((1:Delta)*dt))*dt;   %input is now shifted by one time step 
    ffcontr = ffcontr(:,1:Delta);

    %%% Fast recurrent contribution 
   
    fastsignal = (fastsignal/dt+spikes/dt^2)*this.Fastker((0:Delta)*dt)*dt;
    fastcontr = convn(fastsignal(:,1:Delta),this.Inker((1:Delta)*dt))*dt;  
    fastcontr = this.Wfast*fastcontr(:,1:Delta);

    %%% Slow recurrent contribution
    
    %slowsignal = diag(slowsignal+spikes)*this.Slowker((0:Delta)*dt);
    slowsignalTot = zeros(this.N,Delta+1); % This is because kernel goes 0:Delta below.                                      
    for m = 1:length(this.Slowker)       
        slowsignal{m} = slowsignal{m}/dt*this.Slowker{m}((1:Delta+1)*dt)*dt + spikes/dt*this.Slowker{m}((0:Delta)*dt)*dt;
        slowsignalTot = slowsignalTot + this.Wslow{m}*slowsignal{m};
    end
    slowcontr = convn(slowsignalTot(:,1:Delta),this.Inker((1:Delta)*dt))*dt;  
    slowcontr = slowcontr(:,1:Delta);
    
    %%% Noise
    dnoise = randn(N,Delta)*this.sigmaNoise;
    
    %%% Updating V
    
     if size(this.Inker((1:Delta)*dt),1)==1
        V = V*this.Inker((2:Delta+1)*dt)+ffcontr+fastcontr+slowcontr+dnoise;
    else 
        V = diag(V)*this.Inker((2:Delta+1)*dt)+ffcontr+fastcontr+slowcontr+dnoise;    
     end
    Tsp = find(sum(V>Th*ones(1,Delta),1),1);
    
    %%% No spike:
    
    if isempty(Tsp)
        PastV = [PastV,V];
        THis = [THis,Th*ones(1,Delta)];
        V = V(:,Delta);
        input = input(:,Delta+1:end);
        Tcurr = Tcurr+Delta;
        fastsignal = fastsignal(:,Delta+1); 
        for m = 1:length(this.Slowker)
            g{m} = [g{m},slowsignal{m}(:,2:Delta+1)];
            slowsignal{m} = slowsignal{m}(:,Delta+1);
        end
        if size(this.Inker((1:Delta)*dt),1)==1
            r = [r,r(:,end)*this.Inker((2:Delta+1)*dt)+spikes*this.Inker((1:Delta)*dt)];
        else
            r = [r,diag(r(:,end))*this.Inker((2:Delta+1)*dt)+diag(spikes)*this.Inker((1:Delta)*dt)];
        end
        spikes = zeros(N,1);  
        noise = [noise,dnoise];                  
        Delta = min([2*Delta,size(input,2)-1,1000]); 
        
    %%% If there is a spike
    
    else
        [~,n] = max(V(:,Tsp)-Th);
%        V(:,Tsp) = min(V(:,Tsp),Th+0.001); %To bring neurons to reset potential after they spike
        Tcurr = Tcurr+Tsp;
        PastSpikes = [PastSpikes,[Tcurr;n]];
        nsp(n) = nsp(n)+1;
        PastV = [PastV,V(:,1:Tsp)];
        THis = [THis,Th*ones(1,Tsp)];    
        V = V(:,Tsp);
        fastsignal = fastsignal(:,Tsp+1);
        for m = 1:length(this.Slowker)        
            g{m} = [g{m},slowsignal{m}(:,2:Tsp+1)];
            slowsignal{m} = slowsignal{m}(:,Tsp+1);            
        end
        if size(this.Inker((1:Tsp)*dt),1)==1
            r = [r,r(:,end)*this.Inker((2:Tsp+1)*dt)+spikes*this.Inker((1:Tsp)*dt)];
        else
            r = [r,diag(r(:,end))*this.Inker((2:Tsp+1)*dt)+diag(spikes)*this.Inker((1:Tsp)*dt)];
        end
        spikes = zeros(N,1);
        spikes(n) = 1;
        input = input(:,Tsp+1:end);
        Delta = min(max(50,floor(Delta/2)),size(input,2)-1);
        noise = [noise,dnoise];
    end
  Tcurr %%%% To track the progress
end    

this.T = Th;
this.PastIn = input(:,1);
input = input(:,2:end);

if p.Results.RecordHis
    this.PastSpikes = [this.PastSpikes,PastSpikes];
    this.PastV = [this.PastV,PastV];   
    if isempty(this.slowSignal)
        for m = 1:length(this.Slowker)             
            this.slowSignal{m} = g{m};
        end
    else
        for m = 1:length(this.Slowker)  
            this.slowSignal{m} = [this.slowSignal{m},g{m}];
        end      
    end
else
    if ~isempty(PastSpikes)
        this.PastSpikes = PastSpikes(:,end);
    end
    this.PastV = PastV(:,end);
    for m = 1:length(this.Slowker)
        this.slowSignal{m} = g{m}(:,end);
    end
end
this.fastSignal = fastsignal;

this.Tcurr = Tcurr;
this.filtSpikes = r(:,end);
r = r(:,2:end);    

    
    
    
    
    
    
    
    
    
    


