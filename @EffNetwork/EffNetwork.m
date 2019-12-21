classdef EffNetwork < handle
    
    properties
        dt = 0.001;
        N %number of neurons
        W0 %0/1 connectivity matrix
        Wfast % connectivity, fast connections
        Wslow % connectivity, slow connections
        F     % feedforward connectivity 
        lambdaV = 1;
        lambdaSlow = 1;
        sigmaNoise %amplitude of Gaussian noise in the dynamics of V        
        T %vector of thresholds
        RefrAmp % vector of refractory amplitudes
        Tmemory = 0.01 %in seconds

        Inker
        Fastker
        Slowker
        Refractoryker 
        
        PastV 
        THis
        PastSpikes
        fastSignal
        slowSignal
        PastIn
        filtSpikes
        recentSpikes
        Tcurr %current time at the end of PastV, in timesteps
    end
    
    methods
        % constructor
        function this = EffNetwork(varargin)    
        %%% Neuron's dynamics    
            this.Inker = @(t) (t>0).*exp(-this.lambdaV*(t-this.dt));            
            this.Fastker = @(t) double(t==0);
            %this.Slowker = @(x) exp(-this.lambdaSlow*x);
            %this.Slowker = @(x) exp(-this.lambdaSlow*(x-this.dt))*diag(x>0);
            
           % this.Refractoryker = @(x) -(x>0).*exp(log(1-this.lambdaV*this.dt)*(x-this.dt)/this.dt);
            
        %%% Network size, connectivity and thershold                
            p = inputParser;
            addOptional(p,'N',42);                       
            addOptional(p,'ConnSparseness',1);
            addOptional(p,'dt',[]);
            addOptional(p,'sigmaNoise',0);
            addOptional(p,'Wfast',[], @(W) all(diag(W)<0));
            addOptional(p,'Wslow',[]);
            addOptional(p,'T',[])
            addOptional(p,'Refr',[]);
            addOptional(p,'LambdaSlow',[]);
            addOptional(p,'LambdaFast',10);
            addOptional(p,'W0',[]);
            parse(p,varargin{:});            
            N = p.Results.N;
            this.N = N;
            this.Tcurr = 0;
            if ~isempty(p.Results.dt)
                this.dt = p.Results.dt;
            end
            if isempty(p.Results.LambdaSlow)
                this.lambdaSlow = [];
            else
                this.lambdaSlow = p.Results.LambdaSlow;
            end
            
            defW0 = double(rand(N,N)<p.Results.ConnSparseness);
            defW0 = defW0-diag(diag(defW0))+diag(ones(1,N));
            if isempty(p.Results.W0)
                W0 = defW0;
            else
                W0 = p.Results.W0;
            end
            this.W0 = W0;
            if isempty(p.Results.Wfast)
                defWfast = randn(N,N).*W0;
                defWfast = defWfast-diag(diag(defWfast))-diag(abs(randn(1,N)));
                this.Wfast = defWfast.*W0;
            else
                this.Wfast = p.Results.Wfast.*W0;
            end
            if isempty(p.Results.Wslow)
                for m = 1:length(this.lambdaSlow) 
                    this.Wslow{m} = randn(N,N).*W0; 
                end
            else
                this.Wslow = p.Results.Wslow.*W0;
            end
            if isempty(p.Results.T)
                this.T = 0.1*ones(N,1)+randn(N,1)*0.01; 
            else
                this.T = p.Results.T;
            end
                    
            this.lambdaV = p.Results.LambdaFast;
            this.sigmaNoise = p.Results.sigmaNoise;  
            
            for m = 1:length(this.lambdaSlow)  
                %this.Slowker{m} = @(t) (t>0).*exp(-this.lambdaSlow(m)*(t-this.dt));
                this.Slowker{m} = @(t) (t>=0).*exp(-this.lambdaSlow(m)*t);
            end
        end
        [input,Delta,nsp,g,PastV,r,noise] = runTheNetwork(this,input,varargin) 
        NetReset(this)
        randomInit(this)
        s = spikes(this)
        removeHist(this)
    end
end 
        







































            