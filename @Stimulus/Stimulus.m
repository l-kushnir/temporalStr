classdef Stimulus < handle
    properties
        dt = 0.001;
        dimension
        duration
        isLDS 
        A
        toRepresent
        input
        control
    end
    methods
        function this = Stimulus(varargin)
            p = inputParser;
            addOptional(p,'dimension',4);
            addOptional(p,'duration',100000);
            addOptional(p,'control',[]);
            addOptional(p,'toRepresent',[]);
            addOptional(p,'A',-1);
            addOptional(p,'sensoryInput',[]);
            addOptional(p,'dt',[]);
            parse(p,varargin{:});  
            if ~isempty(p.Results.dt)
                this.dt = p.Results.dt;
            end
            if all(isempty([p.Results.control,p.Results.toRepresent]))
                this.dimension = p.Results.dimension;
                this.duration = p.Results.duration;
                c = smoothrand(this.dimension,'duration',this.duration,'tau',this.duration/50);
                this.duration = size(c,2);
                this.control = c;    
                this.A = p.Results.A;
                this.toRepresent = this.dynVar(c,this.A);
            end
            if ~isempty(p.Results.control)
                this.control = p.Results.control;
                this.A = p.Results.A;
                this.toRepresent = this.dynVar(p.Results.control,this.A);
                this.dimension = size(this.toRepresent,1);
                this.duration = size(this.toRepresent,2);
            end
            if ~isempty(p.Results.toRepresent)
                this.toRepresent = p.Results.toRepresent;
                this.A = p.Results.A;
                this.control = this.inputFor(p.Results.toRepresent);
                this.dimension = size(this.toRepresent,1);
                this.duration = size(this.toRepresent,2);
            end            
            %this.input = this.control;
        end
        x = dynVar(this,c,A);
        c = inputFor(this,x,A);
         
    end
end