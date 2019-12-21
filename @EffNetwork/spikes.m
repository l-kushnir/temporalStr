function o = spikes(this)
o = zeros(this.N,this.Tcurr);
for i = 1:size(this.PastSpikes,2)
    tsp = this.PastSpikes(1,i);
    n = this.PastSpikes(2,i);
    o(n,tsp) = 1;
end