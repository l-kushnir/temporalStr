function [] = removeHist(this)

this.PastSpikes = this.PastSpikes(:,end);
for i = 1:length(this.slowSignal)
    this.slowSignal{i} = this.slowSignal{i}(:,end);
end
this.PastV = this.PastV(:,end);
