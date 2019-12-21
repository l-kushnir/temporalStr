function randomInit(this)
[F,S,~] = svd(this.Wslow{1}+this.Wslow{2});
Dim = sum(diag(S)>0.01*S(1,1)); %0.01 is chosen randomly. Dim is the rank of Omega, the dimensionality of the box.
F = F(:,1:Dim);
d0 = randn(Dim,1);
d0 = d0/sqrt(d0'*d0)*0.5; % half-way between the center and the border of the box

V = F*d0;
g1 = rand(this.N,1)/sqrt(this.N)*2;
g2 = rand(this.N,1)/sqrt(this.N)*2;

if isempty(this.PastV)
    this.PastV = V;
else
    this.PastV(:,end) = V;
end
this.slowSignal = {g1,g2};

