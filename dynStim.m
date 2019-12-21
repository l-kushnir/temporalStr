function stim = dynStim(c0,A,lambda,TL,TR,dt)

c1 = zeros(size(c0,1),TL);
    
c1(:,1:TR) = c0/(dt*TR)*ones(1,TR);
stim1 = Stimulus('dt',dt,'control',c1,'A',A);
c = lambda*stim1.toRepresent;
stim = Stimulus('dt',dt,'control',c,'A',-lambda);
stim.input = stim.control;