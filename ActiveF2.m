function [F,curr,x,spikes] = ActiveF2(c,A,omega,tau,dt,lambda,lambda1)

[J,T] = size(c);
I = eye(J);
ker = [0,exp(-lambda*(0:T/100)*dt)];

curr = [c;zeros(J,T)];
spikes = zeros(1,T);
Tpast = 0;
F = [];
x = zeros(size(curr));
while Tpast<T
    for j = 1:2*J
        y = conv(curr(j,Tpast+1:end),ker)*dt;
        x(j,Tpast+1:end) = y(1:T-Tpast);
    end
    t = find(sqrt(sum(x(:,Tpast+1:end).^2,1))>=omega,1)+Tpast;
    if ~isempty(t)
        d = x(:,t)/norm(x(:,t))*omega;
        F = [F;d'/omega];
        spikes(t) = 1;
        Tpast = t
        %add synaptic current:        
        Ds = (lambda*I+A)*d(1:J,:) + (lambda1*I+A)*tau^(-1)*d(J+1:end,:);
        curr(:,Tpast+1:end) = curr(:,Tpast+1:end) + [-Ds;tau*Ds]*exp(-lambda1*(0:T-Tpast-1)*dt);
    else
        Tpast = T;
    end
end
