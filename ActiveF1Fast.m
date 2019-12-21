function [F,x,spikes] = ActiveF1Fast(c,omega,dt,lambda)

[J,T] = size(c);
I = eye(J);
ker = [0,exp(-lambda*(0:T/100)*dt)];

curr = c;
spikes = zeros(1,T);
Tpast = 0;
F = [];
x = zeros(size(curr));
for j = 1:J
    y = conv(curr(j,:),ker)*dt;
    x(j,:) = y(1:T);
end
while Tpast<T
    x(:,Tpast+1:end) = x(:,Tpast+1:end) - x(:,Tpast+1)*exp(-lambda*(0:T-Tpast-1)*dt);
    t = find(sqrt(sum(x(:,Tpast+1:end).^2,1))>=omega,1)+Tpast;
    if ~isempty(t)
        d = x(:,t)/norm(x(:,t))*omega;
        F = [F;d'/omega];
        spikes(t) = 1;
        Tpast = t
    else
        Tpast = T;
    end
end