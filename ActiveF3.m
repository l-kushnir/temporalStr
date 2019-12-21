function [F,curr,x,spikes] = ActiveF3(c,A,omega,tau1,tau2,tau1B,tau2B,dt,lambda,lambda1,lambda2)

[J,T] = size(c);
I = eye(J);
ker = [0,exp(-lambda*(0:T/100)*dt)];

curr = [c;zeros(2*J,T)];
spikes = zeros(1,T);
Tpast = 0;
F = [];
x = zeros(size(curr));
while Tpast<T
    for j = 1:3*J
        y = conv(curr(j,Tpast+1:end),ker)*dt;
        x(j,Tpast+1:end) = y(1:T-Tpast);
    end
    t = find(sqrt(sum(x(:,Tpast+1:end).^2,1))>=omega,1)+Tpast;
    if ~isempty(t)
        d = x(:,t)/norm(x(:,t))*omega;
        F = [F;d'/omega];
        spikes(t) = 1;
        Tpast = t
        %curr = curr(:,t+1:end);
        %add synaptic current:
        B = (tau1^(-1)*tau2 - tau1B^(-1)*tau2B)^(-1)*(tau1^(-1)*d(J+1:2*J,:) - tau1B^(-1)*d(2*J+1:3*J,:));
        C = (tau2^(-1)*tau1 - tau2B^(-1)*tau1B)^(-1)*(tau2^(-1)*d(J+1:2*J,:) - tau2B^(-1)*d(2*J+1:3*J,:));
        D1 = 1/(lambda2-lambda1)*((lambda2*I+A)*(lambda*I+A)*(d(1:J,:)+B) + (lambda1*I+A)*((lambda+lambda2-lambda1)*I+A)*C);
        D2 = 1/(lambda1-lambda2)*((lambda1*I+A)*(lambda*I+A)*(d(1:J,:)+C) + (lambda2*I+A)*((lambda+lambda1-lambda2)*I+A)*B);
        curr(:,Tpast+1:end) = curr(:,Tpast+1:end) + [-D1;tau1*D1;tau1B*D1]*exp(-lambda1*(0:T-Tpast-1)*dt)+[-D2;tau2*D2;tau2B*D2]*exp(-lambda2*(0:T-Tpast-1)*dt);
    else
        Tpast = T;
    end
end
