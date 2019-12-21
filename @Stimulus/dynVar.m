function x = dynVar(this,c,A)

[N,T] = size(c);
if all(size(A)==[1,1])
    A = A*eye(N);
end

dt = this.dt;
tt = 0:dt:(T*dt);

[V,D] = eig(A);
Vin = V^(-1);
cont = zeros(N,T);
contk = zeros(N,T);

for k = 1:N
    for j = 1:N
        ker = [0,exp(log(1+D(k,k)*dt)*tt/dt)];
        p = find(abs(ker)<0.0001);
        if length(p)>2
         ker = ker(1:p(3));
        end
        contjk = dt*conv(ker,c(j,:));
        contjk = contjk(1:T);
        contk(j,:) = contjk*Vin(k,j);
    end 
    cont(k,:) = sum(contk,1);
end
x = real(V*cont);
end