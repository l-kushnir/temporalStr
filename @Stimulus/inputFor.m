function c = inputFor(this,x,A)
dt = this.dt;
if nargin<3    
    A  = this.A;
end
dx = x(:,2:end)-x(:,1:end-1);
%dx = [dx,dx(end)];
dx = dx/dt;

dx = [dx,dx(:,end)];
c = dx-A*x;
    
end
