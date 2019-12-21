function F = bunch(FF0,Delta)

[N,J] = size(FF0);
dF = [];
for j = 1:J-1;
    dFj = zeros(J,2);
    dFj(j+1,1) = Delta/2;
    dFj(j+1,2) = -Delta/2;
    dF = [dF,dFj];
end
F = [];
for n = 1:N
    F0 = FF0(n,:)';
    TrMatrix = [F0,null(F0')];
    dF = TrMatrix*dF;
    Fn = [F0,F0*ones(1,size(dF,2))+dF];
    F = [F;Fn'];
end

