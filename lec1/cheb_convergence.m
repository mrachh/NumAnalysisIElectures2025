eta = 0.01;
ms = 'x';
f = @(x) eta.^2./(x.^2 + eta^2);
fex = chebfun(f); 
nvec = 1:5:length(fex);
rho = eta + sqrt(eta^2+1);
fcoeffs = chebcoeffs(fex);
figure(1)
semilogy(nvec,abs(fcoeffs(1:5:end)),['k' ms],'MarkerSize',5);
hold on;
semilogy(nvec,fcoeffs(1)./rho.^(nvec-1),['r' ms],'MarkerSize',5);
