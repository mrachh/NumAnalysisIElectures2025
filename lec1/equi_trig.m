f = @(x) exp(cos(x)) + sin(sin(x));

ttest = 0:0.003:2*pi;

% Plot the function
figure(1)
clf
plot(ttest, f(ttest), 'k.');

%%
% Compute trig interpolant for fixed n and test accuracy
n = 60;
ts = (0:2*n)*2*pi/(2*n+1);
ft = f(ts);
fcoefs = fftshift(fft(ft));

% Test accuracy of interpolant

amat = exp(1j*ttest.'*(-n:n));
ftest = amat*fcoefs(:)/(2*n+1); % Evaluate interpolant at test points
fex = f(ttest(:)); % Compute exact function values at test points
err_est = norm(ftest - fex);

% alternate test, using just the coefficients when exact function
% is not known
err_est2 = norm([fcoefs(1:2) fcoefs(end-1:end)]);

figure(2)
clf
semilogy(-n:n,abs(fcoefs),'k.','MarkerSize',15)


%% Convergence plots

nvec = 5:5:100;
errs = zeros(size(nvec));
for i=1:length(nvec)
    n = nvec(i);
    ts = (0:2*n)*2*pi/(2*n+1);
    ft = f(ts);
    fcoefs = fftshift(fft(ft));
    errs(i) = norm([fcoefs(1:2) fcoefs(end-1:end)]);
end

figure(3)
clf
semilogy(nvec, errs, 'k.', 'MarkerSize',15);