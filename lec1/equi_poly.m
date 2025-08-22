f = @(x) exp(x) + sin(5*x);
% f = @(x) exp(abs(sin(x)).^3);


ttest = 0:0.003:1;

% Plot the function
figure(1)
clf
plot(ttest, f(ttest), 'k.');

%%
% Compute trig interpolant for fixed n and test accuracy
n = 8;
ts = (0:n)/(n+1);
vmat = ts(:).^(0:n);
[u,s,v] = svd(vmat);
ft = f(ts(:));
fcoefs = v*(inv(s)*(u'*ft));

% Plot the interpolant
amat = ttest(:).^(0:n);
ftest = amat*fcoefs(:); % Evaluate interpolant at test points
hold on;
plot(ttest,ftest,'r.')

% Evaluate interpolant using barycentric formula
wmat = log(abs(ts - ts.'));
wmat(isinf(wmat)) = 0;
bwts = 1./((-1).^(0:n).*exp(sum(wmat)));

ftest2 = zeros(size(ftest));

[tuse, ia] = setdiff(ttest, ts);
[t2, ia2, ib2] = intersect(ttest, ts);
ftest2(ia2) = ft(ib2);
aa = (ft.*bwts(:))./(tuse - ts.');
bb = bwts(:)./(tuse - ts.');
ftest2(ia) = sum(aa)./sum(bb);


% Test accuracy of interpolant

fex = f(ttest(:)); % Compute exact function values at test points
err_est = norm(ftest - fex);

% alternate test, using just the coefficients when exact function
% is not known
err_est2 = norm(fcoefs(end-1:end));

figure(2)
clf
semilogy(0:n,abs(fcoefs),'k.','MarkerSize',15)


%% Convergence plots

nvec = 1:1:16;
errs = zeros(size(nvec));
for i=1:length(nvec)
    n = nvec(i);
    ts = (0:n)/(n+1);
    vmat = ts(:).^(0:n);
    [u,s,v] = svd(vmat);
    ft = f(ts(:));
    fcoefs = v*(inv(s)*(u'*ft));
    amat = ttest(:).^(0:n);
    ftest = amat*fcoefs(:); % Evaluate interpolant at test points

    errs(i) = norm(ftest - fex);
end

figure(3)
clf
semilogy(nvec, errs, 'k.', 'MarkerSize',15);