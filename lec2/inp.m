
nch = 80;

narms = 3;
amp = 0.2;
chnkr = chunkerfuncuni(@(t) starfish(t,narms,amp),nch); 

figure(1)
clf
plot(chnkr, 'k.'); axis equal

%% single layer representation

S = kernel('laplace', 'single');
Sp = kernel('laplace', 'sprime');

% Construct sprime
A = chunkermat(chnkr, Sp);
A = A + eye(chnkr.npt)/2;

%%
s = svd(A);

A = A + onesmat(chnkr);
%%
s2 = svd(A);

src0 = [3.1; 2.0];
srcinfo = [];
srcinfo.r = src0;
% compute f = \partial G/\partial n(x,src0), x = \partial \Omega
f = Sp.eval(srcinfo, chnkr);

sigma = A\f;

%%
targinfo = [];
targinfo.r = [-0.3, 0.2; -0.6 0.5];

ucomp = chunkerkerneval(chnkr, S, sigma, targinfo);

%%
uex = S.eval(srcinfo, targinfo);