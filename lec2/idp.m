
nch = 80;

narms = 3;
amp = 0.2;
chnkr = chunkerfuncuni(@(t) starfish(t,narms,amp),nch); 

figure(1)
clf
plot(chnkr, 'k.'); axis equal

%% single layer representation

S = kernel('laplace', 'single');
opts = [];
opts.l2scale = true;
A = chunkermat(chnkr, S, opts);
%%
[evs] = eig(A);
%%
figure(2)
semilogy(sort(abs(evs)), 'g.'); 

%% double layer representation
D = kernel('lap', 'd');
B = chunkermat(chnkr, D, opts);
B = B - eye(chnkr.npt)/2;
%%
[evs2] = eig(B);
figure(3)
semilogy(abs(evs2), 'g.');

%% Testing the codes
% matrices regenerated with square-root scaling
A = chunkermat(chnkr, S);
B = chunkermat(chnkr, D); B = B - eye(chnkr.npt)/2;
%%
srcinfo = [];
src_out = [3.5;1.2];
srcinfo.r = src_out;
ubdry = S.eval(srcinfo, chnkr);
%%
sig_s = A\ubdry;
sig_d = B\ubdry;

%%
targinfo = [];
targinfo.r = [0.1;0.3];
uex = S.eval(srcinfo, targinfo);
%%
ucomp_s = chunkerkerneval(chnkr, S, sig_s, targinfo, opts);
ucomp_d = chunkerkerneval(chnkr, D, sig_d, targinfo, opts);


