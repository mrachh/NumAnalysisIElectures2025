n = 600;
ph = 0.5 + rand;
% Define geometry
nosc = 10;
amp = 0.3;
rfun = @(t) 1 + amp*cos(nosc*(t+ph));
rpfun = @(t) -amp*nosc*sin(nosc*(t+ph));
rppfun = @(t) -amp*nosc*nosc*cos(nosc*(t+ph));
xfun = @(t) rfun(t).*cos(t);
yfun = @(t) rfun(t).*sin(t);
dxfun = @(t) rpfun(t).*cos(t) - rfun(t).*sin(t);
dyfun = @(t) rpfun(t).*sin(t) + rfun(t).*cos(t);
d2xfun = @(t) rppfun(t).*cos(t) - 2*rpfun(t).*sin(t) - rfun(t).*cos(t);
d2yfun = @(t) rppfun(t).*sin(t) + 2*rpfun(t).*cos(t) - rfun(t).*sin(t);
dsfun = @(t) sqrt(dxfun(t).^2 + dyfun(t).^2);

% Discretize geometry
alph = rand;
t = alph:2*pi/n:(2*pi + alph - 2*pi/n);
plot(xfun(t), yfun(t), 'k.'); axis equal;

%%
% x_{j}
x = xfun(t); 
y = yfun(t);
% |\gamma'(t_{j})
ds = dsfun(t);
% n(y)
nx = dyfun(t)./dsfun(t);
ny = -dxfun(t)./dsfun(t);

% \kappa
curv = -(dxfun(t).*d2yfun(t) - dyfun(t).*d2xfun(t))./dsfun(t).^3;

kfun = ((x'-x).*nx + (y'-y).*ny)./((x'-x).^2 + (y'-y).^2)/(2*pi);
kfun(1:(n+1):n^2) = curv/(4*pi);

amat = kfun*diag(ds*2*pi/n) - 0.5*eye(n);

%%
x0 = 3.1; y0 = 2.0;
f = log((x-x0).^2  +(y-y0).^2).';
sig = amat\f;

xt = 0.3; yt = 0.2;
uex = log((xt-x0).^2 + (yt - y0).^2);
ktest = ((xt - x).*nx + (yt - y).*ny)./((xt-x).^2 + (yt-y).^2)/2/pi;
ucomp = ktest.*ds*2*pi/n*sig;

norm(ucomp - uex)