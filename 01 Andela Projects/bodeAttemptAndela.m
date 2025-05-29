i = sqrt(-1);
c=5000;
k=70000;
m=1200;
zeta = c/(2*sqrt(k*m));
r = 0:0.05:3;
s = 2*zeta*r*i;
Y = (1+s)./(1 - r.^2 + s);

axis, axis ('normal')
subplot(2,1,1), plot(r, abs(Y))
subplot(2, 1,  2), plot(r, 180*angle(Y)/pi)