theta = linspace(0,pi,100000); 
g = 0.9; 
cdf = ((1-g.^2)./(2*g)).*(1./sqrt(1 + g.^2 - 2.*g.*cos(theta)) - 1./(1 + g));
figure;
plot(theta,1-cdf);
hold on

g = 0.1; 
cdf = ((1-g.^2)./(2*g)).*(1./sqrt(1 + g.^2 - 2.*g.*cos(theta)) - 1./(1 + g));
plot(theta,1-cdf);

g = -0.9; 
cdf = ((1-g.^2)./(2*g)).*(1./sqrt(1 + g.^2 - 2.*g.*cos(theta)) - 1./(1 + g));
plot(theta,1-cdf);

g = 0;
cdf = (cos(theta) + 1)./2;
plot(theta,1-cdf)

% axis([pi/4 pi 0 1])
legend('g=0.9', 'g=0.1','g=-0.9','g=0')