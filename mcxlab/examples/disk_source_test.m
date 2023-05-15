%% test group 4
addpath('/home/stephanie/Documents/mcx-brillouin/');
clear cfg;
close all;
cfg.nphoton=1e7;
cfg.vol=uint8(ones(60,60,60));
cfg.gpuid=1;
cfg.autopilot=1;
cfg.prop=[0 0 1 1;0.56 1.68 0.9 1.35]; % [mua,mus,g,n]
% cfg.prop=[0 0 1 1; 0.005 1 1 1.37];
% this absorbed ~80% vs default that absorbed ~1.5% 
cfg.tstart=0;
cfg.seed=99999;

cfg.srcpos=[30 30 -10];
cfg.srcdir=[0 0 1 0];
cfg.tend=5e-11;
cfg.tstep=5e-11;
cfg.wavelength = 532e-9; 

% a uniform disk source
cfg.srctype='disk';
cfg.srcparam1=[20 0 0 0];
cfg.srcparam2=[0 0 0 0];
cfg.bc='______111110';  % capture photons existing from all faces except z=z_max
% cfg.bc='______000000'; 
% Am I not actually backscattering? No I am, verified
cfg.savedetflag='dpsx';
cfg.debuglevel = 'M';
[flux,detpt,vol,seed,trajectory,~,all_brillouin_angles]=mcxlab(cfg); % Brillouin is actually just theta angle right now
all_brillouin_angles.data(all_brillouin_angles.data == 0) = [];
all_brillouin_angles.data(all_brillouin_angles.data == 5) = [];
all_brillouin_angles.data(all_brillouin_angles.data > 1) = 1;
figure;
histogram(acos(all_brillouin_angles.data),'Binwidth',0.005);
% histogram((all_brillouin_angles.data))
% brillouin_angle = acos(all_brillouin_angles.data);
% cos_alpha = detpt.data(end,:);
% brillouin_angle = acosd(cos_alpha);
% figure(1);
% histogram(brillouin_angle);
%%
cos_alpha = detpt.data(end,:);
greater_one = find(cos_alpha > 1);
cos_alpha(cos_alpha>1) = 1;
brillouin_angle = acos(cos_alpha);

c = 3e8;
lambda = 532e-9;
f = c/lambda;
omega = 2*pi*f; % Should I not be doing this?
vs = 1548; % m/s
brillouin_shift = 1548.*2.*1.35.*omega.*sin(brillouin_angle./2)./c;

figure(2); histogram(brillouin_angle, 'Binwidth', 0.01)
title('Brillouin Angle')
figure(6); histogram(sin(brillouin_angle./2),'Binwidth',0.005)
title('sin(brillouin angle/2)')
figure(3); histogram(brillouin_shift/1e9,'Binwidth',0.05) % This plot looks like a mirror image of what I want?
title('Brillouin Shift')
xlabel('Shift (GHz)')
ylabel('Intensity (a.u)')
set(gca,'YTick', [])

%% Plotting
figure(1)
plot3(detpt.p(:,1),detpt.p(:,2),detpt.p(:,3),'r.');
% view([0 0 -1])
title('a uniform disk source');
xlabel('x');ylabel('y');zlabel('z')
hold on
plot3(0:60,zeros(1,61),zeros(1,61));
plot3(zeros(1,61),0:60,zeros(1,61));
plot3(0:60,ones(1,61)*60,zeros(1,61));
plot3(ones(1,61).*60,0:60,zeros(1,61));

plot3(zeros(1,61),zeros(1,61),0:60);
plot3(ones(1,61)*60,zeros(1,61),0:60);
plot3(0:60,zeros(1,61),ones(1,61)*60);

plot3(zeros(1,61),ones(1,61)*60,0:60);
plot3(ones(1,61)*60,ones(1,61)*60,0:60);
plot3(0:60,ones(1,61)*60,ones(1,61)*60);

plot3(zeros(1,61),0:60,ones(1,61)*60);
plot3(ones(1,61)*60,0:60,ones(1,61)*60);