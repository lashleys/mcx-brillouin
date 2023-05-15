close all 
clear

% addpath(genpath(pwd))
% addpath(genpath('C:\Users\Stephanie Lashley\Documents\matlab-master'))
% addpath(genpath('C:\Users\Stephanie Lashley\Documents\mcx\utils'))

cfg.nphoton=1e7;
cfg.vol=uint8(ones(60,60,60));
cfg.vol(20:40,20:40,10:30)=2;    % add an inclusion
cfg.prop=[0 0 1 1;0.005 1 0 1.37; 0.2 10 0.9 1.37]; % [mua,mus,g,n]
% cfg.prop=[0 0 1 1;0.005 1 0 1.37];

cfg.srctype = 'pattern';
cfg.srcnum = 1;
cfg.srcpattern = zeros(1,60,60);
cfg.srcpattern(1,1,1) = 1;
cfg.srcpattern(1,1,10) = 1;
cfg.srcpattern(1,1,20) = 1;
cfg.srcparam1 = [60 0 0 size(cfg.srcpattern,2)];
cfg.srcparam2 = [0 60 0 size(cfg.srcpattern,3)];
cfg.srcpos=[1 1 1];
cfg.srcdir=[1 0 1]/sqrt(2);
% cfg.srcpos = [1 1 1];
% cfg.srcdir = [1 0 1]/sqrt(2);

cfg.gpuid=1;
% cfg.gpuid='11'; % use two GPUs together
cfg.autopilot=1;
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;
cfg.savedetflag='dspmxvw';
cfg.maxjumpdebug = 1e8;
cfg.detpos = [1 1 55 1; 1 10 55 1; 1 20 55 1];
% cfg.detpos = [55 1 1 1];
%cfg.bc='______111110'; % capture photons existing from all faces except z=z_max

% calculate the flux distribution with the given config
[fluence,detpt,vol,seeds,traj]=mcxlab(cfg);

figure(1)
newtraj=mcxplotphotons_DOT(traj, detpt);title('photon trajectories')
figure(2)
imagesc(squeeze(log(fluence(1).data(:,30,:,1))))
% cfgs(1)=cfg;
% cfgs(2)=cfg;
% cfgs(1).isreflect=0;
% cfgs(2).isreflect=1;
% %cfgs(2).detpos=[30 20 1 1;30 40 1 1;20 30 1 1;40 30 1 1];
% cfgs(2).detpos=[55 0 1 1];
% % calculate the flux and partial path lengths for the two configurations
% [fluxs,detpts,vol,seeds,traj]=mcxlab(cfgs);
% 
% 
% %imagesc(squeeze(log(fluxs(1).data(:,30,:,1)))-squeeze(log(fluxs(2).data(:,30,:,1))));
% %plot3(detpt.p(:,1),detpt.p(:,2),detpt.p(:,3),'r.');
% %view([0 0 -1])
% %view(3)
%%
% cos_alpha = detpt.data(end,:);
% greater_one = find(cos_alpha > 1);
% cos_alpha(cos_alpha>1) = 1;
% brillouin_angle = acos(cos_alpha);
% % brillouin_angle = pi/2 - brillouin_angle;
% % brillouin_angle(brillouin_angle > pi - eps) = []; % removing 0s because they would just add to the input power
% brillouin_angle(brillouin_angle == 0) = [];
% brillouin_angle_sin = sin(brillouin_angle/2);
% c = 3e8;
% lambda = 532e-9;
% f = c/lambda;
% omega = 2*pi*f;
% vs = 1548; % m/s
% % brillouin_shift = 1548.*2.*1.35.*omega.*brillouin_angle_sin./c;
% brillouin_shift = 1548.*2.*1.35.*omega.*brillouin_angle_sin./c;
% 
% figure(1); histogram(brillouin_angle, 'Binwidth', 0.005)
% % figure(2); histogram(sin(brillouin_angle/2),'Binwidth',0.005)
% figure(3); histogram(brillouin_shift/1e9,'Binwidth',0.1) % This plot looks like a mirror image of what I want?
% % figure(3); histogram(brillouin_shift/1e9);
% axis([0 20 0 2500])
% xlabel('Shift (GHz)')
% ylabel('Intensity (a.u)')
% set(gca,'YTick', [])