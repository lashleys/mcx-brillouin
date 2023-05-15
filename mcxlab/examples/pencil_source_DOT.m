%% test group 4
addpath('/home/stephanie/Documents/mcx-brillouin/');
clear cfg;
close all;
% cfg.nphoton=1e10;
cfg.nphoton=1e10;
cfg.vol=uint8(ones(60,60,60));
cfg.vol(1:60,1:60,10:40) = 2;
cfg.gpuid=1;
cfg.autopilot=1;
n = 1.35;
vs = 1548; % milk
% cfg.prop=[0 0 1 1;0.56 1.68 0.9 1.35]; % [mua,mus,g,n] milk
cfg.prop=[0 0 1 1;1.35 0.02 0 1.33;0.56 1.68 0.9 1.35]; % [mua,mus,g,n] milk

cfg.issaveexit=1;
cfg.unitinmm=0.5;
% want 10x10x10mm so this should be changed but doesn't matter right now

% n = 1.33;
% vs = 1498; % water
% cfg.prop=[0 0 1 1;0.0135 0.0002 0 1.33]; % [mua,mus,g,n] water
% cfg.prop=[0 0 1 1;1.35 0.02 0 1.33]; % [mua,mus,g,n] water
cfg.tstart=0;
cfg.seed=99999;

cfg.srcpos=[30 30 0];
% cfg.srcdir=[0 -sind(30) cosd(30)]; % 
cfg.srcdir=[0 1/sqrt(2) 1/sqrt(2)];
cfg.tend=5e-11;
cfg.tstep=5e-11;
cfg.wavelength = 532e-9; 

% a pencil source
cfg.srctype='pencil';
% cfg.bc='______111110';
cfg.detpos = [30 30 0 1]; % x, y, z, radius
% cfg.detpos = [30 30 0 1];
cfg.savedetflag='dpsx';
cfg.debuglevel = 'M';
[flux,detpt,vol,seed,trajectory,~,all_data]=mcxlab(cfg); 


%% This is without multiple materials
% cos_alpha = detpt.data(end-1,:);
% cos_alpha(cos_alpha==5) = [];
% cos_alpha(cos_alpha>1) = 1;
% cos_alpha(cos_alpha<-1) = -1;
% brillouin_angle = acos(cos_alpha);
% 
% n = 1.35;
% % n = 1.33
% lambda = 532e-9;
% vs = 1548; % m/s
% % vs = 1498;
% brillouin_shift = 2*n.*vs.*sin(brillouin_angle./2)./lambda;
% 
% figure(3); H = histogram(brillouin_shift/1e9,'Binwidth',0.05); 
% % axis([7 8 0 6000]);
% title('Brillouin Shift')
% xlabel('Shift (GHz)')
% ylabel('Intensity (a.u)')

%%  This is with multiple materials
material = detpt.data(end,:);
water_ind = find(material == 1);
milk_ind = find(material == 2);
brillouin_shift = ones(1,length(material));
%% This is with multiple materials
cos_alpha = detpt.data(end-1,:);
cos_alpha(cos_alpha==5) = [];
cos_alpha(cos_alpha>1) = 1;
cos_alpha(cos_alpha<-1) = -1;
brillouin_angle = acos(cos_alpha);

n_milk = 1.35;
n_water = 1.33;
lambda = 532e-9;
vs_milk = 1548; % m/s
vs_water = 1498;
% brillouin_shift = 2*n.*vs.*sin(brillouin_angle./2)./lambda;
brillouin_shift(water_ind) = 2*n_water.*vs_water.*sin(brillouin_angle(water_ind)./2)./lambda;
brillouin_shift(milk_ind) = 2*n_milk.*vs_milk.*sin(brillouin_angle(milk_ind)./2)./lambda;

figure(3); H = histogram(brillouin_shift/1e9,'Binwidth',0.05); 
% axis([7 8 0 6000]);
title('Brillouin Shift')
xlabel('Shift (GHz)')
ylabel('Intensity (a.u)')
% set(gca,'YTick', [])
% E = h_bar*omega
% plot brillouin shift vs energy (omega) * number of photons with that
% energy (omega)
% figure;
% photon_nums = interp(H.Values,100);
% shift_vals = interp(H.BinEdges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% 
% figure;
% all_brillouin_data = all_data.data;
% all_brillouin_data(all_brillouin_data==5) = [];
% all_brillouin_data(all_brillouin_data>1) = 1;
% all_brillouin_data(all_brillouin_data<-1) = -1;
% all_brillouin_angles = acos(all_brillouin_data);
% histogram(all_brillouin_angles);
% figure;
% all_brillouin_shift = 2*n.*vs.*sin(all_brillouin_angles./2)./lambda;
% histogram(all_brillouin_shift)
% 
%% Plotting
figure(1)
plot3(detpt.p(:,1),detpt.p(:,2),detpt.p(:,3),'r.');
% view([0 0 -1])
title('a uniform disk source');
xlabel('x');ylabel('y');zlabel('z')
hold on
% whole material
plot3(0:60,zeros(1,61),zeros(1,61),'r');
plot3(zeros(1,61),0:60,zeros(1,61),'r');
plot3(0:60,ones(1,61)*60,zeros(1,61),'r');
plot3(ones(1,61).*60,0:60,zeros(1,61),'r');

plot3(zeros(1,61),zeros(1,61),0:60,'r');
plot3(ones(1,61)*60,zeros(1,61),0:60,'r');
plot3(0:60,zeros(1,61),ones(1,61)*60,'r');

plot3(zeros(1,61),ones(1,61)*60,0:60,'r');
plot3(ones(1,61)*60,ones(1,61)*60,0:60,'r');
plot3(0:60,ones(1,61)*60,ones(1,61)*60,'r');

plot3(zeros(1,61),0:60,ones(1,61)*60,'r');
plot3(ones(1,61)*60,0:60,ones(1,61)*60,'r');
plot3(0:60,ones(1,61)*60,zeros(1,61),'r');

plot3(ones(1,61)*60,0:60,zeros(1,61),'r');
plot3(ones(1,61)*60,0:60,ones(1,61)*60,'r')
plot3(ones(1,61)*60,0:60,ones(1,61)*60,'r');

% inclusion
plot3(0:60,zeros(1,61),ones(1,61)*10,'b');
plot3(zeros(1,61),0:60,ones(1,61)*10,'b');
plot3(0:60,ones(1,61)*60,ones(1,61)*10,'b');
plot3(ones(1,61).*60,0:60,ones(1,61)*10,'b');

plot3(zeros(1,61),zeros(1,61),linspace(10,40,61),'b');
plot3(ones(1,61)*60,zeros(1,61),linspace(10,40,61),'b');
plot3(0:60,zeros(1,61),ones(1,61)*40,'b');

plot3(zeros(1,61),ones(1,61)*60,linspace(10,40,61),'b');
plot3(ones(1,61)*60,ones(1,61)*60,linspace(10,40,61),'b');
plot3(0:60,ones(1,61)*60,ones(1,61)*40,'b');

plot3(zeros(1,61),0:60,ones(1,61)*40,'b');
plot3(ones(1,61)*60,0:60,ones(1,61)*40,'b');
% source and detector

plot3(30,25,0,'k.');
plot3(30,35,0,'k.');
