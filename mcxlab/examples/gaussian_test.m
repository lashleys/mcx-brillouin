% for focus_depth = 5:2:13
clear;
focus_depth = 1;
sprintf('\nRunning focus depth %d\n',focus_depth)
addpath('/home/stephanie/Documents/mcx-brillouin/');
clear cfg;
% close all;

cfg.nphoton=1e8;
cfg.vol=uint8(ones(200,200,200));
cfg.gpuid=1;
cfg.autopilot=1;
n = 1.35;
vs = 1548; % milk

cfg.prop=[0 0 1 1;0.56 1.68 0.9 1.35]; % [mua,mus,g,n] milk
cfg.issaveexit=1;
cfg.unitinmm=0.05;

cfg.tstart=0;
cfg.seed=99999;

cfg.srcpos=[100 100 0];
cfg.srcdir=[0 0 1 focus_depth];
cfg.tend=5e-11;
cfg.tstep=5e-11;
cfg.wavelength = 532e-9; 

% a gaussian source
cfg.srctype='gaussian';
cfg.srcparam1=[1 0 0 0];
cfg.srcparam2=[0 0 0 0];
% cfg.bc='______111110';
cfg.detpos = [100 100 0 1]; % x, y, z, radius

cfg.savedetflag='dpsx';
cfg.debuglevel = 'M';
[flux,detpt,vol,seed,trajectory,~,all_data]=mcxlab(cfg); 

%%
cos_alpha = detpt.data(end-1,:);
cos_alpha(cos_alpha==5) = [];
cos_alpha(cos_alpha>1) = 1;
cos_alpha(cos_alpha<-1) = -1;
brillouin_angle = acos(cos_alpha);


lambda = 532e-9;

brillouin_shift = 2*n.*vs.*sin(brillouin_angle./2)./lambda;

figure(focus_depth); H = histogram(brillouin_shift/1e9,'Binwidth',0.05); 

title('Brillouin Shift')
xlabel('Shift (GHz)')
ylabel('Intensity (a.u)')
% saveas(gcf,sprintf('Focus_depth_%dmicron_1e11_correctshift.fig',focus_depth*50))
% save(sprintf('Focus_depth_%dmicron_1e11_correctshift.mat',focus_depth*50))
% end

%%
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
%% Plotting
% figure(1)
% plot3(detpt.p(:,1),detpt.p(:,2),detpt.p(:,3),'r.');
% % view([0 0 -1])
% title('a uniform disk source');
% xlabel('x');ylabel('y');zlabel('z')
% hold on
% plot3(0:60,zeros(1,61),zeros(1,61));
% plot3(zeros(1,61),0:60,zeros(1,61));
% plot3(0:60,ones(1,61)*60,zeros(1,61));
% plot3(ones(1,61).*60,0:60,zeros(1,61));
% 
% plot3(zeros(1,61),zeros(1,61),0:60);
% plot3(ones(1,61)*60,zeros(1,61),0:60);
% plot3(0:60,zeros(1,61),ones(1,61)*60);
% 
% plot3(zeros(1,61),ones(1,61)*60,0:60);
% plot3(ones(1,61)*60,ones(1,61)*60,0:60);
% plot3(0:60,ones(1,61)*60,ones(1,61)*60);
% 
% plot3(zeros(1,61),0:60,ones(1,61)*60);
% plot3(ones(1,61)*60,0:60,ones(1,61)*60);