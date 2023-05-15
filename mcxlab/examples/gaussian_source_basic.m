clear;
datetime % datetime to keep track of simulation run time

%% Simulation

cfg.nphoton=1e8; % Number of photons simulated
cfg.vol=uint8(ones(60,60,60)); % size of the medium in x, y, z
cfg.gpuid=1;
cfg.autopilot=1;

% Material of the medium
% cfg.prop=[0 0 1 1;0.56 1.68 0.9 1.35]; % [mua,mus,g,n] milk
cfg.prop=[0 0 1 1;1.35 0.02 0 1.33]; % [mua,mus,g,n] water
cfg.issaveexit=13;
cfg.unitinmm=0.05; % dimension of each voxel in mm

cfg.tstart=0;
cfg.seed=99999;

cfg.srcpos=[30 30 0]; % x, y, z
focus_depth = 1;
% Aimed in z direction
cfg.srcdir=[0 0 1 focus_depth]; % focus depth is depth of gaussian focus
cfg.tend=5e-11;
cfg.tstep=5e-11;
cfg.wavelength = 532e-9; 

% a gaussian source
cfg.srctype='gaussian';
cfg.srcparam1=[1 0 0 0];
cfg.srcparam2=[0 0 0 0];
cfg.detpos = [30 30 0 1]; % x, y, z, radius

cfg.savedetflag='dpsx'; % DO NOT CHANGE
cfg.debuglevel = 'M';
[flux,detpt,vol,seed,trajectory,~,all_brillouin_data]=mcxlab(cfg); % MCX code
% detpt has all detected photon data
% all_brillouin_data has the Brillouin shift of every simulated photon, not
% just detected photons
datetime % datetime to track simulation run time

%% Calculating shift
cos_alpha = detpt.data(end-1,:); % Brillouin shift data is stored in 
% detpt.data(end-1,:), this is the cosine of the scattering angle
cos_alpha(cos_alpha==5) = []; % the default value of cos_alpha is 5 so
% a value of 5 means it didn't have a Brillouin event and should be 
% filtered out
cos_alpha(cos_alpha>1) = 1; % This removes floating point error where values
% could be 1.00000001
cos_alpha(cos_alpha<-1) = -1; % removes floating point error for -1
brillouin_angle = acos(cos_alpha); % take acos to get actual scattering angle

lambda = cfg.wavelength;
% vs = 1548; % milk speed of sound
% n = 1.35; % milk refractive index
vs = 1498; % water speed of sound
n = 1.33; % water refractive index
brillouin_shift = 2*n.*vs.*sin(brillouin_angle./2)./lambda; % calculate shift

% Plot histogram of results
figure; H = histogram(brillouin_shift/1e9,'Binwidth',0.05); 

title('Brillouin Shift')
xlabel('Shift (GHz)')
ylabel('Number of photons')
% Save off figure as .fig and data as .mat (comment unless an important
% run)
saveas(gcf,sprintf('change_figure_name.fig'))
save(sprintf('change_mat_name.mat'))


%% Plot of medium with location of source/detector for 
% visualization
figure;
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
xlabel('x'); ylabel('y'); zlabel('z');
plot3(30,30,0,'ro') % Source/detector location
view(3)