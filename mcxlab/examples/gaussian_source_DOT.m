clear;
datetime; % datetime to keep track of simulation run time

%% Simulation

cfg.nphoton=1e9; % Number of photons simulated
cfg.vol=uint8(ones(60,60,60)); % size of the medium in x, y, z
cfg.vol(1:60,1:60,15:55) = 2; % inclusion of different material type
cfg.gpuid=1;
cfg.autopilot=1;

cfg.prop=[0 0 1 1;1.35 0.02 0 1.33;0.4 1.8 0.7 1.33]; % [mua,mus,g,n] water
% and a fake material
cfg.issaveexit=13;
cfg.unitinmm=0.5; % dimension of each voxel in mm

cfg.tstart=0;
cfg.seed=99999;

cfg.srcpos=[30 30 0]; % x, y, z
cfg.srcdir = [0 0 1 1]; % aimed in z direction focused at depth of 1 voxel
cfg.tend=5e-11;
cfg.tstep=5e-11;
cfg.wavelength = 532e-9; 

% a gaussian source
cfg.srctype='gaussian';
cfg.srcparam1=[1 0 0 0];
cfg.srcparam2=[0 0 0 0];
det_pos = 28; % Separating source and detector
cfg.detpos = [30 det_pos 0 1]; % x, y, z, radius

cfg.savedetflag='dpsx'; % DO NOT CHANGE
cfg.debuglevel = 'M';
[flux,detpt,vol,seed,trajectory,~,all_brillouin_data]=mcxlab(cfg); 
% detpt has all detected photon data
% all_brillouin_data has the Brillouin shift of every simulated photon, not
% just detected photons
datetime % datetime to track simulation run time%% This is with multiple materials

%% Calculating shift
cos_alpha = detpt.data(end-1,:); % Brillouin shift data is stored in 
% detpt.data(end-1,:), this is the cosine of the scattering angle
ind = find(cos_alpha == 5); % default value of cos_alpha is 5
cos_alpha(ind) = []; % filter out any photons without Brillouin data
cos_alpha(cos_alpha>1) = 1; % fix floating point error
cos_alpha(cos_alpha<-1) = -1; % fix floating point error
brillouin_angle = acos(cos_alpha); % take acos to get actual scattering angle

material = detpt.data(end,:); % The material type (1 or 2) is stored in
% detpt.data(end,:)
material(ind) = []; % filter out material at indices where the photons did
% not have a Brillouin event
fake_ind = find(material == 2); % indices for fake material (or material 2)
water_ind = find(material == 1); % indices for material 1 (water)
brillouin_shift = ones(1,length(material)); % make vector to store shifts

n_water = 1.33; % refractive index of water
n_fake = 1.33; % refractive index of fake material
lambda = cfg.wavelength;
vs_water = 1498; % speed of sound of water
% speed of sound of glycerin (chose something with a large value)
vs_fake = 1904; %http://koski.ucdavis.edu/BRILLOUIN/glycerin/glycerin.html
% fake shift should be just less than 10.5GHz if backscattering
brillouin_shift(fake_ind) = 2*n_fake.*vs_fake.*sin(brillouin_angle(fake_ind)./2)./lambda;
brillouin_shift(water_ind) = 2*n_water.*vs_water.*sin(brillouin_angle(water_ind)./2)./lambda;

figure(det_pos); H = histogram(brillouin_shift/1e9,'Binwidth',0.05); 
title('Brillouin Shift')
xlabel('Shift (GHz)')
ylabel('Number of photons)')

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
plot3(30,30,0,'ro') % Source
plot3(30,det_pos,0,'bo') % Detector
view(3)