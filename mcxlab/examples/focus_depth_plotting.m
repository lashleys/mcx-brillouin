%%
clear;

load('Focus_depth_650micron_1e11_correctshift.mat')
bs650 = brillouin_shift;
ns650 = detpt.nscat;
load('Focus_depth_550micron_1e11_correctshift.mat')
bs550 = brillouin_shift;
ns550 = detpt.nscat;
load('Focus_depth_450micron_1e11_correctshift.mat')
bs450 = brillouin_shift;
ns450 = detpt.nscat;
load('Focus_depth_350micron_1e11_correctshift.mat')
bs350 = brillouin_shift;
ns350 = detpt.nscat;
load('Focus_depth_250micron_1e11_correctshift.mat')
bs250 = brillouin_shift;
ns250 = detpt.nscat;
load('Focus_depth_150micron_1e11_correctshift.mat')
bs150 = brillouin_shift;
ns150 = detpt.nscat;
load('Focus_depth_50micron_1e11_correctshift.mat')
bs50 = brillouin_shift;
ns50 = detpt.nscat;
load('Focus_depth_1050micron_1e11_correctshift.mat')
bs1050 = brillouin_shift;
ns1050 = detpt.nscat;
%% Figure 1c
figure(1);

% 
% [values, edges] = histcounts(brillouin_shift650_correctshift)%,'Normalization','probability','BinWidth',1e6);
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values)
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% hold on
% [values, edges] = histcounts(brillouin_shift550_correctshift);
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values)
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% 
% [values, edges] = histcounts(brillouin_shift450_correctshift);
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values)
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% 
% [values, edges] = histcounts(brillouin_shift350_correctshift);
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values)
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% 
% [values, edges] = histcounts(brillouin_shift250_correctshift);
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values)
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% 
% [values, edges] = histcounts(brillouin_shift150_correctshift);
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values)
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);
% 
% [values, edges] = histcounts(brillouin_shift50,'Normalization','probability');
% centers = (edges(1:end-1)+edges(2:end))/2;
% % plot(centers, values, 'k-')
% photon_nums = interp(values,100);
% shift_vals = interp(edges(1:end-1)+0.025,100);
% plot(shift_vals,photon_nums.*shift_vals);

figure(1)
% subplot(7,1,1)
histogram(bs50/1e9,'BinWidth',0.05)
title('50')
a = axis;
% axis([7 8 0 2500])
% subplot(7,1,2)
figure(3)
histogram(bs150/1e9,'BinWidth',0.05)
title('150')
% y_val = ylim;
axis(a)
ylim([0 7500]);
% subplot(7,1,3)
figure(4)
histogram(bs250/1e9,'BinWidth',0.05)
title('250')
% y_val = ylim;
axis(a)
ylim([0 7500]);
% axis([7e9 8e9 0 0.25])
% subplot(7,1,4)
figure(5)
histogram(bs350/1e9,'BinWidth',0.05)
title('350')
% y_val = ylim;
axis(a)
ylim([0 7500]);
% axis([7e9 8e9 0 0.25])
% subplot(7,1,5)
figure(6)
histogram(bs450/1e9,'BinWidth',0.05)
title('450')
% y_val = ylim;
axis(a)
ylim([0 7500]);
% axis([6.2e9 8e9 0 30])
% subplot(7,1,6)
figure(7)
histogram(bs550/1e9,'BinWidth',0.05)
title('550')
% y_val = ylim;
axis(a)
ylim([0 7500]);
% axis([7e9 8e9 0 0.25])
% subplot(7,1,7)
figure(8)
histogram(bs650/1e9,'BinWidth',0.05)
title('650')
% y_val = ylim;
axis(a)
ylim([0 7500]);
% axis([7e9 8e9 0 0.25])
%,'Bandwidth',0.2

%% Quantifying shift left and tail
ind = 1;
peak_ratio = ones(1,7);
for i = 50:100:650
h = histogram(eval(sprintf('bs%d',i)));%,'BinWidth',0.1);
highest = max(h.Values);
s_highest = max(h.Values(h.Values<max(h.Values)));
peak_ratio(ind) = highest/(length(eval(sprintf('bs%d',i))));
ind = ind + 1;
end
figure;
plot(peak_ratio);

%%
figure;
x = linspace(6.2, 8, 500);
% bw = 0.08;
[f,xi] = ksdensity(bs50/1e9, x, 'bandwidth', 0.052); 
[m,i] = max(f);
shift(1) = xi(i);
h = plot(xi,f);
set(h,'linewidth',2.5)
hold on
[f,xi] = ksdensity(bs150/1e9, x, 'bandwidth', 0.045);
plot(xi,f);
[m,i] = max(f);
shift(2) = xi(i);
[f,xi] = ksdensity(bs250/1e9, x, 'bandwidth', 0.055);
plot(xi,f);
[m,i] = max(f);
shift(3) = xi(i);
[f,xi] = ksdensity(bs350/1e9, x, 'bandwidth', 0.06);
plot(xi,f);
[m,i] = max(f);
shift(4) = xi(i);
[f,xi] = ksdensity(bs450/1e9, x, 'bandwidth', 0.06);
plot(xi,f);
[m,i] = max(f);
shift(5) = xi(i);
[f,xi] = ksdensity(bs550/1e9, x, 'bandwidth', 0.06);
plot(xi,f);
[m,i] = max(f);
shift(6) = xi(i);
[f,xi] = ksdensity(bs650/1e9, x, 'bandwidth', 0.06);
h = plot(xi,f);
[m,i] = max(f);
shift(7) = xi(i);
set(h,'linewidth',2.5)

legend('50','150','250','350','450','550','650','Location','northwest')
axis([7.2 8 0 2.5]')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])

figure;
plot(shift);

%%
H = histogram(bs50);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);
hold on;

H = histogram(bs150);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);

H = histogram(bs250);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);

H = histogram(bs350);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);

H = histogram(bs450);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);

H = histogram(bs550);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);

H = histogram(bs650);
photon_nums = interp(H.Values,100);
photon_nums = photon_nums./max(photon_nums);
shift_vals = interp(H.BinEdges(2:end)+0.01,100);
plot(shift_vals,photon_nums.*shift_vals);

%%
x = linspace(6.2, 8, 500);
bw = 0.035;
[y50, ~, bw1] = ksdensity(bs50 / 1e9, x, 'bandwidth', bw);
[y150, ~, bw2] = ksdensity(bs150 / 1e9, x, 'bandwidth', bw);
[y250, ~, bw3] = ksdensity(bs250 / 1e9, x, 'bandwidth', bw);
[y350, ~, bw4] = ksdensity(bs350 / 1e9, x, 'bandwidth', bw);
[y450, ~, bw5] = ksdensity(bs450 / 1e9, x, 'bandwidth', bw);
[y550, ~, bw6] = ksdensity(bs550 / 1e9, x, 'bandwidth', bw);
[y650, ~, bw7] = ksdensity(bs650 / 1e9, x, 'bandwidth', bw);

y50 = y50 / max(y50);
y150 = y150 / max(y150);
y250 = y250 / max(y250);
y350 = y350 / max(y350);
y450 = y450 / max(y450);
y550 = y550 / max(y550);
y650 = y650 / max(y650);

figure
hold on
plot(x, y50)
plot(x, y150)
plot(x, y250)
plot(x, y350)
plot(x, y450)
plot(x, y550)
plot(x, y650)
legend({'50','150','250','350','450','550','650'})

%% Figure 3a (unpolarized - depolarized) subtract out everything that scattered more than once

sub250 = bs250(nscat250 == 1);
sub450 = bs450(nscat450 == 1);
sub650 = bs650(nscat650 == 1);
depol250 = bs250(nscat250 > 1);
depol450 = bs450(nscat450 > 1);
depol650 = bs650(nscat650 > 1);

%%

% [f,xi] = ksdensity(brillouin_shift50/1e9);
% h = plot(xi,f);
% set(h,'linewidth',2)
% ksdensity(brillouin_shift150/1e9);
figure;
ksdensity(bs250/1e9);
hold on
ksdensity(sub250/1e9);
ksdensity(depol250/1e9);
legend('Unpolarized','Subtracted','Depolarized')
% ksdensity(brillouin_shift350/1e9);
figure;
ksdensity(bs450/1e9);
hold on
ksdensity(sub450/1e9)
ksdensity(depol450/1e9)
legend('Unpolarized','Subtracted','Depolarized')
% ksdensity(brillouin_shift550/1e9);
figure;
[f,xi] = ksdensity(bs650/1e9);
h = plot(xi,f);
hold on
ksdensity(sub650/1e9);
ksdensity(depol650/1e9)
% set(h,'linewidth',2)
legend('Unpolarized','Subtracted','Depolarized')
% axis([7 8 0 3]')

%% 
% figure(1);
bw = 0.05;
figure;
subplot(1,3,1);
histogram(bs250/1e9,'Binwidth',bw);
title('Unpolarized, full shift, 250um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
a = axis;
% figure(2);
subplot(1,3,2);
histogram(depol250/1e9,'Binwidth',bw);
title('Depolarized, only multiple scattering, 250um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)
subplot(1,3,3)
histogram(sub250/1e9,'Binwidth',bw);
title('Subtracted, without multiple scattering, 250um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)

figure;
subplot(1,3,1);
histogram(bs450/1e9,'Binwidth',bw);
title('Unpolarized, full shift, 450um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)
% a = axis;
subplot(1,3,2);
histogram(depol450/1e9,'Binwidth',bw);
title('Depolarized, only multiple scattering, 450um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a);
subplot(1,3,3)
histogram(sub450/1e9,'Binwidth',bw);
title('Subtracted, without multiple scattering, 450um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)

figure;
subplot(1,3,1);
histogram(bs650/1e9,'Binwidth',bw);
title('Unpolarized, full shift, 650um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)
% a = axis;
subplot(1,3,2);
histogram(depol650/1e9,'Binwidth',bw);
title('Depolarized, only multiple scattering, 650um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)
subplot(1,3,3);
histogram(sub650/1e9,'Binwidth',bw);
title('Subtracted, without multiple scattering, 650um')
xlabel('Frequency Shift (GHz)')
ylabel('Intensity (a.u.)')
set(gca,'YTick', [])
axis(a)


%%
figure;







figure;
ksdensity(bs50/1e9)
hold on
ksdensity(bs150/1e9);
ksdensity(bs250/1e9)
%% This is still a histogram for everything except 50?
figure(1)
[val,be] = histcounts(bs50);%,'BinWidth',0.001);
h = plot(be(2:end),val.*be(2:end)./max(val));
set(h,'linewidth',2.5)
hold on
[val,be] = histcounts(bs150);%,'BinWidth',0.001);
plot(be(2:end),val.*be(2:end)./max(val));
[val,be] = histcounts(bs250);%,'BinWidth',0.1);
plot(be(2:end),val.*be(2:end)./max(val));
[val,be] = histcounts(bs350);%,'BinWidth',0.001);
plot(be(2:end),val.*be(2:end)./max(val));
[val,be] = histcounts(bs450);%,'BinWidth',0.001);
plot(be(2:end),val.*be(2:end)./max(val));
[val,be] = histcounts(bs550);%,'BinWidth',0.001);
plot(be(2:end),val.*be(2:end)./max(val));
[val,be] = histcounts(bs650);%,'BinWidth',0.001);
h = plot(be(2:end),val.*be(2:end)./max(val));
set(h,'linewidth',2.5)
legend('50','150','250','350','450','550','650','Location','northwest')
% axis([7.2e9 8e9 0 2.5])
% xlabel('Frequency Shift (GHz)')
% ylabel('Intensity (a.u.)')
% set(gca,'YTick', [])