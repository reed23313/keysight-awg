%% demo with two amplfiers 405
load('Z:\SC\Measurements\SPG716\demo\demo_right_value\demo\SPG716_demo_demo_demo_right_value_2022-04-01 18-00-48.mat')
% i = 400;
i = 100;
gain = 1;
figure
subplot(2,1,1)
plot(C1x(i,:)*1e6+5, C1y(i,:)/gain, 'DisplayName', ' SNSPD')
ylabel('voltage (V)')
% xlim([0 8])
set(gca,'XTickLabel',[])
ylim([-0.2 1])
legend
subplot(2,1,2)
plot(C2x(i,:)*1e6+5, (C2y(i,:)+0.4)/gain, 'DisplayName', ' 1st nT')
hold on
plot(C3x(i,:)*1e6+5, (C3y(i,:)+0.2)/gain, 'DisplayName', ' 2nd nT')
% y = C4y(i,:);
plot(C4x(i,:)*1e6+5, C4y(i,:)/gain, 'DisplayName', ' 3rd nT')
xlabel('time (\mus)')
ylabel('voltage (V)')
% xlim([0 8])
ylim([-0.08 0.5])
legend




%% snspd counter light on perfect 405nm 
load('Z:\SC\Measurements\SPG716\SNSPD_200_150n_counter_2bit_WORKING\snspd_counter_2bit\snspd_2000sq_150nm_2bitcounter_working_1000\SPG716_snspd_2000sq_150nm_2bitcounter_working_1000_SNSPD_200_150n_counter_2bit_WORKING_snspd_counter_2bit_2022-03-23 00-33-56.mat')
% 
% for i = 300:10:1000
%     figure(1)
%     i
%     plot(C1x(i,:), C1y(i,:)*0.2, 'DisplayName', ' SNSPD')
%     hold on
%     plot(C2x(i,:), C2y(i,:)-0.005, 'DisplayName', ' 1st nT')
%     plot(C3x(i,:), C3y(i,:)-0.012, 'DisplayName', ' 2nd nT')
%     plot(C4x(i,:), C4y(i,:)-0.018, 'DisplayName', ' 3nd nT')
%     hold off
%     xlabel('time (s)')
%     ylabel('voltage (V)')
%     legend
%     pause(10)
% end
i = 300;
gain = 31.6;
figure
subplot(2,1,1)
plot(C1x(i,:)*1e6+6, C1y(i,:)*1e3/gain, 'DisplayName', ' SNSPD')
ylabel('voltage (mV)')
xlim([0 11])
set(gca,'XTickLabel',[])
ylim([-0.25 1.5])
legend
subplot(2,1,2)
y = bandstop(C2y(i,:),[0.95e6 1.1e6],1/(C2x(i,2)-C2x(i,1)));
plot(C2x(i,:)*1e6+6, (y+0.017)*1e3/gain, 'DisplayName', ' 1st nT')
hold on
plot(C3x(i,:)*1e6+6, (C3y(i,:)+0.008)*1e3/gain, 'DisplayName', ' 2nd nT')
% y = C4y(i,:);
y = bandstop(C4y(i,:),[0.95e6 1.1e6],1/(C4x(i,2)-C4x(i,1)));
plot(C4x(i,:)*1e6+6, y*1e3/19.95, 'DisplayName', ' 3rd nT')
xlabel('time (\mus)')
ylabel('voltage (mV)')
xlim([0 11])
ylim([-0.1 0.8])
legend


%% snspd counter light on perfect 1550nm 
load('Z:\SC\Measurements\SPG716\SNSPD_200_150n_counter_2bit_WORKING2\snspd_counter_2bit_2light\working_light_max_1550_pulsed_long\SPG716_working_light_max_1550_pulsed_long_SNSPD_200_150n_counter_2bit_WORKING2_snspd_counter_2bit_2light_2022-03-23 01-36-38.mat')

% for i = 10:10:1000
%     figure(1)
%     plot(C1x(i,:), C1y(i,:)*0.2, 'DisplayName', ' SNSPD')
%     hold on
%     plot(C2x(i,:), C2y(i,:)-0.005, 'DisplayName', ' 1st nT')
%     plot(C3x(i,:), C3y(i,:)-0.012, 'DisplayName', ' 2nd nT')
%     plot(C4x(i,:), C4y(i,:)-0.018, 'DisplayName', ' 3nd nT')
%     hold off
%     xlabel('time (s)')
%     ylabel('voltage (V)')
%     legend
%     pause(10)
% end
i = 80;
gain = 31.6;
figure
subplot(2,1,1)
plot(C1x(i,:)*1e6+18, C1y(i,:)*1e3/gain, 'DisplayName', ' SNSPD')
ylabel('voltage (mV)')
xlim([0 22])
set(gca,'XTickLabel',[])
ylim([-0.25 1.5])
legend
subplot(2,1,2)
y = bandstop(C2y(i,:),[0.95e6 1.1e6],1/(C2x(i,2)-C2x(i,1)));
plot(C2x(i,:)*1e6+18, (y+0.012)*1e3/gain, 'DisplayName', ' 1st nT')
hold on
plot(C3x(i,:)*1e6+18, (C3y(i,:)+0.006)*1e3/gain, 'DisplayName', ' 2nd nT')
% y = C4y(i,:);
y = bandstop(C4y(i,:),[0.95e6 1.1e6],1/(C4x(i,2)-C4x(i,1)));
plot(C4x(i,:)*1e6+18, y*1e3/19.95, 'DisplayName', ' 3rd nT')
xlabel('time (\mus)')
ylabel('voltage (mV)')
xlim([0 22])
ylim([-0.1 0.8])
legend

%% snspd 1 bit counter light on perfect 1550nm 
load('Z:\SC\Measurements\SPG716\SNSPD_5000_100n_counter_1bit\snspd_counter_1bit_1550\working_light_max_1550_pulsed\SPG716_working_light_max_1550_pulsed_SNSPD_5000_100n_counter_1bit_snspd_counter_1bit_1550_2022-03-23 02-19-59.mat')

for i = 10:10:1000
    figure(1)
    plot(C1x(i,:), C1y(i,:)*0.2, 'DisplayName', ' SNSPD')
    hold on
    plot(C2x(i,:), C2y(i,:)-0.005, 'DisplayName', ' 1st nT')
    plot(C3x(i,:), C3y(i,:)-0.012, 'DisplayName', ' 2nd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(10)
end


%% sweep 2bit counter snspd perfect fine
clear all
myFolder = 'Z:\SC\Measurements\SPG716\SNSPD_200_150n_counter_2bit_WORKING\snspd_counter_2bit\prob_snspd\'; % Define your working folder
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);

right_prob2 = zeros(1, length(matFiles));

min_peak_height = 0.002;
time_tol = 0.1/10e6;

baseFileName = matFiles.name;
fullFileName = fullfile(myFolder, baseFileName);

base = 2;
right_prob2 = compute_prob2bit_snspd(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'For V_c1 = %.3f, V_c2 = %.3f, V_g = %.3f, freq = %.3f\n', v_c1(k), v_c2(k), v_g(k), freq(k));
fprintf(1, 'the probablity to work in base 2 is %f\n', right_prob2);
% 
% 
% figure(100)
% hold on
% plot(v_c1*1e6/10e3, right_prob2, '-o','DisplayName', 'base 2', 'LineWidth', 1)
% 
% legend
% ylabel('probability')
% xlabel('I_{bias1} (\muA)')
% title('Probability to count with zero bit errors in base 2')
% box on





%% counter perfect no snspd
load('Z:\SC\Measurements\SPG716\SNSPD_counter_1bit\SNSPD_counter1bit_osc\2bcounter_v1_798_v2_1160_v3_1800_vg_900\SPG716_2bcounter_v1_798_v2_1160_v3_1800_vg_900_SNSPD_counter_1bit_SNSPD_counter1bit_osc_2022-03-21 20-09-55.mat')
i = 900;
gain = 31.6;
figure
subplot(2,1,1)
plot(C1x(i,:)*1e6+0.16, C1y(i,:)*1e6/(9e3*2), 'DisplayName', ' Input 1st nT')
ylabel('current (\muA)')
xlim([0 2])
set(gca,'XTickLabel',[])
ylim([-10 40])
legend
subplot(2,1,2)
y = bandstop(C2y(i,:),[0.95e6 1.1e6],1/(C2x(i,2)-C2x(i,1)));
plot(C2x(i,:)*1e6+0.16, (y+0.030)*1e3/gain, 'DisplayName', ' 1st nT')
hold on
plot(C3x(i,:)*1e6+0.16, (C3y(i,:)+0.015)*1e3/gain, 'DisplayName', ' 2nd nT')
% y = C4y(i,:);
y = bandstop(C4y(i,:),[0.95e6 1.1e6],1/(C4x(i,2)-C4x(i,1)));
plot(C4x(i,:)*1e6+0.16, y*1e3/gain, 'DisplayName', ' 3rd nT')
xlabel('time (\mus)')
ylabel('voltage (mV)')
xlim([0 2])
ylim([-0.2 1.5])
legend

%%
z = linspace(1,10,1000);
x = exp(2*z-1);
y = bandstop(x,[100 200],1e3);
figure
plot(y)
%% snspd counter light on 
load('Z:\SC\Measurements\SPG716\counter_1bit\snspd_counter_1bit\snspd1bcounter_light4_3\SPG716_snspd1bcounter_light4_3_counter_1bit_snspd_counter_1bit_2022-03-22 02-45-18.mat')

for i = 380:10:1000
    figure(1)
    plot(C1x(i,:), C1y(i,:), 'DisplayName', ' SNSPD')
    hold on
    plot(C2x(i,:), C2y(i,:)-0.013, 'DisplayName', ' 1st nT')
    plot(C3x(i,:), C3y(i,:)-0.024, 'DisplayName', ' 2nd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(2)
end

%% 1bit counter and snspd no snspd output
load('Z:\SC\Measurements\SPG716\counter_1bit\snspd_counter_1bit\snspd1bcounter_light\SPG716_snspd1bcounter_light_counter_1bit_snspd_counter_1bit_2022-03-22 01-50-23.mat')


for i = 150:10:1000
    figure(1)
    plot(C2x(i,:), C2y(i,:)-0.013, 'DisplayName', ' 1st nT')
    hold on
    plot(C3x(i,:), C3y(i,:)-0.024, 'DisplayName', ' 2nd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(2)
end

%% osc 1bit count snspd
load('Z:\SC\Measurements\SPG716\counter_1bit\snspd_counter_1bit\osc_snspd1bcounter_v1_366_v2_945_v3_112\SPG716_osc_snspd1bcounter_v1_366_v2_945_v3_112_counter_1bit_snspd_counter_1bit_2022-03-21 23-20-16.mat')

for i = 150:10:1000
    figure(1)
    plot(C2x(i,:), C2y(i,:)-0.013, 'DisplayName', ' 1st nT')
    hold on
    plot(C3x(i,:), C3y(i,:)-0.024, 'DisplayName', ' 2nd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(2)
end

%% 2bit counter snspd light
load('Z:\SC\Measurements\SPG716\counter_2bit\snspd_counter_2bit\snspd2bcounter_v1_370_v2_1245_v3_1320_v4_1000\SPG716_snspd2bcounter_v1_370_v2_1245_v3_1320_v4_1000_counter_2bit_snspd_counter_2bit_2022-03-21 23-06-05.mat')
for i = 18:2:200
    figure(1)
    plot(C1x(i,:), C1y(i,:), 'DisplayName', ' SNSPD')
    hold on
    plot(C2x(i,:), C2y(i,:)-0.02, 'DisplayName', ' 1st nT')
    plot(C3x(i,:), C3y(i,:)-0.03, 'DisplayName', ' 2nd nT')
    plot(C4x(i,:), C4y(i,:)-0.04, 'DisplayName', ' 3rd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(2)
end

%% 2bit counter snspd osc
load('Z:\SC\Measurements\SPG716\SNSPD_counter\SNSPD_counter2bit_osc\SNSPD_2bcounter_osc\SPG716_SNSPD_2bcounter_osc_SNSPD_counter_SNSPD_counter2bit_osc_2022-03-21 18-33-24.mat')
for i = 18:2:200
    figure(1)
    plot(C1x(i,:), C1y(i,:), 'DisplayName', ' SNSPD')
    hold on
    plot(C2x(i,:), C2y(i,:)-0.02, 'DisplayName', ' 1st nT')
    plot(C3x(i,:), C3y(i,:)-0.03, 'DisplayName', ' 2nd nT')
    plot(C4x(i,:), C4y(i,:)-0.04, 'DisplayName', ' 3rd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(2)
end


%% 1bit counter snspd osc
load('Z:\SC\Measurements\SPG716\counter_1bit\snspd_counter_1bit\osc_snspd1bcounter_v1_366_v2_945_v3_112\SPG716_osc_snspd1bcounter_v1_366_v2_945_v3_112_counter_1bit_snspd_counter_1bit_2022-03-21 23-20-16.mat')
for i = 48:2:200
    figure(1)
    plot(C1x(i,:), C1y(i,:), 'DisplayName', ' SNSPD')
    hold on
    plot(C2x(i,:), C2y(i,:)-0.015, 'DisplayName', ' 1st nT')
    plot(C3x(i,:), C3y(i,:)-0.03, 'DisplayName', ' 2nd nT')
    hold off
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend
    pause(2)
end

%% sweep SNSPD 2bit counter perfect  -------- new biterror
clear all
myFolder = 'Z:\SC\Measurements\SPG716\demo\demo_right_value\demo\';
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);

right_prob2 = zeros(1, length(matFiles));
detect_prob2 = zeros(1, length(matFiles));
% right_prob3 = zeros(1, length(matFiles));
% right_prob4 = zeros(1, length(matFiles));
% right_prob5 = zeros(1, length(matFiles));
% right_prob6 = zeros(1, length(matFiles));
% right_prob7 = zeros(1, length(matFiles));
freq = zeros(1, length(matFiles));
v_c1 = zeros(1, length(matFiles));
v_c2 = zeros(1, length(matFiles));
v_g = zeros(1, length(matFiles));

min_peak_height = 0.024;
time_tol = 4e-8;

for k = 1:length(matFiles)
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  base = 2;
  [right_prob2(k), detect_prob2(k), light_power(k), v_b1(k), v_b2(k), v_b3(k), v_g(k), v_snspd(k)] = compute_prob2bit_snspd_biterr(fullFileName, min_peak_height, time_tol, base);
  fprintf(1, 'For P_light = %.3f, v_b1 = %.3f, v_b2 = %.3f, v_b3 = %.3f, v_g = %.3f, v_snspd = %.3f\n', light_power(k), v_b1(k), v_b2(k), v_b3(k), v_g(k), v_snspd(k));
  fprintf(1, 'the probablity to detect all the photons is %f\n', detect_prob2(k));
  fprintf(1, 'the probablity to work in base 2 is %f\n', right_prob2(k));
  
%   base = 3;
%   [right_prob3(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 3 is %f\n', right_prob3(k));
%   
%   base = 4;
%   [right_prob4(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 4 is %f\n', right_prob4(k));
%   
%   base = 5;
%   [right_prob5(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 5 is %f\n', right_prob5(k));
%   
%   base = 6;
%   [right_prob6(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 6 is %f\n', right_prob6(k));
%   
%   base = 7;
%   [right_prob7(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 7 is %f\n', right_prob7(k));

end

% figure
% hold on
% plot(v_c1*1e6/10e3, right_prob2, '-o','DisplayName', 'base 2', 'LineWidth', 1)
% plot(v_c1*1e6/10e3, detect_prob2, '-o','DisplayName', 'base 2', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob3, '-o', 'DisplayName', 'base 3','LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob4, '-o', 'DisplayName', 'base 4', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob5, '-o', 'DisplayName', 'base 5', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob6, '-o', 'DisplayName', 'base 6', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob7, '-o', 'DisplayName', 'base 7', 'LineWidth', 1)
% 
% legend
% ylabel('probability')
% xlabel('I_{bias1} (\muA)')
% title('Probability to count with zero bit errors in different bases (30 spikes burst)')
% box on

%% sweep SNSPD 2bit counter perfect  -------- new biterror - fine sweep
clear all
myFolder = 'Z:\SC\Measurements\SPG716\fine_sweep_snspd_cnt\snspd_cnt\demo\';
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);

right_prob2 = zeros(1, length(matFiles));
detect_prob2 = zeros(1, length(matFiles));
% right_prob3 = zeros(1, length(matFiles));
% right_prob4 = zeros(1, length(matFiles));
% right_prob5 = zeros(1, length(matFiles));
% right_prob6 = zeros(1, length(matFiles));
% right_prob7 = zeros(1, length(matFiles));
freq = zeros(1, length(matFiles));
v_c1 = zeros(1, length(matFiles));
v_c2 = zeros(1, length(matFiles));
v_g = zeros(1, length(matFiles));

min_peak_height = 1.5e-3;
time_tol = 1e-7;


for k = 1:length(matFiles)
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  base = 2;
  [right_prob2(k), detect_prob2(k), light_power(k), v_b1(k), v_b2(k), v_b3(k), v_g(k), v_snspd(k)] = compute_prob2bit_snspd_biterr(fullFileName, min_peak_height, time_tol, base);
  fprintf(1, 'For P_light = %.3f, v_b1 = %.3f, v_b2 = %.3f, v_b3 = %.3f, v_g = %.3f, v_snspd = %.3f\n', light_power(k), v_b1(k), v_b2(k), v_b3(k), v_g(k), v_snspd(k));
  fprintf(1, 'the probablity to detect all the photons is %f\n', detect_prob2(k));
  fprintf(1, 'the probablity to work in base 2 is %f\n', right_prob2(k));
  
%   base = 3;
%   [right_prob3(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 3 is %f\n', right_prob3(k));
%   
%   base = 4;
%   [right_prob4(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 4 is %f\n', right_prob4(k));
%   
%   base = 5;
%   [right_prob5(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 5 is %f\n', right_prob5(k));
%   
%   base = 6;
%   [right_prob6(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 6 is %f\n', right_prob6(k));
%   
%   base = 7;
%   [right_prob7(k), freq3, v_c13, v_c23, v_g3] = compute_prob2bit(fullFileName, min_peak_height, time_tol, base);
%   fprintf(1, 'the probablity to work in base 7 is %f\n', right_prob7(k));

end

% figure
% hold on
% plot(v_c1*1e6/10e3, right_prob2, '-o','DisplayName', 'base 2', 'LineWidth', 1)
% plot(v_c1*1e6/10e3, detect_prob2, '-o','DisplayName', 'base 2', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob3, '-o', 'DisplayName', 'base 3','LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob4, '-o', 'DisplayName', 'base 4', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob5, '-o', 'DisplayName', 'base 5', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob6, '-o', 'DisplayName', 'base 6', 'LineWidth', 1)
% % plot(v_c1*1e6/10e3, right_prob7, '-o', 'DisplayName', 'base 7', 'LineWidth', 1)
% 
% legend
% ylabel('probability')
% xlabel('I_{bias1} (\muA)')
% title('Probability to count with zero bit errors in different bases (30 spikes burst)')
% box on

%%
load('snspd_cnt.mat');

v1 = unique(v_b1);
v2 = unique(v_b2);    
for i = 1 : numel(v_b2)
    for k = 1 : numel(v2)
        for j = 1 : numel(v1)
            if v_b1(i)==v1(j) && v_b2(i)==v2(k)
                    prob_out2(k, j) = right_prob2(i);
                    det(k, j) = detect_prob2(i);
            end
        end
    end
end

figure
imagesc(v1*1e6/9e3, v2*1e6/9e3, 1-prob_out2)
tit = sprintf('BER');
title(tit)
xlabel('I_{bias1} (\muA)')
ylabel('I_{bias2} (\muA)')
set(gca,'YDir','normal')
colormap(flip(parula))
axis square
colorbar


figure
imagesc(v1*1e6/9e3, v2*1e6/9e3, det)
tit = sprintf('Detection Efficiency');
title(tit)
xlabel('I_{bias1} (\muA)')
ylabel('I_{bias2} (\muA)')
set(gca,'YDir','normal')
axis square
colorbar



%%



