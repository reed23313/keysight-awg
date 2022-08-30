clear all
%fullFileName = 'C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_10Mhzinput_2022-08-16_223250';

%%good one base2
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep_2022-08-16_183744.mat'; 
%good one base3
%fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base3_2022-08-16_225131.mat';
% %base3 less fine
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base3_2022-08-16_230236.mat'; 
%base2 less fine
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base3_2022-08-16_230509.mat'; 
%best base2
fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base2fine_2022-08-17_001230.mat'; 
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base2fine_2022-08-17_004604.mat'; 
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base2fine_2022-08-17_012857.mat'; 
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base2fine_2022-08-17_014657.mat'; 
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base2fine_2022-08-17_015435.mat'; %%%%10MHz
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter2bit_sweep035_base2fine_2022-08-17_020108.mat'; %%40MHz
% fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter1bit_2022-08-17_171949.mat'; %%%%%2bit counter 50 MHz
fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter1bit_2022-08-17_174737.mat'; %%%%good 10MHz delay


detection = 0;
% time_tol = 0.8e-7;
base = 2;
[detect, BER_base2, v_in, v_b1, v_b2, v_b3, max_N_pulses] = compute_prob2bit_biterr(fullFileName, base, detection);

VIN = v_in(1);
VB3 = v_b3(1);

v1 = unique(v_b1);
v2 = unique(v_b2);    
for i = 1 : numel(v_b2)
    if v_in(i) == VIN && v_b3(i) == VB3
        for k = 1 : numel(v2)
            for j = 1 : numel(v1)
                if v_b1(i)==v1(j) && v_b2(i)==v2(k)
                    B(k, j) = BER_base2(i);
                    det(k, j) = detect(i);
                end
            end
        end
    end
end

figure
R = 5.11e3;
imagesc(2*v1*1e6/R, 2*v2*1e6/R, B)
tit = sprintf('BER');
title(tit)
xlabel('I_{bias1} (\muA)')
ylabel('I_{bias2} (\muA)')
set(gca,'YDir','normal')
colormap(flip(parula))
set(gca,'ColorScale','log')
caxis([1/max_N_pulses 1])
axis square
colorbar

if detection
    figure
    imagesc(2*v1*1e6/R, 2*v2*1e6/R, det)
    tit = sprintf('Detection Efficiency');
    title(tit)
    xlabel('I_{bias1} (\muA)')
    ylabel('I_{bias2} (\muA)')
    set(gca,'YDir','normal')
    set(gca,'ColorScale','log')
    axis square
    colorbar
end

%% 1bit

clear all

fullFileName ='C:\Users\Demo\Documents\keysight-awg\counter1bit_2022-08-17_021349.mat'; %%50MHz

time_tol = 1e-8;
base = 2;
[detect, BER_base2, v_in, v_b1, v_b2, max_N_pulses] = compute_prob1bit_biterr(fullFileName, time_tol, base);
base = 3;
[detect, BER_base3, v_in, v_b1, v_b2, max_N_pulses] = compute_prob1bit_biterr(fullFileName, time_tol, base);
base = 4;
[detect, BER_base4, v_in, v_b1, v_b2, max_N_pulses] = compute_prob1bit_biterr(fullFileName, time_tol, base);
base = 5;
[detect, BER_base5, v_in, v_b1, v_b2, max_N_pulses] = compute_prob1bit_biterr(fullFileName, time_tol, base);
base = 6;
[detect, BER_base6, v_in, v_b1, v_b2, max_N_pulses] = compute_prob1bit_biterr(fullFileName, time_tol, base);
base = 7;
[detect, BER_base7, v_in, v_b1, v_b2, max_N_pulses] = compute_prob1bit_biterr(fullFileName, time_tol, base);

figure
R = 5.11e3;
plot(2*v_b1*1e6/R, BER_base2)
hold on
plot(2*v_b1*1e6/R, BER_base3)
plot(2*v_b1*1e6/R, BER_base4)
plot(2*v_b1*1e6/R, BER_base5)
plot(2*v_b1*1e6/R, BER_base6)
plot(2*v_b1*1e6/R, BER_base7)
legend

% figure
% R = 5.11e3;
% imagesc(2*v_b1*1e6/R, 2*v2*1e6/R, B)
% tit = sprintf('BER');
% title(tit)
% xlabel('I_{bias1} (\muA)')
% ylabel('I_{bias2} (\muA)')
% set(gca,'YDir','normal')
% colormap(flip(parula))
% set(gca,'ColorScale','log')
% caxis([1/max_N_pulses 1])
% axis square
% colorbar
% 
% figure
% imagesc(2*v1*1e6/R, 2*v2*1e6/R, det)
% tit = sprintf('Detection Efficiency');
% title(tit)
% xlabel('I_{bias1} (\muA)')
% ylabel('I_{bias2} (\muA)')
% set(gca,'YDir','normal')
% set(gca,'ColorScale','log')
% axis square
% colorbar