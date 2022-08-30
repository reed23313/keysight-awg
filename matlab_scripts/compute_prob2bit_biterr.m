function [detect, BER, V_in, V_b1, V_b2, V_b3, max_N_pulses] = compute_prob2bit_biterr(mat, base, detection)
    load(mat)
    
%     s1 = find(peaks1 == -2);
%     s2 = find(peaks2 == -2);
%     s3 = find(peaks3 == -2);
%     s4 = find(peaks4 == -2);
%     N = numel(s1)-1;
    N = size(peaks1,1);
    max_N_pulses = 0;
    time_tol = 1/(2*double(freq));

    for i = 1 : N
        p1 = peaks1(i,:);
        p2 = peaks2(i,:);
        p3 = peaks3(i,:);
        p4 = peaks4(i,:);

        p1(p1==0) = [];
        p2(p2==0) = [];
        p3(p3==0) = [];
        p4(p4==0) = [];
        
        if numel(p1) > max_N_pulses
            max_N_pulses = numel(p1);
        end

        t_peak_1 = tvec_awg(p1);
        t_peak_2 = tvec_daq(p2);
        t_peak_3 = tvec_daq(p3);
        t_peak_4 = tvec_daq(p4);
        if detection
            [detect(i), no_peaks_det] = compare_peaks_detection_biterr(t_peak_1, t_peak_2, time_tol);
        else
            detect(i) = 1;
            no_peaks_det = 0;
        end
        [right(i), no_peaks] = compare_peaks2bit_biterr(t_peak_1, t_peak_3, t_peak_4, base, time_tol);

        V_in = vin_bias123(:,1);
        V_b1 = vin_bias123(:,2);
        V_b2 = vin_bias123(:,3);
        V_b3 = vin_bias123(:,4);

        if rem(i,max(1,(N/100))) == 0
            fprintf('%i\n', i*100/N)
        end
        
    end 
    BER = 1-right;


end


