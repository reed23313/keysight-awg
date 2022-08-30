function [detect, BER, V_in, V_b1, V_b2, max_N_pulses] = compute_prob1bit_biterr(mat, time_tol, base)
    load(mat)
    
%     s1 = find(peaks1 == -2);
%     s2 = find(peaks2 == -2);
%     s3 = find(peaks3 == -2);
%     s4 = find(peaks4 == -2);
%     N = numel(s1)-1;
    N = size(peaks1,1);
    max_N_pulses = 0;

    for i = 1 : N
        p1 = peaks1(i,:);
        p2 = peaks2(i,:);
        p3 = peaks3(i,:);

        p1(p1==0) = [];
        p2(p2==0) = [];
        p3(p3==0) = [];
        
        if numel(p1) > max_N_pulses
            max_N_pulses = numel(p1);
        end

        t_peak_1 = tvec_awg(p1);
        t_peak_2 = tvec_daq(p2);
        t_peak_3 = tvec_daq(p3);
        
        [detect(i), no_peaks_det] = compare_peaks_detection_biterr(t_peak_1, t_peak_2, time_tol);
        [right(i), no_peaks] = compare_peaks1bit_biterr(t_peak_1, t_peak_3, base, time_tol);

        V_in = vin_bias123(:,1);
        V_b1 = vin_bias123(:,2);
        V_b2 = vin_bias123(:,3);
        
        
    end 
    BER = 1-right;


end


