function [right, no_peaks]  = compare_peaks_detection_biterr(t_peak_1, t_peak_2, time_tol)

    right = 0;
    no_peaks = 0;
    

    seq = zeros(1, numel(t_peak_1));
    
%     numel(t_peak_1)

%     figure(1) 
%     subplot(2,1,1)
%     plot(t1,nt1); hold on; plot(t_peak_1,peak_1); hold off;
%     subplot(2,1,2)
%     plot(t2,nt2); hold on; plot(t_peak_2,peak_2); hold off;
%     pause(1)
%     
    for i = 1 : numel(t_peak_1)
        same_t = find((abs(t_peak_2 - t_peak_1(i)) < time_tol)); 
        if numel(same_t) > 0
            t_peak_2(same_t(1)) =  [];
            seq(i) = 1;
        end
    end

    
    if not(any(seq)) || numel(t_peak_1) == 0
        no_peaks = 1;
    else
        right = mean(seq);
    end
end

