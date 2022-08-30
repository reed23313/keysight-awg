function [right, no_peaks]  = compare_peaks2bit(t1, nt1, t2, nt2, t3, nt3, base, min_peak_height_1, min_peak_height_2, time_tol)

    right = 0;
    
    if abs(max(nt1)) > 0.015
        min_peak_height_1 = max(nt1)/1.2;
    else
        min_peak_height_1 = 1000*max(nt1);
    end
    if abs(max(nt2)) > 0.025 
        min_peak_height_2 = max(nt2)/1.4;
    else
        min_peak_height_2 = 1000*max(nt2);
    end
    if abs(max(nt3)) > 0.025
        min_peak_height_3 = max(nt3)/1.4;
    else
        min_peak_height_3 = 1000*max(nt3);
    end
    
    

    [peak_1, i_peak_1] = findpeaks(nt1, 'MinPeakProminence', min_peak_height_1, 'MinPeakDistance', 3e-8/(t1(2)- t1(1)));
    t_peak_1 = t1(i_peak_1);

    [peak_2, i_peak_2] = findpeaks(nt2, 'MinPeakProminence', min_peak_height_2,'MinPeakDistance', 3e-8/(t2(2)- t2(1)));
    t_peak_2 = t2(i_peak_2);
    
    [peak_3, i_peak_3] = findpeaks(nt3, 'MinPeakProminence', min_peak_height_3,'MinPeakDistance', 3e-8/(t3(2)- t3(1)));
    t_peak_3 = t3(i_peak_3);
    
    
    seq = zeros(1, numel(t_peak_1));
    seq2 = zeros(1, numel(t_peak_1));
    
%     numel(t_peak_1)
    
    for i = 1 : numel(t_peak_1)
        same_t = (abs(t_peak_2 - t_peak_1(i)) < time_tol); 
        seq(i) = any(same_t);
        
        same_t_2 = (abs(t_peak_3 - t_peak_1(i)) < time_tol); 
        seq2(i) = any(same_t_2);
    end
    
%     figure(1) 
%     subplot(3,1,1)
%     plot(t1,nt1); hold on; plot(t_peak_1,peak_1); hold off;
%     subplot(3,1,2)
%     plot(t2,nt2); hold on; plot(t_peak_2,peak_2); hold off;
%     subplot(3,1,3)
%     plot(t3,nt3); hold on; plot(t_peak_3,peak_3); hold off;
%     pause(1)
    no_peaks = 0;
    err_one2 = 0;
    N = numel(t_peak_1);
    right_vect = zeros(1, N);
    right_vect2 = zeros(1, N);
    base2 = base^2;
    %|| sum(seq) >= 1.1*numel(seq)/base || sum(seq) <= 0.9*numel(seq)/base
    if N>=4
        if not(any(seq)) && not(any(seq2))
            right = 0; %finireeeee
        else
            right_seq = zeros(1, N);
            right_seq2 = zeros(1, N);
            one = find(seq == 1);
            one2 = find(seq2 == 1);
            if numel(one) == 0
                one = base+1;
            end
            if numel(one2) == 0
                one2 = base2+1;
            end

            if one(1)<=base 
                right_seq(one(1) : base : numel(right_seq)) = 1;
            else
                right_seq(base : base : numel(right_seq)) = 1;
            end

            yaxis = (one(1):base:one(1)+base2-base);
            [~, idx] = min(abs(one2(1)-yaxis.'));
            right_seq2(yaxis(idx) : base2 : numel(right_seq)) = 1;

            num_one = 0;
            for i = 1 : numel(right_seq)
                if num_one == base
                    num_one = 0;
                end
                if seq(i) == 1
                    num_one = num_one + 1;
                end

                right_vect(i) = not(xor(right_seq(i), seq(i)));

                if right_vect(i) == 0 && seq(i) == 1
                    right_seq(i+1:end) = zeros(1, N-i);
                    right_seq(i+base : base : numel(right_seq)) = 1;

                    if err_one2 == 0
                        right_seq2(i:end) = zeros(1, N-i+1);
                        right_seq2(i+base*(base-num_one) : base2 : numel(right_seq2)) = 1;  
                    end
                end

                if right_vect(i) == 0 && seq(i) == 0
                    right_seq(i+1:end) = zeros(1, N-i);
                    right_seq(i+1 : base : numel(right_seq)) = 1;

                    if err_one2 == 0
                        right_seq2(i:end) = zeros(1, N-i+1);
                        right_seq2(i+1+base*(base-num_one-1) : base2 : numel(right_seq)) = 1;
                    end

                end

                right_vect2(i) = not(xor(right_seq2(i), seq2(i)));

                err_one2 = 0;
                if right_vect2(i) == 0 && seq2(i) == 1
                    right_seq2(i+1:end) = zeros(1, N-i);
                    right_seq2(i+base2 : base2 : numel(right_seq2)) = 1;
                end

                if right_vect2(i) == 0 && seq2(i) == 0
                    right_seq2(i+1:end) = zeros(1, N-i);
                    right_seq2(i+1 : base2 : numel(right_seq2)) = 1;
                    err_one2 = 1;
                end

            end
        end
        if not(any(seq))
            right_vect = zeros(1, N);
        end
        if not(any(seq2)) 
            right_vect2 = zeros(1, N);
        end
        if N~=0
            right = mean([right_vect right_vect2]);
        end
    else
        no_peaks = 1;
    end
end

