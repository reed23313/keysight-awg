function [right, no_peaks]  = compare_peaks2bit_biterr(t_peak_1, t_peak_2, t_peak_3, base, time_tol)

%     figure(1) 
%     subplot(3,1,1)
%     plot(t1,nt1); hold on; plot(t_peak_1,peak_1); hold off;
%     subplot(3,1,2)
%     plot(t2,nt2); hold on; plot(t_peak_2,peak_2); hold off;
%     subplot(3,1,3)
%     plot(t3,nt3); hold on; plot(t_peak_3,peak_3); hold off;
% %     pause(1)
    right = 0;
    no_peaks = 0;
    seq = zeros(1, numel(t_peak_1));
    seq2 = zeros(1, numel(t_peak_1));
    same_t_1 = [];
    same_t_2 = [];

    for i = 1 : numel(t_peak_1)
        if i == numel(t_peak_1)
            same_t_1 = (t_peak_2 > t_peak_1(i) - time_tol) & (t_peak_2 < t_peak_1(i) + time_tol); 
            seq(i) = any(same_t_1);
    
            same_t_2 = (t_peak_3 > t_peak_1(i) - time_tol) & (t_peak_3 < t_peak_1(i) + time_tol); 
            seq(i) = any(same_t_1);
        else
            same_t_1 = (t_peak_2 > t_peak_1(i) - time_tol) & (t_peak_2 < t_peak_1(i+1) - time_tol); 
            seq(i) = any(same_t_1);
    
            same_t_2 = (t_peak_3 > t_peak_1(i) - time_tol) & (t_peak_3 < t_peak_1(i+1) - time_tol); 
            seq(i) = any(same_t_1);
        end
       
    end
        
    
    err_one2 = 0;
    N = numel(t_peak_1);
    right_vect = zeros(1, N);
    right_vect2 = zeros(1, N);
    base2 = base^2;
        
    if N>4
        if not(any(seq)) && not(any(seq2))
            right = 0;
        else
            right_seq = zeros(1, N);
            right_seq(base : base : numel(right_seq)) = 1;
            right_seq2 = zeros(1, N);
            right_seq2(base2 : base2 : numel(right_seq)) = 1;
            
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
        right = mean([right_vect right_vect2]);  
    else
        no_peaks = 1;
    end

end

