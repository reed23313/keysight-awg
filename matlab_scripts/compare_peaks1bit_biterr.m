function [right, no_peaks]  = compare_peaks1bit_biterr(t_peak_1, t_peak_2, base, time_tol)

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

    for i = 1 : numel(t_peak_1)
        same_t_1 = [];
        same_t_1 = find((abs(t_peak_2 - t_peak_1(i)) < time_tol)); 
        
        if numel(same_t_1) > 0
%             t_peak_2(same_t_1(1)) =  [];
            seq(i) = 1;
        end
    end
        
    
    N = numel(t_peak_1);
    right_vect = zeros(1, N);
        
    if N>4
        if not(any(seq))
            right = 0;
        else
            right_seq = zeros(1, N);
            right_seq(base : base : numel(right_seq)) = 1;
            
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
                   
                end
                
                if right_vect(i) == 0 && seq(i) == 0
                    right_seq(i+1:end) = zeros(1, N-i);
                    right_seq(i+1 : base : numel(right_seq)) = 1;
                    
                end          
              
                        
            end
        end
        if not(any(seq))
            right_vect = zeros(1, N);
        end
        right = mean(right_vect);  
    else
        no_peaks = 1;
    end

end

