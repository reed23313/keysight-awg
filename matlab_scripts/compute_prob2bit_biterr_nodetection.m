function [right_prob, freq, v_c1, v_c2, vg] = compute_prob2bit_biterr(mat, min_peak_height, time_tol, base)
    load(mat)
%     t1 = C2x;
    t1 = C1x;
    t2 = C3x;
    t3 = C4x;
%     nt1 = C2y;
    nt1 = C1y;
    nt2 = C3y;
    nt3 = C4y;
    right = zeros(1, size(t1, 1));
    right_noshift = zeros(1, size(t1, 1));
%     min_peak_height_1 = (v_c1/1.25)*min_peak_height;
%     min_peak_height_2 = (v_c2/1.25)*min_peak_height;
    min_peak_height_1 = min_peak_height*100;
    min_peak_height_2 = min_peak_height;
    
    
    for i = 1 : size(t1, 1)
%         C = max(nt1(i,:))/max(nt2(i,:));
        C = 1;
        right(i) = compare_peaks2bit_biterr(t1(i,:), nt1(i,:), t2(i,:), nt2(i,:)*C, t3(i,:), nt3(i,:)*C, base,  min_peak_height_1, min_peak_height_2, time_tol);

%         right_1(i) = r1(1);
%         right_2(i) = r2(1);
%         right_noshift_1(i) = r1(2);
%         right_noshift_2(i) = r2(2);
        
  
    end 
    
%     C1 = max(nt1(10,:))/max(nt2(10,:));
%     C2 = max(nt1(10,:))/max(nt3(10,:));
    
%     figure(1)
%     plot(t1(10,:), nt1(10,:))
%     hold on
%     plot(t2(10,:), nt2(10,:)*C1)
%     plot(t3(10,:), nt3(10,:)*C2)
%     hold off


    
    right_prob = mean(right);
end


