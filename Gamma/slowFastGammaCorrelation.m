% Calculate the correlation between slow and fast gamma
% Code from Laura Colgin.
function [slow_fast_corrcoeff,slow_fast_corrcoeff_p] = slowFastGammaCorrelation(non_overlap_peak_ind_unique_slow,non_overlap_peak_ind_unique_fast,kstart_save, kstop_save)

N = length(kstart_save);

W_ID_slow = zeros(1,N);
W_ID_fast = zeros(1,N);
for m = 1:N
    for n = 1:length(non_overlap_peak_ind_unique_slow)
        if non_overlap_peak_ind_unique_slow(n) > kstart_save(m) && non_overlap_peak_ind_unique_slow(n) < kstop_save(m)
            W_ID_slow(m) = 1;
        end
    end
      
    for n = 1:length(non_overlap_peak_ind_unique_fast)
        if non_overlap_peak_ind_unique_fast(n) > kstart_save(m) && non_overlap_peak_ind_unique_fast(n) < kstop_save(m)
            W_ID_fast(m) = 1;
        end
    end
end


index_slow = find(W_ID_slow == 1);

index_fast = find(W_ID_fast == 1);


index_union_fast_slow = union(index_fast,index_slow);

W_ID_fast_union_only = W_ID_fast(index_union_fast_slow);
W_ID_slow_union_only = W_ID_slow(index_union_fast_slow);


[slow_fast_corrcoeff,slow_fast_corrcoeff_p] = corrcoef(W_ID_fast_union_only,W_ID_slow_union_only);



