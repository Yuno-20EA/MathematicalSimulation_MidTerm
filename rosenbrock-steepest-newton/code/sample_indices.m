function idx = sample_indices(m, n)
%SAMPLE_INDICES  Lay ~n chi so tu 1..m, doi xung theo thang log, giu nguyen 15 diem dau.
    idx = unique([1:min(m, 15), round(logspace(0, log10(m), n)), m]);
end
