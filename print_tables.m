function print_tables(x0, Rsd, Rnt, opts)
%PRINT_TABLES  In ket qua theo bo cuc: moi phuong phap mot bang
%   Iteration | f(x) | Step Length, roi dong ket luan.
    fprintf('\n--- STARTING FROM POINT [%.1f, %.1f] ---\n\n', x0);
    print_one('Steepest Descent (Backtracking Line Search)', Rsd, opts);
    print_one('Newton''s Method (Backtracking Line Search)', Rnt, opts);
    if Rnt.iters < Rsd.iters
        fprintf('Newton converged in fewer iterations.\n');
    end
    fprintf('Summary: Steepest Descent %d iters | Newton %d iters | x_SD = (%.6f, %.6f), x_Newton = (%.6f, %.6f)\n', ...
        Rsd.iters, Rnt.iters, Rsd.x, Rnt.x);
    if Rnt.fallbacks > 0
        fprintf('Note: Newton direction replaced by -grad %d time(s) (Hessian not positive definite).\n', Rnt.fallbacks);
    end
end

function print_one(name, R, opts)
    fprintf('%s\n', name);
    fprintf('Iteration      f(x)   |   Step Length\n');
    n = R.iters;
    for k = 1:min(n, opts.printMax), print_row(k, R, opts); end
    if n > opts.printMax
        nOmit = n - opts.printMax - 1;
        if nOmit > 0
            fprintf('      ...   (%d iterations omitted)   ...\n', nOmit);
        end
        print_row(n, R, opts);
    end
    fprintf('\n');
end

function print_row(k, R, opts)
    fval = R.f(k+1);                                  % f(x_k) sau k vong
    if fval >= 1e-3, fs = sprintf('%.4f', fval); else, fs = sprintf('%.1e', fval); end
    fmt = sprintf('%%9d | %%9s | %%.%df\n', opts.stepDecimals);
    fprintf(fmt, k, fs, R.s(k));
end
