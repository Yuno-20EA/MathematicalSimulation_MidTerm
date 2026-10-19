function [s, nEval] = backtracking(x, p, f0, g0, opts)
%BACKTRACKING  Thuat toan 3.1: tim kiem duong di lui (backtracking line search).
%   Bat dau tu s = s0, nhan lien tiep voi rho cho toi khi thoa dieu kien Armijo
%       f(x + s p) <= f(x) + c * s * g' p        (g' p < 0)
%   Dau vao : x (diem hien tai), p (huong giam), f0 = f(x), g0 = grad f(x)
%   Dau ra  : s (chieu dai buoc chap nhan), nEval (so lan tinh f da dung)
    slope = g0.' * p;
    if slope >= 0
        error('backtracking:notDescent', 'p khong phai huong giam (g''p = %g >= 0).', slope);
    end
    s = opts.s0;  nEval = 0;
    while true
        nEval = nEval + 1;
        if rosenbrock(x + s*p) <= f0 + opts.c * s * slope
            return;                                   % Armijo thoa
        end
        s = opts.rho * s;
        if s < opts.sMin
            warning('backtracking:tinyStep', 'Buoc nho hon sMin, dung backtracking.');
            return;
        end
    end
end
