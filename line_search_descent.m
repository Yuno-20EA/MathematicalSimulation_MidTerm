function R = line_search_descent(x0, method, opts)
%LINE_SEARCH_DESCENT  Phuong phap giam (steepest descent / Newton) + backtracking.
%   R = line_search_descent(x0, 'steepest' | 'newton', opts)
%   Hai phuong phap dung CHUNG mot vong lap; chi khac nhau o huong tim kiem p_k,
%   nen moi khac biet ve ket qua deu quy duoc ve huong.
%
%   R.path     : 2 x (iters+1)  cac diem lap x_0..x_k
%   R.f, R.gnorm : f(x_k), ||grad f(x_k)|| (k = 0..iters)
%   R.s        : chieu dai buoc s_k (k = 1..iters)
%   R.fallbacks: so lan Newton phai thay bang -grad (Hessian khong xac dinh duong)
%   R.nEval    : tong so lan tinh f trong backtracking
    x = x0(:);
    [f, g] = rosenbrock(x);
    R.method = method;
    R.path = x;  R.f = f;  R.gnorm = norm(g);  R.s = [];
    R.fallbacks = 0;  R.nEval = 0;  k = 0;

    while R.gnorm(end) > opts.tol && k < opts.maxit
        switch lower(method)
            case 'steepest'
                p = -g;
            case 'newton'
                [~, ~, H] = rosenbrock(x);
                w = warning('off', 'all');
                p = -(H \ g);                         % giai H p = -g, khong tinh inv(H)
                warning(w);
                if ~all(isfinite(p)) || g.' * p >= 0  % khong phai huong giam
                    p = -g;  R.fallbacks = R.fallbacks + 1;
                end
            otherwise
                error('Phuong phap khong hop le: %s', method);
        end

        [s, ne] = backtracking(x, p, f, g, opts);     % Thuat toan 3.1
        x = x + s*p;
        [f, g] = rosenbrock(x);
        k = k + 1;

        R.path(:, end+1) = x;   R.f(end+1) = f;
        R.gnorm(end+1) = norm(g);  R.s(end+1) = s;  R.nEval = R.nEval + ne;
    end
    R.iters = k;  R.x = x;  R.converged = R.gnorm(end) <= opts.tol;
end
