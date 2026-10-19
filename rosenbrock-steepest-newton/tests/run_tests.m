function run_tests()
%RUN_TESTS  Kiem thu tu dong (khong mo hinh). Chay:  >> run_tests
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'src'));
    opts = default_options();
    nFail = 0;

    nFail = nFail + check('gradient khop sai phan huu han', @() t_gradient());
    nFail = nFail + check('Hessian khop sai phan huu han',  @() t_hessian());
    nFail = nFail + check('Hessian tai x* xac dinh duong, kappa ~ 2.5e3', @() t_hess_star());
    nFail = nFail + check('backtracking thoa Armijo va s = rho^m', @() t_backtracking(opts));
    nFail = nFail + check('SD hoi tu ve (1,1) tu ca hai diem', @() t_converge(opts, 'steepest'));
    nFail = nFail + check('Newton hoi tu ve (1,1) tu ca hai diem', @() t_converge(opts, 'newton'));
    nFail = nFail + check('f giam don dieu, moi s <= s0', @() t_monotone(opts));
    nFail = nFail + check('Newton it vong lap hon SD it nhat 100 lan', @() t_newton_faster(opts));

    if nFail == 0, fprintf('\nALL TESTS PASSED\n'); else, error('%d test that bai', nFail); end
end

function n = check(name, fn)
    try, fn(); fprintf('PASS  %s\n', name); n = 0;
    catch err, fprintf('FAIL  %s -> %s\n', name, err.message); n = 1; end
end

function t_gradient()
    rng_state = [0.3 -0.7; 1.2 1.2; -1.2 1; 0.9 0.8];
    for r = 1:size(rng_state, 1)
        x = rng_state(r,:).';  [~, g] = rosenbrock(x);  gn = zeros(2,1);  h = 1e-6;
        for k = 1:2, e = zeros(2,1); e(k) = h;
            gn(k) = (rosenbrock(x+e) - rosenbrock(x-e)) / (2*h); end
        assert(norm(g - gn) / max(1, norm(g)) < 1e-6, 'sai lech gradient');
    end
end

function t_hessian()
    x = [-1.2; 1];  [~, ~, H] = rosenbrock(x);  Hn = zeros(2);  h = 1e-6;
    for k = 1:2, e = zeros(2,1); e(k) = h;
        [~, gp] = rosenbrock(x+e);  [~, gm] = rosenbrock(x-e);
        Hn(:,k) = (gp - gm) / (2*h); end
    assert(norm(H - Hn, 'fro') / norm(H, 'fro') < 1e-6, 'sai lech Hessian');
end

function t_hess_star()
    [~, ~, H] = rosenbrock([1; 1]);  ev = eig(H);
    assert(all(ev > 0), 'Hessian khong xac dinh duong');
    kappa = max(ev) / min(ev);
    assert(kappa > 2400 && kappa < 2600, 'kappa = %g', kappa);
end

function t_backtracking(opts)
    x = [-1.2; 1];  [f0, g0] = rosenbrock(x);  p = -g0;
    s = backtracking(x, p, f0, g0, opts);
    assert(rosenbrock(x + s*p) <= f0 + opts.c*s*(g0.'*p), 'Armijo bi vi pham');
    m = log(s) / log(opts.rho);  assert(abs(m - round(m)) < 1e-9, 's khong phai rho^m');
    if s < opts.s0   % buoc lon gap doi phai vi pham Armijo
        s2 = s / opts.rho;
        assert(rosenbrock(x + s2*p) > f0 + opts.c*s2*(g0.'*p), 'co buoc qua tay');
    end
end

function t_converge(opts, method)
    for i = 1:size(opts.starts, 2)
        R = line_search_descent(opts.starts(:,i), method, opts);
        assert(R.converged, 'khong hoi tu (iters=%d)', R.iters);
        assert(norm(R.x - [1;1]) < 1e-4, 'sai nghiem: %s', mat2str(R.x, 6));
    end
end

function t_monotone(opts)
    for i = 1:size(opts.starts, 2)
        for m = {'steepest', 'newton'}
            R = line_search_descent(opts.starts(:,i), m{1}, opts);
            assert(all(diff(R.f) <= 0), 'f tang');
            assert(all(R.s <= opts.s0 + eps) && all(R.s > 0), 'buoc ngoai (0, s0]');
        end
    end
end

function t_newton_faster(opts)
    for i = 1:size(opts.starts, 2)
        a = line_search_descent(opts.starts(:,i), 'steepest', opts);
        b = line_search_descent(opts.starts(:,i), 'newton',   opts);
        assert(a.iters >= 100 * b.iters, 'SD %d vs Newton %d', a.iters, b.iters);
    end
end
