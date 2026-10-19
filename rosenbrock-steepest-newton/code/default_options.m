function opts = default_options()
%DEFAULT_OPTIONS  Tham so mac dinh cho toan bo du an.
    % --- diem xuat phat (moi cot la mot diem)
    opts.starts = [1.2, 1.2; -1.2, 1.0].';

    % --- Thuat toan 3.1 (backtracking line search)
    opts.s0   = 1;        % chieu dai buoc ban dau (theo de bai)
    opts.rho  = 0.5;      % he so co buoc,  0 < rho < 1
    opts.c    = 1e-4;     % hang so Armijo, 0 < c < 1
    opts.sMin = 1e-16;    % chan an toan: dung co buoc neu s qua nho

    % --- vong lap ngoai
    opts.tol   = 1e-6;    % dung khi ||grad f|| <= tol
    opts.maxit = 50000;   % so vong lap toi da

    % --- hien thi
    opts.printMax      = 20;    % so vong dau in ra moi bang (con lai: ..., vong cuoi)
    opts.stepDecimals  = 6;     % so chu so thap phan cua chieu dai buoc
    opts.zoomIters     = 15;    % so vong dau cua SD trong cua so thu phong (B)
    opts.zoomHalfWidth = 0.05;  % nua canh cua so thu phong quanh x* (C)

    % --- luu ket qua
    opts.saveFigures = true;
    opts.outDir      = 'results';
end
