function fig = plot_surface3d(starts, R, ~)
%PLOT_SURFACE3D  Hinh 1: be mat 3D Rosenbrock + 4 quy dao + chu thich (callout).
    zmax = 200;
    fig = figure('Name', 'Rosenbrock 3D', 'Position', [60 60 1000 650], 'Color', 'w');
    x1 = linspace(-2, 2, 220);  x2 = linspace(-0.5, 3, 220);
    [X1, X2] = meshgrid(x1, x2);
    F = 100*(X2 - X1.^2).^2 + (1 - X1).^2;

    hsurf = surf(X1, X2, min(F, zmax), 'EdgeColor', 'none');  hold on;
    try, set(hsurf, 'FaceAlpha', 0.7); catch, end   % be mat hoi trong suot de thay quy dao
    try, colormap(parula); catch, colormap(jet); end
    floor_contours(x1, x2, F, logspace(-1, 3.3, 45), 0);

    style = {'bs-', 'ro-'; 'cd-', 'g^-'};            % (diem xuat phat, phuong phap)
    names = {'Steepest Descent', 'Newton''s Method'};
    h = zeros(1, 4);  n = 0;  lbl = cell(1, 4);
    for i = 1:size(starts, 2)
        for j = 1:2
            n = n + 1;
            P = R{i,j}.path;
            idx = sample_indices(size(P, 2), 120);   % SD co ~13000 diem -> lay mau
            z = min(R{i,j}.f(idx), zmax - 5) + 4;    % nang nhe de khong bi be mat che
            h(n) = plot3(P(1,idx), P(2,idx), z, style{i,j}, 'LineWidth', 1.3, ...
                         'MarkerSize', 4, 'MarkerFaceColor', style{i,j}(1));
            lbl{n} = sprintf('%s from [%.1f, %.1f]', names{j}, starts(:,i));
        end
    end

    % --- chu thich co mui ten
    for i = 1:size(starts, 2)
        x0 = starts(:,i);  f0 = min(rosenbrock(x0), zmax - 5) + 4;
        callout([1.8 - 3.5*(i == 2), 2.6, 200], [x0; f0], ...
                sprintf('Start point [%.1f, %.1f]', x0));
    end
    [i, b] = first_backtracked(R);                  % buoc Newton dau tien bi co
    if ~isempty(i)
        q = R{i,2}.path(:, b);
        callout([q(1) + 0.9, q(2) + 1.3, 200], [q; min(R{i,2}.f(b), zmax) + 4], ...
                'Backtracking occurring here');
    end
    plot3(1, 1, 4, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
    callout([1.9, 1.0, 40], [1 1 4], 'Ideal - [1, 1]');

    xlabel('x_1'); ylabel('x_2'); zlabel('f(x)');
    zlim([0 zmax]);  view(-30, 50);  grid on;
    title('Rosenbrock Function (Banana Function) Optimization: Steepest Descent vs. Newton', 'FontSize', 10);
    legend(h, lbl, 'Location', 'northwest', 'FontSize', 8);
end

function callout(textPos, target, str)
    plot3([textPos(1) target(1)], [textPos(2) target(2)], [textPos(3) target(3)], ...
          'k-', 'LineWidth', 1);
    text(textPos(1), textPos(2), textPos(3), str, 'FontWeight', 'bold', ...
         'BackgroundColor', 'w', 'EdgeColor', 'k', 'HorizontalAlignment', 'center');
end

function [i, b] = first_backtracked(R)
    i = [];  b = [];
    for r = 1:size(R, 1)
        k = find(R{r,2}.s < 1, 1);
        if ~isempty(k), i = r; b = k; return; end
    end
end
