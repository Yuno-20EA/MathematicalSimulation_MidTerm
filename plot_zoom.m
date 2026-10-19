function fig = plot_zoom(starts, R, opts)
%PLOT_ZOOM  Hinh 2: thu phong tung quy dao, moi hang la mot diem xuat phat.
%   (A) toan canh, co khung (B), (C)
%   (B) thu phong vao zoomIters vong dau cua Steepest Descent
%   (C) thu phong quanh nghiem x* = (1,1)
    nS = size(starts, 2);
    fig = figure('Name', 'Zoom comparison', 'Position', [40 40 1500 440*nS], 'Color', 'w');
    for i = 1:nS
        Rs = R{i,1};  Rn = R{i,2};
        winB = padded_box(Rs.path(:, 1:min(opts.zoomIters + 1, end)), 0.25);
        hw = opts.zoomHalfWidth;
        winC = [1-hw 1+hw; 1-hw 1+hw];
        winA = padded_box([Rs.path, Rn.path, [1; 1]], 0.12);

        subplot(nS, 3, (i-1)*3 + 1);
        draw_panel(winA, Rs, Rn, 0, 0, ...
                   sprintf('x_0 = [%.1f, %.1f] - (A) overview', starts(:,i)), {winB, 'B'; winC, 'C'}, i == 1);
        subplot(nS, 3, (i-1)*3 + 2);
        draw_panel(winB, Rs, Rn, opts.zoomIters + 1, 12, ...
                   sprintf('(B) zoom: first %d iterations of Steepest Descent', opts.zoomIters), {}, false);
        subplot(nS, 3, (i-1)*3 + 3);
        draw_panel(winC, Rs, Rn, 0, 0, '(C) zoom: neighbourhood of x^* = [1, 1]', {}, false);
    end
end

function win = padded_box(P, frac)
    lo = min(P, [], 2);  hi = max(P, [], 2);
    pad = frac * max(max(hi - lo), 1e-3);
    win = [lo - pad, hi + pad];
end

function draw_panel(win, Rs, Rn, nMarkSD, nLabelSD, ttl, boxes, withLegend)
    xl = win(1,:);  yl = win(2,:);
    [X1, X2] = meshgrid(linspace(xl(1), xl(2), 160), linspace(yl(1), yl(2), 160));
    F = 100*(X2 - X1.^2).^2 + (1 - X1).^2;
    lv = logspace(log10(max(min(F(:)), 1e-10)), log10(max(F(:))), 28);
    [~, hc] = contour(X1, X2, F, lv);  hold on;
    set(hc, 'LineColor', [0.75 0.75 0.75]);

    hs = plot(Rs.path(1,:), Rs.path(2,:), 'b-', 'LineWidth', 1);
    if nMarkSD > 0
        m = min(size(Rs.path, 2), nMarkSD);
        plot(Rs.path(1,1:m), Rs.path(2,1:m), 'bs', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    end
    hn = plot(Rn.path(1,:), Rn.path(2,:), 'ro-', 'LineWidth', 1.4, ...
              'MarkerSize', 5, 'MarkerFaceColor', 'r');
    plot(Rs.path(1,1), Rs.path(2,1), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot(1, 1, 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');

    dx = 0.015 * diff(xl);  dy = 0.015 * diff(yl);
    for k = 1:min(size(Rn.path, 2) - 1, 15)          % danh so vong lap Newton
        label_point(Rn.path(:,k+1), k, 'r', xl, yl, dx, dy);
    end
    for k = 1:min(size(Rs.path, 2) - 1, nLabelSD)    % danh so vong lap SD
        label_point(Rs.path(:,k+1), k, 'b', xl, yl, dx, -dy);
    end
    for b = 1:size(boxes, 1)                         % khung thu phong
        w = boxes{b,1};
        rectangle('Position', [w(1,1) w(2,1) diff(w(1,:)) diff(w(2,:))], ...
                  'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.2);
        text(w(1,2), w(2,2), [' ' boxes{b,2}], 'FontWeight', 'bold');
    end
    set(gca, 'XLim', xl, 'YLim', yl);
    xlabel('x_1'); ylabel('x_2'); title(ttl, 'FontSize', 9);
    if withLegend
        legend([hs hn], {'Steepest Descent', 'Newton''s Method'}, 'Location', 'southeast');
    end
end

function label_point(q, k, col, xl, yl, dx, dy)
    if q(1) >= xl(1) && q(1) <= xl(2) && q(2) >= yl(1) && q(2) <= yl(2)
        text(q(1) + dx, q(2) + dy, num2str(k), 'Color', col, 'FontSize', 8, 'FontWeight', 'bold');
    end
end
