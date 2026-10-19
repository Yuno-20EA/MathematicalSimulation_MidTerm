function floor_contours(x, y, F, levels, z0)
%FLOOR_CONTOURS  Ve duong dong muc cua F tai mat phang z = z0 (dung contourc).
    C = contourc(x, y, F, levels);
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        plot3(C(1, k+1:k+n), C(2, k+1:k+n), z0*ones(1, n), ...
              'Color', [0.3 0.3 0.3], 'LineWidth', 0.4);
        k = k + n + 1;
    end
end
