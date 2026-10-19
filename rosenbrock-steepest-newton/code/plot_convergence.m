function fig = plot_convergence(starts, R, ~)
%PLOT_CONVERGENCE  Hinh 3: f(x_k), ||grad f(x_k)|| va chieu dai buoc s_k theo vong lap
%   (truc k theo thang log vi Steepest Descent can ~10^4 vong, Newton chi ~10).
    nS = size(starts, 2);
    fig = figure('Name', 'Convergence', 'Position', [40 40 1500 400*nS], 'Color', 'w');
    floorv = 1e-30;
    for i = 1:nS
        Rs = R{i,1};  Rn = R{i,2};
        subplot(nS, 3, (i-1)*3 + 1);
        loglog(1:Rs.iters, max(Rs.f(2:end), floorv), 'b-', 'LineWidth', 1.2); hold on;
        loglog(1:Rn.iters, max(Rn.f(2:end), floorv), 'ro-', 'LineWidth', 1.2, 'MarkerFaceColor', 'r');
        grid on;  xlabel('iteration k'); ylabel('f(x_k)');
        title(sprintf('x_0 = [%.1f, %.1f]: f(x_k)  (f^* = 0)', starts(:,i)));
        legend('Steepest Descent', 'Newton''s Method', 'Location', 'southwest');

        subplot(nS, 3, (i-1)*3 + 2);
        loglog(1:Rs.iters, max(Rs.gnorm(2:end), floorv), 'b-', 'LineWidth', 1.2); hold on;
        loglog(1:Rn.iters, max(Rn.gnorm(2:end), floorv), 'ro-', 'LineWidth', 1.2, 'MarkerFaceColor', 'r');
        grid on;  xlabel('iteration k'); ylabel('||\nabla f(x_k)||');
        title('gradient norm');

        subplot(nS, 3, (i-1)*3 + 3);
        loglog(1:Rs.iters, Rs.s, 'b.', 'MarkerSize', 4); hold on;
        loglog(1:Rn.iters, Rn.s, 'ro-', 'LineWidth', 1.2, 'MarkerFaceColor', 'r');
        ylim([3e-4 2]);  grid on;  xlabel('iteration k'); ylabel('step length s_k');
        title('step length s_k (log scale; s_0 = 1)');
    end
end
