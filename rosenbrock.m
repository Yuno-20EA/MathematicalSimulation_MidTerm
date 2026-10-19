function [f, g, H] = rosenbrock(x)
%ROSENBROCK  Ham Rosenbrock (3.27), gradient va Hessian.
%   f(x) = 100*(x2 - x1^2)^2 + (1 - x1)^2
%   [f, g, H] = rosenbrock(x) voi x la vector cot 2x1.
%   Chi tinh g, H khi duoc yeu cau (nargout) de cac lan goi chi can f ranh hon.
    x1 = x(1);  x2 = x(2);
    f = 100*(x2 - x1^2)^2 + (1 - x1)^2;
    if nargout > 1
        g = [-400*x1*(x2 - x1^2) - 2*(1 - x1); ...
              200*(x2 - x1^2)];
    end
    if nargout > 2
        H = [1200*x1^2 - 400*x2 + 2, -400*x1; ...
             -400*x1,                 200];
    end
end
