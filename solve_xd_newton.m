function xd = solve_xd_newton(k_sup, k, L, y, z, tol, max_iter)
    % 输入参数
    % k_sup: 波数 k_sup(f)，标量
    % k: 波数 k(f)，标量
    % L: 光栅长度向量 [L1, L2, ..., LM]
    % y: y 坐标向量 [y1, y2, ..., yM]
    % z: 相位环绕数向量 [z1, z2, ..., zM]
    % tol: 收敛阈值（例如 1e-6）
    % max_iter: 最大迭代次数（例如 100）

    % 初始值：选为光栅的几何中心
    xd = mean(L)+0.15;
    
    % 方程数
    M = length(L) - 1; % 共 M-1 个方程

    % 牛顿迭代
    for iter = 1:max_iter
        % 计算 g(xd) 和 g'(xd)
        g = zeros(M, 1);       % g(xd) 向量
        g_prime = 0;          % g'(xd) 累加值（雅可比矩阵按行累加）
        
        for m = 1:M
            % 当前光栅单元的参数
            Lm = L(m); Lm1 = L(m+1);
            ym = y(m); ym1 = y(m+1);
            zm = z(m); zm1 = z(m+1);
            
            % 计算 g_m(xd)
            term1 = k_sup * (Lm1 - Lm);
            term2 = k * (sqrt((xd - Lm1)^2 + ym1^2) - sqrt((xd - Lm)^2 + ym^2));
            g(m) = term1 + term2 - 2*pi * (zm - zm1);
            
            % 计算 g_m'(xd)
            grad1 = k * (xd - Lm1) / sqrt((xd - Lm1)^2 + ym1^2);
            grad2 = k * (xd - Lm) / sqrt((xd - Lm)^2 + ym^2);
            g_prime = g_prime + (grad1 - grad2); % 累加导数
        end
        
        % 更新 xd
        xd_new = xd - sum(g) / g_prime;
        
        % 检查收敛条件
        if abs(xd_new - xd) < tol
            fprintf('收敛至解 xd = %.6f after %d iterations.\n', xd_new, iter);
            xd = xd_new;
            return;
        end
        
        % 更新迭代值
        xd = xd_new;
    end
    
    % 未收敛警告
    warning('牛顿法未在最大迭代次数内收敛，最后的 xd = %.6f', xd);
end
