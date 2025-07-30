function [y, W_out] = SISO_FxNLMS(u, e, mu_in)
    persistent W X S_hat L mu eps alpha

    if isempty(W)
        % 参数初始化
        L = 128;
        W = zeros(L, 1);
        X = zeros(L, 1);
        eps = 1e-6;
        alpha = 1e-4;
        mu = mu_in;  % 将传入的 mu 保存在内部
        load('2ndPath.mat', 'w');
        if isrow(w), w = w'; end
        if length(w) ~= L
            error('次级路径估计长度错误');
        end
        S_hat = w;
    end

    % 输入缓冲
    X = [u; X(1:end-1)];

    % Filtered-x
    x_f = conv(X, S_hat, 'same');

    % 输出
    y = W' * x_f;

    % 权重更新
    norm_factor = x_f' * x_f + eps;
    W = (1 - alpha) * W + mu * x_f * e / norm_factor;

    if nargout > 1
        W_out = W;
    end
end
