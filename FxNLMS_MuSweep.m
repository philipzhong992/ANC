clear; clc; close all;

% 参数设置
mu_list = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.4];
fs = 1500;
nfft = 8192;
results = [];  % 存储结果：mu, 收敛点, 降噪量

% Step 1: 数据加载
data = readmatrix('Ref & Err.txt');
ref = data(:, 3);
d = data(:, 4);
N = length(ref);

% Step 1.5: 数据读取
hp_order = 64;  % 滤波器阶数（可根据实际情况调整）
cutoff_freq = 50;  % 截止频率 (Hz)
Wn = cutoff_freq / (fs/2);  % 归一化截止频率
b_hp = fir1(hp_order, Wn, 'high');  % 设计高通FIR滤波器
d = filtfilt(b_hp, 1, d);  % 零相位滤波器作用于 d

% 加载真实次级路径
load('2ndPath.mat', 'w');
if isrow(w), w = w'; end
h_secondary = w;

% 循环测试不同步长
for i = 1:length(mu_list)
    mu = mu_list(i);
    fprintf('\n=== 测试步长 mu = %.3f ===\n', mu);

    % 初始化
    clear SISO_FxNLMS
    y_array = zeros(N, 1);
    e_array = zeros(N, 1);

    for n = 1:N
        u = ref(n);
        [y, ~] = SISO_FxNLMS(u, 0, mu);
        y_array(n) = y;

        % 模拟经过次级路径
        if n >= length(h_secondary)
            y_filtered = sum(y_array(n - length(h_secondary) + 1:n) .* flipud(h_secondary));
        else
            y_padded = [zeros(length(h_secondary) - n, 1); y_array(1:n)];
            y_filtered = sum(y_padded .* flipud(h_secondary));
        end

        % 更新误差
        e_array(n) = d(n) + y_filtered;

        % 再次更新权重
        SISO_FxNLMS(u, e_array(n), mu);
    end

    %% 评估主峰降噪量
    window = hamming(length(d));
    D = abs(fft(d .* window, nfft));
    E = abs(fft(e_array .* window, nfft));
    D_dB = 20*log10(D(1:nfft/2) + 1e-8);
    E_dB = 20*log10(E(1:nfft/2) + 1e-8);
    f = linspace(0, fs/2, nfft/2);

    [peak_val, peak_idx] = max(D_dB);
    peak_freq = f(peak_idx);
    val_D = D_dB(peak_idx);
    val_E = E_dB(peak_idx);
    dn_dB = val_D - val_E;

    %% 评估收敛速度
    e_sq = e_array.^2;
    e_mov = movmean(e_sq, 200);
    final_val = mean(e_sq(end-200+1:end));
    threshold = final_val * 1.05;
    converge_idx = find(e_mov < threshold, 1);
    converge_time = converge_idx / fs;

    % 记录
    results = [results; mu, converge_idx, converge_time, dn_dB];

    % 可选：画出每条曲线
    figure(10);
    plot(f, D_dB, 'b--'); hold on;
    plot(f, E_dB, 'r'); title('频谱对比');
    legend('ANC前','ANC后'); xlim([0 300]); grid on;
    title(sprintf('mu = %.3f, 降噪量 %.2f dB', mu, dn_dB));
    ylim([-20 40])
    
    % 标记主峰位置
    plot(peak_freq, val_D, 'bo', 'MarkerSize', 6, 'LineWidth', 1.5);
    plot(peak_freq, val_E, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    % 添加注释
    text(peak_freq + 10, val_D, ...
         sprintf('%.1f Hz: %.1f dB (前)', peak_freq, val_D), ...
         'Color', 'b', 'FontSize', 10);
    text(peak_freq + 10, val_E, ...
         sprintf('%.1f Hz: %.1f dB (后)', peak_freq, val_E), ...
         'Color', 'r', 'FontSize', 10);

    pause(0.5);
    hold off
end

% 输出结果表格
fprintf('\n====== 汇总结果 ======\n');
fprintf('mu\t收敛点\t时间(s)\t主峰降噪(dB)\n');
for i = 1:size(results, 1)
    fprintf('%.3f\t%d\t%.2f\t%.2f\n', results(i,1), results(i,2), results(i,3), results(i,4));
end

% 可视化对比图
figure;
subplot(2,1,1);
plot(results(:,1), results(:,2), '-o');
xlabel('\mu'); ylabel('收敛样本数'); title('步长对收敛速度影响'); grid on;

subplot(2,1,2);
plot(results(:,1), results(:,4), '-o');
xlabel('\mu'); ylabel('主峰降噪量 (dB)'); title('步长对降噪性能影响'); grid on;
