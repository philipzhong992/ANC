clear all; clc; close all;

% Step 1: 数据读取
data = readmatrix('Ref & Err.txt');
ref = data(:, 3);      % 参考信号（u）
d = data(:, 4);        % 原始误差信号（ANC未启用）
fs = 1500;

% Step 2: 加载主路径（用于生成 d，可用于验证）、次级路径
load('1stPathFilter.mat', 'w');   % 主路径
if isrow(w), w = w'; end
h_primary = w;

load('2ndPathFilter.mat', 'w');   % 次级路径估计（用于ANC建模）
if isrow(w), w = w'; end
h_secondary = w;

% Step 3: 初始化
N = length(ref);
y_array = zeros(N, 1);
e_array = zeros(N, 1);

% 清空滤波器状态
clear SISO_FxNLMS

% Step 4: ANC仿真
for n = 1:N
    u = ref(n);   % 当前参考输入

    % 运行ANC算法，得到扬声器输出 y(n)
    y = SISO_FxNLMS(u, 0);  % 初始误差信号设置为0（仅用于初始化）

    % 记录 y，注意后续才真正用于计算误差
    y_array(n) = y;

    % 计算经过次级路径的滤波输出（反馈噪声）
    if n >= length(h_secondary)
        y_segment = y_array(n - length(h_secondary) + 1 : n);
        y_filtered = sum(y_segment .* flipud(h_secondary));  % 卷积结果 y*S_hat
    else
        y_padded = [zeros(length(h_secondary) - n, 1); y_array(1:n)];
        y_filtered = sum(y_padded .* flipud(h_secondary));
    end

    % 当前误差信号为：
    e_array(n) = d(n) + y_filtered;  % 注意是加，因为扬声器发出的是y

    % 用这个误差信号更新自适应滤波器
    SISO_FxNLMS(u, e_array(n));
end

% Step 5: 可视化对比
figure;
subplot(2,1,1);
plot(d, 'b'); hold on; plot(e_array, 'r');
legend('原始误差 d(n)', 'ANC开启后 e(n)');
title('误差信号时域波形对比');
xlabel('样本点'); ylabel('幅值');

subplot(2,1,2);

% 设置参数
nfft = 8192;
f = linspace(0, fs/2, nfft/2);  % 实际频率坐标

% 计算 FFT
window = hamming(length(d));   % 加窗防泄漏
D = abs(fft(d .* window, nfft));
E = abs(fft(e_array .* window, nfft));

% 转 dB，并加小值防止 log(0)
D_dB = 20*log10(D(1:nfft/2) + 1e-8);
E_dB = 20*log10(E(1:nfft/2) + 1e-8);

% 画图
plot(f, D_dB, 'b', 'LineWidth', 1); hold on;
plot(f, E_dB, 'r', 'LineWidth', 1);
legend('ANC前','ANC后');
xlabel('频率 (Hz)');
ylabel('幅度 (dB)');
title('误差信号频谱对比');
ylim([-20 40]);
xlim([0 300]);
grid on;

% ----- 自动识别 ANC 前的主峰频率 -----
[peak_val, peak_idx] = max(D_dB);     % ANC前的频谱峰值索引
peak_freq = f(peak_idx);              % 峰值对应的频率 (Hz)

% 获取 ANC 后在同一频率的 dB 值
val_D = D_dB(peak_idx);
val_E = E_dB(peak_idx);

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

% Step 6: 收敛速度评估
e_sq = e_array.^2;  % 误差信号能量
window_size = 200;  % 移动平均窗口长度
e_moving_avg = movmean(e_sq, window_size);  % 平滑曲线

% 找到初始平均值和稳定平均值
initial_energy = mean(e_sq(1:window_size));
final_energy = mean(e_sq(end - window_size + 1:end));

% 收敛阈值设置为最终能量 + 5%
threshold = final_energy * 1.05;

% 查找从头开始，第一个低于阈值的点
converge_idx = find(e_moving_avg < threshold, 1);

if isempty(converge_idx)
    fprintf('未在信号长度内收敛\n');
else
    fprintf('收敛至稳定区间所需样本点数：%d（约 %.2f 秒）\n', ...
        converge_idx, converge_idx/fs);
end

% 可视化误差能量趋势
figure;
plot(e_moving_avg, 'k'); hold on;
yline(threshold, 'r--', '阈值');
xline(converge_idx, 'b--', sprintf('收敛点 = %d', converge_idx));
xlabel('样本点'); ylabel('误差能量 (平方均值)');
title('误差信号能量收敛过程');
grid on;
