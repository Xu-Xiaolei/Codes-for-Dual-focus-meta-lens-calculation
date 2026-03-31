clc;clear all; close all
%% 画能带图
data1 = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\1_l.csv');
data2 = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\2_l.csv');
data3 = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\3_l.csv');

n = 40;  m = 9; %n为布里渊扫描对数，m为计算阶数
x1 = data1(1:n, 1);
x2 = data2(1:n, 1);
x3 = data3(1:n, 1);

% 提取每个绘图的数据
y1 = zeros(n, m); 
y2 = zeros(n, m); 
y3 = zeros(n, m); 
for i = 1:m
    y1(:, i) = data1((i-1)*n+1:i*n, 2);
end
for i = 1:m
    y2(:, i) = data2((i-1)*n+1:i*n, 2);
end
for i = 1:m
    y3(:, i) = data3((i-1)*n+1:i*n, 2);
end

% 创建图形和坐标轴
fig = figure('Position', [500 500 370 570]);   % 设置图像尺寸
axes1 = axes('parent', fig, 'fontname', 'Times New Roman');

% 使用循环绘制所有 y 数据
hold on
colorAll = [134 237 142]*0.58/256;
for i = 1:m
    plot(x1, y1(:, i)/1e3, '-', 'LineWidth', 1.5,Color=colorAll);
    hold on
end
for i = 1:m
    plot(x2, y2(:, i)/1e3, '-', 'LineWidth', 1.5,Color=colorAll);
    hold on
end
for i = 1:m
    plot(x3, y3(:, i)/1e3, '-', 'LineWidth', 1.5,Color=colorAll);
    hold on
end
% 设置Y轴刻度和范围
% set(gca, 'YTick', 0:100:1000);  % 设置纵坐标刻度间隔为15
% axis([0 3 0 1000]);            % 设置坐标轴范围

% 设置X轴刻度和标签
fontsize = 16;
set(gca, 'XTick', [0 1 2 3]);
xtixlabel = {'M', 'Г', 'X', 'M'};
set(gca, 'XTickLabel', xtixlabel, 'FontName', 'Times New Roman', 'FontSize', fontsize);

% 添加坐标轴标签
xlabel('Wavenumber/(kπ/a)', 'FontSize', fontsize);   % X轴标签
ylabel('Frequency/(kHz)', 'FontSize', fontsize);      % Y轴标签

y=[y1;y2;y3];
x=[x1;x2;x3];

% 使用 fill 函数突出显示带隙
max1 = max(y(:, 6));   % 第二阶最小值
min1 = min(y(:, 7));    % 第一阶最大值
xx = [0, 3, 3, 0];
yy = [min1, min1, max1, max1]/1e3;
fill(xx, yy, 'k', 'EdgeColor', (1-colorAll)*0.2, 'FaceColor', (1-colorAll)*0.7,'LineWidth', 0.5, 'FaceAlpha', 0.3);  % 填充带隙区域
text(1.5,max1/1e3+4,strcat(num2str(round(max1/1e3,1)),' ~ ',num2str(round(min1/1e3,1)),'kHz'), ...
    'FontName','Times New Roman','FontSize',fontsize,'FontWeight','bold','Color',(1-colorAll)*0.6);

% 设置图形边框和图层
box on
set(gca, 'Layer', 'top', 'LineWidth', 1,'xgrid','on');

% 以期望的分辨率保存图像
print(gcf, '-dpng', '-r600', './img.png');

hold off

