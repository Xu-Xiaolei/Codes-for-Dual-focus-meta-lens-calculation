clc;clear all;close all
%% 假设以下参数已知（根据实际情况修改）
xf = 0.25;              % 焦点的x坐标
yf = 0;                 % 焦点的y坐标
f = 120e3;              % 频率（示例值，需根据实际情况确定）
ki = 250;               % ki(f)的值（示例值，需根据实际情况确定）
k = 210;                % k的值（示例值，需根据波长等实际情况确定）
n = 2;                  % 不同的n值（可根据需要调整范围）
Ain = 1;
C = 1;
d = 2e-3;               % 衍射缝隙宽度
t = 10e-6;
%% 计算干涉系数alpha
alpha = pi / (2 * xf);
yy = -0.057:0.006:0.057;
M = length(yy);
for m = 1:M 
    % 计算不同n值下的Lm长度
    pf = [xf, yf]; % 假设焦点在x轴上，y坐标为0
    H_pf = k*sqrt((ki^2 - k^2) * (yf - yy(m))^2 + ki^2 * xf^2 + 4 * n * pi * ki * xf + 4 * n^2 * pi^2);
    Lm(m) = (k^2 * xf + 2 * ki * n * pi - H_pf) / (k^2 - ki^2);
end
%% 计算总的衍射场
m = 0;n = 0;j = sqrt(-1);
xres = 0.05:0.001:0.65;
yres = -0.06:0.001:0.06;
A = zeros(length(xres),length(yres));
for x = 0.05:0.001:0.65
    m = m+1;
    n = 0;
    for y = -0.06:0.001:0.06
        n = n+1;        
        for i = 1:M   
            dst = sqrt((x-Lm(i))^2+(y-yy(i))^2);
            A(m,n) = A(m,n)+exp(j*(-ki*Lm(i)-k*dst))/sqrt(dst);
        end
        A(m,n) = A(m,n)*Ain*C*exp(j*2*pi*f*t)*d*(1/M);
        
        if x < 0.10
            A(m,n) = 0;
        % elseif x < xf  
        %     A(m,n) = A(m,n)*0.9996^(((x-xf)^2)*15000);
        end
        % if y > yf+0.001 || y < yf-0.001
        %     A(m,n) = A(m,n)*0.99^(((y-yf)^2)*15000);
        % end
    end
end
%% 使用imagesc()函数绘制图像
figure(Position=[500 500 450 400])
subplot(2,1,1)
imagesc(xres,yres, abs(A.').^4*1e5);
% clim([0,3]*1e-4);
% colorbar; % 显示颜色条
colormap hot;
set(gca,'xtick',[]);
hold on
for i = 1:length(Lm)
    plot([0,Lm(i)],[yy(i),yy(i)],LineWidth=3,Color=[239 199 39]*0.9/256);
end

r = xf/2;
t = 0:0.01:2*pi; % 参数t的范围，步长为0.01
x = r*cos(t); 
y = r*sin(t);
plot(x,y,LineStyle="-",LineWidth=2.5,Color=[0.0 0.2 0.9]);
hold on
r = 2*r; % 圆的半径

t = 0:0.01:2*pi; % 参数t的范围，步长为0.01
x = r*cos(t); 
y = r*sin(t);
plot(x,y,LineStyle=":",LineWidth=2.5,Color=[0.1 0.8 0.1]);
hold on
r = 2*r; % 圆的半径

t = 0:0.01:2*pi; % 参数t的范围，步长为0.01
x = r*cos(t); 
y = r*sin(t);
plot(x,y,LineStyle="--",LineWidth=2.5,Color=[0.0 0.7 0.8]);
hold on


ylabel('Y(m)')
set(gca,'FontName','Times New Roman','FontSize',16)
%% 计算干涉系数alpha
clear all;
xf = 0.25;              % 焦点的x坐标
yf = 0;                 % 焦点的y坐标
f = 120e3;              % 频率（示例值，需根据实际情况确定）
ki = 250;               % ki(f)的值（示例值，需根据实际情况确定）
k = 210;                % k的值（示例值，需根据波长等实际情况确定）
n = 2;                  % 不同的n值（可根据需要调整范围）
Ain = 1;
C = 1;
d = 2e-3;               % 衍射缝隙宽度
t = 10e-6;
alpha = pi / (2 * xf);
yy = -0.057:0.006:0.057;
M = length(yy);
% 计算不同n值下的Lm长度
for m = 1:M 
    % 计算不同n值下的Lm长度
    pf = [xf, yf]; % 假设焦点在x轴上，y坐标为0
    H_pf = k*sqrt((ki^2 - k^2) * (yf - yy(m))^2 + ki^2 * xf^2 + 4 * n * pi * ki * xf + 4 * n^2 * pi^2);
    Lm(m) = (k^2 * xf + 2 * ki * n * pi - H_pf) / (k^2 - ki^2);
end
%% 计算总的衍射场
m = 0; n = 2; j = sqrt(-1);
rres = 0.05:0.005:0.65;  % r坐标范围
thetares = -pi/2:0.005:pi/2;  % θ坐标范围
A = zeros(length(rres), length(thetares));

for r = 0.05:0.005:0.65
    m = m + 1;
    n = 0;
    for theta = -pi/2:0.005:pi/2
        n = n + 1;
        % 计算x和y的值
        x = r * cos(theta);
        y = r * sin(theta);
        % 计算从每个Lm(i)到(r,theta)的距离dst
        for i = 1:M   
            dst = sqrt((x - Lm(i))^2 + (y - yy(i))^2);
            A(m,n) = A(m,n) + exp(j*(-ki*Lm(i)-k*dst))/sqrt(dst);
        end
        A(m,n) = A(m,n) * Ain * C * exp(j*2*pi*f*t) * d * (1/M);
        if r < 0.10
            A(m,n) = 0;
        elseif r < xf-0.05 || r> xf+0.05
             A(m,n) = A(m,n)*0.999^(((r-xf)^2)*15000);
        end
    end
end
%% 使用imagesc()函数绘制图像
AA=abs(A.').^1.5*1e5*0.8;
subplot(2,1,2)
imagesc(rres, thetares, AA);  % 使用极坐标数据
colorbar;  % 显示颜色条
% clim([0,8]*1e-2);
colormap hot;  % 使用jet颜色映射
hold on
plot([xf xf],[-2 2],LineStyle=":",LineWidth=2.5,Color=[0.1 0.8 0.1]);
hold on
plot([xf/2 xf/2],[-2 2],LineStyle="-",LineWidth=2.5,Color=[0.0 0.2 0.9]);
hold on
plot([xf*2 xf*2],[-2 2],LineStyle="--",LineWidth=2.5,Color=[0.0 0.7 0.8]);
ylabel('\it{θ}\rm{/(rad)}')
xlabel('\it{r}\rm{/(m)}')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18)
%% 截出
figure(Position=[300 300 450 220])
k = find(rres==xf);
plot(thetares,AA(:,k),LineStyle=":",LineWidth=2,Color=[0.1 0.8 0.1]*0.8);
hold on
k = find(rres==0.5*xf);
plot(thetares,AA(:,k),LineStyle="-",LineWidth=1,Color=[0.0 0.2 0.9]);
hold on
k = find(rres==2*xf);
plot(thetares,AA(:,k),LineStyle="--",LineWidth=1,Color=[0.0 0.7 0.8]);

xlim([thetares(1) thetares(end)])
set(gca,'XGrid','on','FontName','Times New Roman','FontSize',18)

ylabel('\it{E}\rm{/(mJ)}')
xlabel('\it{θ}\rm{/(rad)}')
