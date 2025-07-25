clc; clear all; close all;
%% calculate the transmission rate using FEM method in 1mm
tar = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\targetSurface_1mm.dat');
cont = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\contrastSurface_1mm.dat');

tr = tar(:,3)./cont(:,3);
trr = reshape(tr,91,46);

figure(Position=[300 300 380 290])
xdata = tar(1:91:4186,1);
ydata = tar(1:91,2)/1e3;
imagesc(xdata, ydata, trr);
set(gca,'YDir','Normal')
clim([0 1])
colorbar;
ylabel('\it{f}\rm{/(kHz)}')
xlabel('\it{L}\rm{_m/(mm)}')
set(gca,'FontName','Times New Roman','FontSize',18)

%% calculate the transmission rate using FEM method in 2mm
tar = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\targetSurface_2mm.dat');
cont = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\contrastSurface_2mm.dat');

tr = tar(:,3)./cont(:,3);
trr = reshape(tr,91,46);

figure(Position=[300 400 380 290])
xdata = tar(1:91:4186,1);
ydata = tar(1:91,2)/1e3;
imagesc(xdata, ydata, trr);
set(gca,'YDir','Normal')
clim([0 1])
colorbar;
ylabel('\it{f}\rm{/(kHz)}')
xlabel('\it{L}\rm{_m/(mm)}')
set(gca,'FontName','Times New Roman','FontSize',18)

%% calculate the transmission rate using FEM method in 3mm
tar = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\targetSurface_3mm.dat');
cont = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\contrastSurface_3mm.dat');

tr = tar(:,3)./cont(:,3);
trr = reshape(tr,91,46);

figure(Position=[300 500 380 290])
xdata = tar(1:91:4186,1);
ydata = tar(1:91,2)/1e3;
imagesc(xdata, ydata, trr);
set(gca,'YDir','Normal')
clim([0 1])
colorbar;
ylabel('\it{f}\rm{/(kHz)}')
xlabel('\it{L}\rm{_m/(mm)}')
set(gca,'FontName','Times New Roman','FontSize',18)

%% calculate the transmission rate using FEM method in 1mm-PLA
tar = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\targetSurface_ptfe.dat');
cont = load('D:\项目文件\8. 导波_论文\5. 超分辨率系列\2. 超材料透镜\3. 过程数据\透射率计算\contrastSurface_ptfe.dat');

tr = tar(:,3)./cont(:,3);
trr = reshape(tr,91,46);

figure(Position=[400 500 380 290])
xdata = tar(1:91:4186,1);
ydata = tar(1:91,2)/1e3;
imagesc(xdata, ydata, trr);
set(gca,'YDir','Normal')
clim([0 1])
colorbar;
ylabel('\it{f}\rm{/(kHz)}')
xlabel('\it{L}\rm{_m/(mm)}')
set(gca,'FontName','Times New Roman','FontSize',18)
