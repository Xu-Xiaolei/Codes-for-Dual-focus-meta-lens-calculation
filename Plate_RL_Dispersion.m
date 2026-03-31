%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       %
%       平板Rayleigh Lamb的频散曲线      %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
global E;                        %定义杨氏模量为全局变量
global NU;                       %定义泊松比为全局变量
global rho;                      %定义密度为全局变量
global d;                        %定义板的厚度为全局变量
%材料参数初始化
%E=70300000000;                   %杨氏模量单位（Pa）(视需要修改)
% NU=0.33;                         %泊松比(视需要修改)
% rho=2700;                        %密度单位（Kg/m^5）(视需要修改)
% d=0.003;                         %板的厚度单位(m)(视需要修改)

E=70.3e9;                       %杨氏模量单位（Pa）(视需要修改)
NU=0.33;                         %泊松比(视需要修改)
rho=2700;                        %密度单位（Kg/m^5）(视需要修改)
d=0.003;                         %板的厚度单位(m)(视需要修改)
%求解模式的频率上限与相速度上限
f0=500;
fstep=500;
fend=1000000;                    %频率单位(Hz)
Cpstep=10;                       %速度单位（m/s)
Cpend=25000;                     %欲求解更多模式可以扩大速度和频率上限
i=0;
for f=f0:fstep:fend              %求解每一频率下的频散方程解
    fprintf('\n当前计算进度：%4.2f%%',f*100/fend);
    i=i+1;
    fd=f*d;
    s=0;
    Cp0=50;
    while Cp0<Cpend                     %二分法扫掠
        Cpdata=0;
        Cp1=Cp0;
        Cp2=Cp1+Cpstep;         
        y1=Symmetric(Cp1,fd);                %调用对称模式下的频散方程
        y2=Symmetric(Cp2,fd);                %调用对称模式下的频散方程
        if y1==0
            Cpdata=Cp1;
        elseif y2==0
            Cpdata=Cp2;
        elseif y1*y2<0
            while(Cp2-Cp1)>0.001     %二分法终值
                Cptem=(Cp1+Cp2)/2;
                ytem=Symmetric(Cptem,fd);
                y2=Symmetric(Cp2,fd);
                if ytem*y2>0
                    Cp2=Cptem;
                else
                    Cp1=Cptem;
                end
            end
            Cpdata=Cptem;
        end
        if Cpdata~=0
           s=s+1;
           Cp1data(s,i)=Cpdata;           %将结果存储为二维数组
        end
    Cp0=Cp0+Cpstep;
    end
    ph(i)=fd;                          %坐标横轴赋值
end
[m,n]=size(Cp1data);                      %数组维度
for u=1:m
    for v=1:n
        if Cp1data(u,v)==0
           Cp1data(u,v)=inf;
        end
    end
end
count=0;
for u=1:m                             %数据处理(将一个fd对应两个Cp的模态,前面数据会视为两个模态,所以要进行处理)
     for v=1:n-1
        s=0;
        if Cp1data(u,v+1)==inf &&Cp1data(u,v)~=inf
            count=count+1;
            xaxis(count)=u;
            yaxis(count)=v;
            s=0;
            for t=v:-1:1
                if Cp1data(u,t)~=inf
                    s=s+1;
                end
            end
            for k=1:s
                abcp1(count,k)=Cp1data(u,v-k+1);
                abph1(count,k)=ph(v-k+1);
                Cp1data(u,v-k+1)=inf;
            end
            for i=1:n-1
                if Cp1data(u-1,i)~=inf
                    for j=1:n-i+1;
                        abcp2(count,j)=Cp1data(u-1,i+j-1);
                        abph2(count,j)=ph(i+j-1);
                        %data(u-1,i+j-1)=inf;       %不无穷大化，用于计算群速度
                    end
                    break
                end
            end
        end
     end
end
if count~=0
    Cp2data=[abcp1 abcp2];%归并模态
    Cp2ph=[abph1 abph2];%归并模态对应fd
end
%计算群速度
for u=1:m
    for v=1:n-1
        if Cp1data(u,v)~=inf
            Cg(u,v)=Cp1data(u,v)^2/(Cp1data(u,v)-ph(v)*(Cp1data(u,v+1)-Cp1data(u,v))/(ph(v+1)-ph(v)));
        else
            Cg(u,v)=inf;
        end
    end
end
%计算波数
for u=1:m
    for v=1:n
        kk(u,v) = 2*pi*ph(v)/(Cp1data(u,v)*(1000*d));
    end
end
%反对称模式求解，各设置与上相同
i=0;
for f=f0:fstep:fend
    i=i+1;
    fd=f*d;
    s=0;
    Cp0=0.05;
    while Cp0<Cpend
        Cpdata=0;
        Cp1=Cp0;
        Cp2=Cp1+Cpstep;
        y1=Antisymmetric(Cp1,fd);   %调用反对称模式下的频散方程
        y2=Antisymmetric(Cp2,fd);   %调用反对称模式下的频散方程
        if y1==0
            Cpdata=Cp1;
        elseif y2==0
            Cpdata=Cp2;
        elseif y1*y2<0
            while(Cp2-Cp1)>0.000001
                Cptem=(Cp1+Cp2)/2;
                ytem=Antisymmetric(Cptem,fd);
                y2=Antisymmetric(Cp2,fd);
                if ytem*y2>0
                    Cp2=Cptem;
                else
                    Cp1=Cptem;
                end
            end
            Cpdata=Cptem;
        end
        if Cpdata~=0
            s=s+1;
            data1(s,i)=Cpdata;
        end
        Cp0=Cp0+Cpstep;
        
    end
end
[m,n]=size(data1);
for u=1:m
    for v=1:n
        if data1(u,v)==0
            data1(u,v)=inf;
        end
    end
end
for u=1:m
    for v=1:n-1
        if data1(u,v)~=inf&&data1(u,v+1)~=inf
            Cg1(u,v)=data1(u,v)^2/(data1(u,v)-ph(v)*(data1(u,v+1)-data1(u,v))/(ph(v+1)-ph(v)));
    
        else
            Cg1(u,v)=inf;
        end
    end
end
for u=1:m
    for v=1:n
        kk1(u,v) = 2*pi*ph(v)/(data1(u,v)*(1000*d));
    end
end
%% 画出结果
figure('Position', [500, 500, 500, 300]);
plot(ph(1:length(Cp1data))/1000,Cp1data/1000,'black','linewidth',1);
hold on;
plot(Cp2ph/1000,Cp2data/1000,'black','linewidth',1);
hold on;
plot(ph(1:length(data1))/1000,data1/1000,'r--','linewidth',1);
xlabel('f·d(MHz·mm)'),
ylabel('Cp(km/s)');
% grid on;
% title('Phase Velocity Dispersion Curves')
set(gca,"fontsize",15,'fontname','Times New Roman')
ylim([0 8])
xlim([0 3])
saveas(gcf, 'cp.svg');
%axis([0 20 0 20]);

figure('Position', [600, 500, 500, 300]);
plot(ph(1:length(Cg))/1000,Cg/1000,'black','linewidth',1);
hold on;
plot(ph(1:length(Cg))/1000,Cg1/1000,'r--','linewidth',1)
xlabel('f·d(MHz·mm)'),
ylabel('Cg(km/s)');
% grid on;
% title('Group Velocity Dispersion Curves')
set(gca,"fontsize",15,'fontname','Times New Roman')
ylim([0 8])
xlim([0 3])
saveas(gcf, 'cg.svg');
%axis([0 20 0 7]);

figure('Position', [600, 500, 500, 300]);
plot(ph(1:length(kk))/(1000*d),kk,'black','linewidth',1);
hold on;
plot(ph(1:length(kk1))/(1000*d),kk1,'r--','linewidth',1)
xlabel('f (kHz)'),
ylabel('k (rad/mm)');
xlim([0 300])
set(gca,"fontsize",15,'fontname','Times New Roman')
%axis([0 20 0 7]);