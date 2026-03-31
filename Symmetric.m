%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        %
%       平板的频散曲线的对称模态公式      %
%                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=Symmetric(Cp,fd)       %Cp相速度,fd平板的频率*厚度  (对称模态公式)
global E;                         %定义杨氏模量为全局变量
global NU;                        %定义泊松比为全局变量
global rho;                       %定义密度为全局变量
mu=E/(2*(1+NU));                  %剪切模量
lambda=NU*E/((1+NU)*(1-2*NU));    %Lame常数
Cl=sqrt((lambda+2*mu)/rho);       %纵波速度
Ct=sqrt(mu/rho);                  %横波速度
a=sqrt(1/Ct^2-1/Cp^2);            %对称与反对称模态公式中的q/w
b=sqrt(1/Cl^2-1/Cp^2);            %对称与反对称模态公式中的p/w
if Cp<Ct
    a=abs(a);
    b=abs(b);
    y=(a^2+1/Cp^2)^2*cosh(b*pi*fd)*sinh(a*pi*fd)-4*a*b*sinh(b*fd*pi)*cosh(a*fd*pi)/Cp^2;
elseif Ct<=Cp&&Cp<Cl
    b=abs(b);
    y=(a^2-1/Cp^2)^2*cosh(b*pi*fd)*sin(a*pi*fd)-4*a*b*sinh(b*pi*fd)*cos(a*fd*pi)/Cp^2;
else
    y=(a^2-1/Cp^2)^2*sin(a*pi*fd)*cos(b*pi*fd)+4*a*b*sin(b*pi*fd)*cos(a*pi*fd)/Cp^2;
end
if real(y)==0
    y=imag(y);
else
    y=real(y);
end
    
