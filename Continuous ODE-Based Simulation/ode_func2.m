%% 定义微分方程组的右端函数
function [dfdt,out2] = ode_func2(f)
   
    x = f(1); 
    y = f(2); 
    z = f(3); 
    M = f(4); 
    F = f(5); 
    s = f(6); 
    
    Kx =15000000;
    Ks =500000;
    Ky =0.1*(x);
    Kz =0.15*(M+F);
    Kl =0.1*(s);
    
    r=0.15;
    r1=0.2;
    
    qz =0.24*(1-z/Kz);
    kpr=1;
    
%     mx=0.1;
%     mmf=0.2;
%     ms=0.2;
%     alpha =0.56+(0.78-0.56)/(1+exp(0.6*(y)/Ky));
    mx=0.0000005;
    mmf=0.0001;
    ms=0.00005;
    alpha =0.7;
    
    dy=0.61;
    dg=0.25;
    ds=0.54;
    dz=0.1;
    dm=0.5;
    df=0.5;
    
%     dy=0.61;
%     dg=0.1;
%     ds=0.54;
%     dz=0.1;
%     dm=0.1;
%     df=0.1;
   
    x_ =x.*r*(1-x/Kx)-y.*(mx*x);
    y_ =y.*(-dy-dg)+30*F*exp(-y/Ky);
    %y_ =y.*(-dy-dg)+300*F*exp(1-y/Ky);
    z_ =z.*(qz-dz);
    M_ =(alpha).*dg*y*exp(-(M+F)/Kl)-(mmf*(M+F)).*z-dm*M;
    F_ =y.*(1-alpha).*dg*exp(-(M+F)/Kl)-(mmf*(M+F)).*z-df*F;
    %s_ =s.*r1.*exp(1-s/Ks)-kpr.*ds.*(M+F).*(ms*s);  
    s_ =s.*r1.*exp(-s/Ks)-kpr.*ds.*(M+F).*(ms*s);  
    
    dfdt = [x_;y_;z_;M_;F_;s_]; % 状态向量的导数
    out2  = [x,y,z,M,F,s];
end