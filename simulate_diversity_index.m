clear;
clc;
k=10;
for i=1:k
x(1,i) =15000000;
y(1,i) =100000;
% z(1) =200; 
z(1,i) =800; 
M(1,i)=10000;
F(1,i)=7000;
% M(1)=70000;
% F(1)=70000;
s(1,i)=500000;
end
    
    Kx =15000000;
    Ks =500000;
    
    r=0.15;
    r1=0.2;
    
    kpr=1;
    
%     mx=0.1;
%     mmf=0.2;
%     ms=0.2;
%     alpha =0.56+(0.78-0.56)/(1+exp(0.6*(y)/Ky));
    mx=0.0000005;
    mmf=0.0001;
    ms=0.00005;

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
for i=1:k
   for n=1:100
    Ky =0.1*(x(n,i));
    Kz =0.15*(M(n,i)+F(n,i));
    Kl =0.1*(s(n,i));
    qz =0.24*(1-z(n,i)/Kz);
    %alpha =0.56+(0.78-0.56)/(1+exp(0.6*(y(n))/Ky));    
    alpha =0.56+(0.78-0.56)*i/k;
    x_ =x(n,i)*r*(1-x(n,i)/Kx)-y(n,i)*(mx*x(n,i));
    y_ =y(n,i)*(-dy-dg)+30*F(n,i)*exp(-y(n,i)/Ky);
    %y_ =y.*(-dy-dg)+300*F*exp(1-y/Ky);
    z_ =z(n,i)*(qz-dz);
    M_ =(alpha).*dg*y(n,i)*exp(-(M(n,i)+F(n,i))/Kl)-(mmf*(M(n,i)+F(n,i)))*z(n,i)-dm*M(n,i);
    F_ =y(n,i)*(1-alpha)*dg*exp(-(M(n,i)+F(n,i))/Kl)-(mmf*(M(n,i)+F(n,i)))*z(n,i)-df*F(n,i);
    %s_ =s.*r1.*exp(1-s/Ks)-kpr.*ds.*(M+F).*(ms*s);  
    s_ =s(n,i)*r1*exp(-s(n,i)/Ks)-kpr*ds*(M(n,i)+F(n,i))*(ms*s(n,i));
    
    sum=x(n,i)+y(n,i)+M(n,i)+F(n,i)+z(n,i)+s(n,i);
    P1=x(n,i)/sum;
    P2=(y(n,i)+M(n,i)+F(n,i))/sum;
    P3=z(n,i)/sum;
    P4=s(n,i)/sum;
    sn(n,i)=-P1*log(P1)-P2*log(P2)-P3*log(P3)-P4*log(P4);

    x(n+1,i)=x(n,i)+x_;
    y(n+1,i)=y(n,i)+y_;
    z(n+1,i)=z(n,i)+z_;
    M(n+1,i)=M(n,i)+M_;
    F(n+1,i)=F(n,i)+F_;
    s(n+1,i)=s(n,i)+s_;
   end
end
for i=1:k-3
    plot(sn(:,i));
    hold on;
end
legend('1','2','3','4','5','6','7');