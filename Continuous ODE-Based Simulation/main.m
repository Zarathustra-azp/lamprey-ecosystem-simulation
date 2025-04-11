clc;clear;
close all;

T=100;dt=1;n=0;
%% 定义求解区间和初始值向量
x(1) =15000000;
y(1) =100000;
% z(1) =200; 
z(1) =800; 
M(1)=10000;
F(1)=7000;
% M(1)=70000;
% F(1)=70000;
S(1)=500000;
f = [x(1), y(1), z(1), M(1), F(1), S(1)];

for qwq=1:T/dt
    [dfdt, out(qwq,:)]= ode_func(f(qwq,:));
    f(qwq+1,:) = f(qwq,:)+dt.*(dfdt');
    n=n+1;
end


T=100;dt=1;n=0;
%% 定义求解区间和初始值向量
x =15000000;
y =100000;
z =800; 
M=10000;
F=7000;
S=500000;
f = [x(1), y(1), z(1), M(1), F(1), S(1)];

for qwq=1:T/dt
    [dfdt, out1(qwq,:)]= ode_func1(f(qwq,:));
    f(qwq+1,:) = f(qwq,:)+dt.*(dfdt');
    n=n+1;
end


T=100;dt=1;n=0;
%% 定义求解区间和初始值向量
x =15000000;
y =100000;
z =800; 
M=10000;
F=7000;
S=500000;
f = [x(1), y(1), z(1), M(1), F(1), S(1)];

for qwq=1:T/dt
    [dfdt, out2(qwq,:)]= ode_func2(f(qwq,:));
    f(qwq+1,:) = f(qwq,:)+dt.*(dfdt');
    n=n+1;
end