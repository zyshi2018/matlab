# matlab离散系统李雅普诺夫指数
```
clear all;clc
%% 主要是计算
b1=1;
b2=1;
d=0.5;
c1=1.08;
c2=1;
z=0.5;
m=0.25;
p=0.5;
T=0.1;
S=0.1;
k2=0.15;

N1 = 100;
N2 = 400;
k10 = 0:0.001:0.4;
f = zeros(N1+N2,length(k10));
f2 = zeros(N1+N2,length(k10));


for kk=1:length(k10)
    x0=2.5;
    y0=2.35;
    k1=k10(kk)
    L1=0;L2=0;
    %     s=zeros(2,1);
    J1=eye(2);
    for j = 1:N1+N2
        x1=x0+ k1.*x0.*(a-2*b1.*x0+d.*y0+b1.*(c1+(1-m).*z-S));
        y1=y0+ k2.*y0.*(a-2*b2.*y0+d.*x0+b2.*(c2-p.*m.*z+T));
        f(j,kk) =x1;
        f2(j,kk) =y1;
        j1=1 + a*k1 + d*k1*y1 + b1*k1 *(c1 - S - 4*x1 + z - m*z);
        j2=d*k1*x1;
        j3=d*k2*y1;
        j4=1 + a*k2 + d *k2 *x1 + b2 *k2 *(c2 + T - 4 *y1 - m *p* z);
        J=[j1 j2;j3 j4];
        x0=x1;y0=y1;
    end
end

hold on
f = f(N1+1:end,:);
plot(k10,f,'b.','MarkerSize',1)
ylabel('Price');

hold on
f2 = f2(N1+1:end,:);
plot(k10,f2,'r.','MarkerSize',1)
xlabel('\mu');
ylabel('Price');

% plot(k10,Lm2,'r','linewidth',0.5)
plot(k10,Lm1,'b')
line([0 0.4],[0 0])
```
