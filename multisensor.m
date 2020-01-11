close
clear
clc

%Multisensor switching control strategy with fault tolerance guarantees
%% system parameters
tspan = 10;
t = 0:0.01:10;
A = [0.75 0.1
     0 0.75];
B = [0 0.1]';
E = B;

uref = 1./(1+exp(-0.5*t));
% uref = sin(2*pi*t./10);
u = uref;
eu = zeros(size(uref));
w = randn(1,length(t));
w = w / std(w);
w = w - mean(w);
w = 0 + 0.01 * w; % mean 0.1 variance 0.0001  
w = 1/2 * w; % |w|<0.04
% w = zeros(1,length(t));

figure(1)
plot(t,w,'linewidth',1)
grid on
xlabel('t / s')
ylabel('amplitude')
legend('white noise')

x1 = zeros(2,length(t));
x1up = x1;
y1 = x1;
ex1 = x1;
x2 = x1;
x2up = x2;
y2 = x2;
xref = x1;
ex2 = x2;
%% sensor parameters
s1 = zeros(size(x1));
s1up = s1;
As1 = 0.6065;
Bs1 = 0.5;
Cs1 = 0.7869;
n1 = 0.005*sin(10*pi*t);
% n1 = zeros(1,length(t));
s1ref = s1;
es1 = s1;
sm1 = s1;

s2 = zeros(size(s1));
s2up = s2;
% As2 = 0.8187;
% Bs2 = 0.5;
% Cs2 = 0.3625;
As2 = 0.6065;
Bs2 = 0.5;
Cs2 = 0.7869;
n2 = 0.005*cos(10*pi*t);
% n2 = zeros(1,length(t));
s2ref = s2;
sm2 = s2;

C = [1 0];
L1 = [1.1064 0.3010]';
% L2 = [1.1945 0.2985]';
L2 = [1.1064 0.3010]';
Ls1 = 1.3089;
% Ls2 = 2.8405;
Ls2 = 1.3089;
es2 = s2;

M = [A zeros(2,1);
     Bs1*C As1]\[L1; Ls1];
M1 = M(1:2);
Ms1 = M(3);
M = [A zeros(2,1);
     Bs2*C As2]\[L2; Ls2];
M2 = M(1:2);
Ms2 = M(3);

Q = [0.1007 0
    0 6.3187];
R = 7.2598;
K = [0.1118 0.0091];
P = [0.0236 0.8602
     0.8602 3.0214];

Al1 = [A zeros(2,1); Bs1*C As1] - [L1; Ls1]*[zeros(1,2) Cs1];
Al2 = [A zeros(2,1); Bs2*C As2] - [L2; Ls2]*[zeros(1,2) Cs2];
error1 = zeros(3,length(t));
error2 = error1;
esm1 = sm1;
esm2 = sm2;
z1 = x1;
z2 = x2;
ess1 = s1;
ess2 = s2;
z1up = x1;
z2up = x2;
r1 = es1;
r2 = es2;
x = x1;
zw1 = x1;
zw2 = x2;
flag = zeros(1,length(t));
K = [0.0118 0.0091];
%% Iteration
for i = 1:length(t) - 1
    if i>=200 && i<= 400
        s1(:,i) = 0;
    end
    if i>=600 && i<=800
        s2(:,i) = 0;
    end
    %control law
    eu(i) = u(i) - uref(i);
    esm1(:,i) = sm1(:,i) - s1ref(:,i);
    esm2(:,i) = sm2(:,i) - s2ref(:,i);
    z1(:,i) = x1(:,i) - xref(:,i);
    z2(:,i) = x2(:,i) - xref(:,i);
    ess1(:,i) = s1(:,i) - s1ref(:,i);
    ess2(:,i) = s2(:,i) - s2ref(:,i);
    z1up(:,i) = x1up(:,i) - xref(:,i);
    z2up(:,i) = x2up(:,i) - xref(:,i);
    r1(:,i) = M1 .* Cs1 .* es1(:,i) + M1 * n1(i);
    r2(:,i) = M2 .* Cs2 .* es2(:,i) + M2 * n2(i);
    
    zw1(i) = z1(:,i)' * P * z1(:,i);
    zw2(i) = z2(:,i)' * P * z2(:,i);
    if zw1(i) < zw2(i)
        u(i) = uref(i) - K *  (x1(:,i) - xref(:,i) + r1(:,i));
        flag(i) = 1;
    else
        u(i) = uref(i) - K * (x2(:,i) - xref(:,i) + r2(:,i));
        flag(i) = 2;
    end
%     u(i) = uref(i);
    %reference system
    xref(:,i+1) = A * xref(:,i) + B * uref(i);
    x(:,i+1) = A * x(:,i) + B * u(i) + E * w(i);
  
    %reference sensor 1
    s1ref(:,i+1) = As1 * s1ref(:,i) + Bs1 * C * xref(:,i);
    
    %reference sensor 2
    s2ref(:,i+1) = As2 * s2ref(:,i) + Bs2 * C * xref(:,i);
    
    %sensor 1 measurement
    sm1(:,i+1) = As1 * sm1(:,i) + Bs1 * C * x(:,i);
    y1(:,i) = Cs1 * sm1(:,i) + n1(i);
    
    %sensor 2 meansurement
    sm2(:,i+1) = As2 * sm2(:,i) + Bs2 * C * x(:,i);
    y2(:,i) = Cs2 * sm2(:,i) + n2(i);
    
    %sensor 1 state observer
    x1(:,i+1) = A * x1(:,i) + B * u(i) + L1 .* (y1(:,i) - Cs1 * s1(:,i));
    s1(:,i+1) = As1 * s1(:,i) + Bs1 * C * x1(:,i) + Ls1 .* (y1(:,i) - Cs1 * s1(:,i));
    x1up(:,i+1) = x1(:,i) + M1 .* (y1(:,i) - Cs1 * s1(:,i));
    s1up(:,i+1) = s1(:,i) + Ms1 .* (y1(:,i) - Cs1 * s1(:,i));
    
    %sensor 2 state observer
    x2(:,i+1) = A * x2(:,i) + B * u(:,i) + L2 .* (y2(:,i) - Cs1 * s2(:,i));
    s2(:,i+1) = As2 * s2(:,i) + Bs2 * C * x2(:,i) + Ls2 .* (y2(:,i) - Cs2 * s2(:,i));
    x2up(:,i+1) = x2(:,i) + M2 .* (y2(:,i) - Cs2 * s2(:,i));
    s2up(:,i+1) = s2(:,i) + Ms2 .* (y2(:,i) - Cs2 * s2(:,i));
    
    %error calculation
    error1(:,i+1) = Al1 * error1(:,i) + [E;0] * w(i) - [L1;Ls1] * n1(i);
    ex1(:,i+1) = error1(1:2);
    es1(i+1) = error1(3);
    error2(:,i+1) = Al2 * error2(:,i) + [E;0] * w(i) - [L2;Ls2] * n2(i);
    ex2(:,i+1) = error2(1:2);
    es2(i+1) = error2(3);
end

%% Plotting
figure(2)
subplot(411)
plot(t,uref,t,u,'--','linewidth',2)
grid on
ylabel('amplitude')
legend('u_{ref}','u')

subplot(412)
scatter(t(1:length(t)-1),flag(1:length(t)-1))
grid on
ylabel('switching law')

subplot(413)
plot(t,y1(1,:),t,s1(1,:),'--','linewidth',2)
grid on
ylabel('amplitude')
legend('x_{s1}','s_1')

subplot(414)
plot(t,y2(1,:),t,s2(1,:),'--','linewidth',2)
grid on
xlabel('t / s')
ylabel('amplitude')
legend('x_{s2}','s_2')