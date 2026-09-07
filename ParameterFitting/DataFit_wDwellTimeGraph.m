%% Load force length data
forceLengthData = readmatrix('ForceLength.dat');
delF = forceLengthData(2,1)-forceLengthData(1,1);
clrs = linspecer(4);
%% Plot the force length data

hold on
box on
set (gca, 'fontsize',15)
xlabel("Force (pN)")
ylabel("\langle\it{L}\rangle \rm(nm)")
xlim([0.07,5])

plot(forceLengthData(:,1),forceLengthData(:,2),"LineWidth",6,"Color",clrs(1,:))

%% Fit with 4 parameters

frwdSteps=csvread("NoOfForwardSteps.csv");
bcwdSteps=csvread("NoOfReverseSteps.csv");
dwellTimeFrwd=csvread("DwellTimeForward.csv");
dwellTimeRevr=csvread("DwellTimeReverse.csv");

deltaR = 22;
kBT = 4.114;
x1 = frwdSteps(:,1);
y1 = (frwdSteps(:,2)/(2.4*16.6+1.02*9.95)); %Convert to per second
x2 = bcwdSteps(:,1);
y2 = (bcwdSteps(:,2)/(2.4*16.6+1.02*9.95)); %Convert to per second
x3 = dwellTimeRevr(1:11,1); %Average forward and backward dwell times at 0.4 pN
delX = mean(diff(x3));
y3 = (dwellTimeRevr(1:11,2)./(sum(dwellTimeRevr(:,2))*delX));
x4 = dwellTimeFrwd(1:11,1);
delX = mean(diff(x4));
y4 = (dwellTimeFrwd(1:11,2)./sum(dwellTimeFrwd(:,2))*delX);
fForce = 0.4;
xdata=[x1; x2; x3; x4];
ydata=[y1; y2; (y3).^0.1; (y4).^0.1]; 
A = eye(4)*-1;
b = [0;0;0;-1];
%*(fit(1)*exp((1-fit(3)*f*deltaR))+fit(2)*exp(3*deltaR/(2*b)))



modelfun = @(fit,xdata) [((k2Plus(fit(3),x1,forceLengthData,delF)+k1Minus(fit(2),fit(4),x1,deltaR,kBT))./...
    (k2Plus(fit(3),x1,forceLengthData,delF)+k1Minus(fit(2),fit(4),x1,deltaR,kBT)+k1Plus(fit(1),fit(4),x1,deltaR,kBT)).*k1Plus(fit(1),fit(4),x1,deltaR,kBT));...
    ((k1Plus(fit(1),fit(4),x2,deltaR,kBT))./...
    (k2Plus(fit(3),x2,forceLengthData,delF)+k1Minus(fit(2),fit(4),x2,deltaR,kBT)+k1Plus(fit(1),fit(4),x2,deltaR,kBT)).*k1Minus(fit(2),fit(4),x2,deltaR,kBT));...
    ((k1Minus(fit(2),fit(4),fForce,deltaR,kBT)+k2Plus(fit(3),fForce,forceLengthData,delF)).*exp(-(k1Minus(fit(2),fit(4),fForce,deltaR,kBT)+k2Plus(fit(3),fForce,forceLengthData,delF)).*x3)).^0.1;...
    (k2Plus(fit(3),fForce,forceLengthData,delF)./(k2Plus(fit(3),fForce,forceLengthData,delF)+k1Minus(fit(2),fit(4),fForce,deltaR,kBT)).*...
    dwFwd1(k1Plus(fit(1),fit(4),fForce,deltaR,kBT),k1Minus(fit(2),fit(4),fForce,deltaR,kBT),k2Plus(fit(3),fForce,forceLengthData,delF),x4)...
    +...
    k1Minus(fit(2),fit(4),fForce,deltaR,kBT)/(k2Plus(fit(3),fForce,forceLengthData,delF)+k1Minus(fit(2),fit(4),fForce,deltaR,kBT))*...
    k1Plus(fit(1),fit(4),fForce,deltaR,kBT).*exp(-k1Plus(fit(1),fit(4),fForce,deltaR,kBT).*x4)).^0.1];
bestResnorm = Inf;
bestParams = [];
bestOutput = [];
nFits = 20;
for i=1:nFits
    [params, resnorm, residual, exitflag, output] = lsqcurvefit(modelfun,[rand,rand,rand,rand],xdata,ydata);
    if resnorm<bestResnorm && isreal(params) &&  all(params(:) > 0)
        bestParams = params;
        bestOutput = output;
    end
end


fit = bestParams

%% Plot with data

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

xForce =(0.1:0.05:5)';
yForce = (k2Plus(fit(3),xForce,forceLengthData,delF)+k1Minus(fit(2),fit(4),xForce,deltaR,kBT))./...
    (k2Plus(fit(3),xForce,forceLengthData,delF)+k1Minus(fit(2),fit(4),xForce,deltaR,kBT)+k1Plus(fit(1),fit(4),xForce,deltaR,kBT)).*k1Plus(fit(1),fit(4),xForce,deltaR,kBT);
%yForce = fit(1)*fit(4)*exp((1-2*fit(2))*xForce*deltaR/kBT);
yForce
hold on
box on

xlabel("Force (pN)")
ylabel(["Number of forward steps","in 1 second"])
set (gca, 'fontsize',15)
set(gca, "LineWidth",4)
plot(xForce,yForce,"LineWidth",6,"Color",clrs(1,:))
scatter(x1,y1,100,"filled","Color",clrs(2,:))

%% Plot with data

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

xForce = (0.1:0.05:5)';
%yForce = fit(1)*fit(4)*exp((1-2*fit(2))*xForce*deltaR/kBT)+fit(1)*fit(3)*exp(-fit(2)*xForce*deltaR/kBT*3*deltaR/(2*b));
yForce = (k1Plus(fit(1),fit(4),xForce,deltaR,kBT))./...
    (k2Plus(fit(3),xForce,forceLengthData,delF)+k1Minus(fit(2),fit(4),xForce,deltaR,kBT)+k1Plus(fit(1),fit(4),xForce,deltaR,kBT)).*k1Minus(fit(2),fit(4),xForce,deltaR,kBT);
hold on
box on
set (gca, 'fontsize',15)
set(gca, "LineWidth",4)
xlabel("Force (pN)")
ylabel(["No. of reverse steps","per second"])
plot(xForce,yForce,"LineWidth",6,"Color",clrs(1,:))
scatter(x2,y2,100,"filled","Color",clrs(2,:))

%% Plot dwell time of reverse step

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

dwellTimeFrwd=csvread("DwellTimeForward.csv");
dwellTimeRevr=csvread("DwellTimeReverse.csv");

xForce = 0.4
k1plus = fit(1)*exp(-fit(4)*xForce*deltaR/kBT);
k2plus = fit(3)*exp(0.58*(SampleLofF(forceLengthData,xForce,delF)-52));
k1minus = fit(2)*exp((1-fit(4))*xForce*deltaR/kBT)
t = 0:0.1:100;
prob = (k1minus+k2plus).*exp(-(k1minus+k2plus).*t);
hold on
box on
set(gca,"linewidth",4)
set (gca, 'fontsize',15)
xlabel("Time (s)")
ylabel("P(\tau_R)")
plot(t,prob,"LineWidth",6,"Color",clrs(1,:))
x3 = dwellTimeRevr(:,1); %Average forward and backward dwell times at 0.4 pN
y3 = dwellTimeRevr(:,2)./sum(dwellTimeRevr(:,1).*dwellTimeRevr(:,2));
scatter(x3,y3,100,"filled","Color",clrs(2,:))
%% Plot dwell time of forward step

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

dwellTimeFrwd=csvread("DwellTimeForward.csv");
dwellTimeRevr=csvread("DwellTimeReverse.csv");

xForce = 0.4;
k1plus = fit(1)*exp(-fit(4)*xForce*deltaR/kBT);
k2plus = fit(3)*exp(0.58*(SampleLofF(forceLengthData,xForce,delF)-60));
k1minus = fit(2)*exp((1-fit(4))*xForce*deltaR/kBT);
t = 0:0.1:100;
prob = (k1plus*(k1minus+k2plus)/(k1minus-k1plus+k2plus)*(exp(-k1plus*t)-exp(-(k1minus+k2plus)*t)))*(k2plus)/(k2plus+k1minus)+(k1minus)/(k2plus+k1minus)*k1plus*exp(-k1plus*t);
hold on
box on
set(gca,"linewidth",4)
set (gca, 'fontsize',15)
xlabel("Time (s)")
ylabel("P(\tau_F)")
plot(t,prob,"LineWidth",6,"Color",clrs(1,:))
x4 = dwellTimeFrwd(:,1);
delX = mean(diff(x4));
y4 = dwellTimeFrwd(:,2)./(sum(dwellTimeFrwd(:,2))*delX);
scatter(x4,y4,100,"filled","Color",clrs(2,:))
%% Plot lifetime of II

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

xForce = (0.1:0.01:5)';
%yForce = fit(1)*fit(4)*exp((1-2*fit(2))*xForce*deltaR/kBT)+fit(1)*fit(3)*exp(-fit(2)*xForce*deltaR/kBT*3*deltaR/(2*b));
yForce = 1./(k1Minus(fit(2),fit(4),xForce,deltaR,kBT)+k2Plus(fit(3),xForce,forceLengthData,delF));
hold on
box on
set(gca,"linewidth",4)
set (gca, 'fontsize',15)
xlabel("Force (pN)")
ylabel(["\tau_{X_1} (s)"])
xlim([0,5])
plot(xForce,yForce,"LineWidth",6,"Color",clrs(1,:))

%% Plot stationary probability of II
clrs = linspecer(5);

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

xForce = (0.1:0.01:5)';
%yForce = fit(1)*fit(4)*exp((1-2*fit(2))*xForce*deltaR/kBT)+fit(1)*fit(3)*exp(-fit(2)*xForce*deltaR/kBT*3*deltaR/(2*b));
yForce = (k1Minus(fit(2),fit(4),xForce,deltaR,kBT)+k2Plus(fit(3),xForce,forceLengthData,delF))./(k1Plus(fit(1),fit(4),xForce,deltaR,kBT)+k1Minus(fit(2),fit(4),xForce,deltaR,kBT)+k2Plus(fit(3),xForce,forceLengthData,delF));
hold on
box on
set(gca,"linewidth",4)
set (gca, 'fontsize',15)
xlabel("Force (pN)")
ylabel("Stationary Probability")
xlim([0,5])
plot(xForce,yForce,"LineWidth",6,"Color",clrs(1,:))
plot(xForce,1-yForce,"LineWidth",6,"Color",clrs(2,:))
legend("\pi_{X_0}","\pi_{X_1}")

%% Plot k2p*piII

xForce = (0.1:0.01:5)';
%yForce = fit(1)*fit(4)*exp((1-2*fit(2))*xForce*deltaR/kBT)+fit(1)*fit(3)*exp(-fit(2)*xForce*deltaR/kBT*3*deltaR/(2*b));
yForce = k2Plus(fit(3),xForce,forceLengthData,delF).*k1Plus(fit(1),fit(4),xForce,deltaR,kBT)./(k1Plus(fit(1),fit(4),xForce,deltaR,kBT)+k1Minus(fit(2),fit(4),xForce,deltaR,kBT)+k2Plus(fit(3),xForce,forceLengthData,delF));
hold on
box on
set (gca, 'fontsize',15)
xlabel("Force (pN)")
ylabel(["k_2^+\cdot\pi_{II}"])
xlim([0,5])
plot(xForce,yForce,"LineWidth",6)


%% Plot k_2^+/k_1^-

xForce = 0:0.1:5;
yForce = ((fit(1)*fit(2)*exp((1-2*fit(4))*xForce*deltaR/kBT))...
          ./(fit(3)+fit(1)*exp(-fit(4)*xForce*deltaR/kBT)+fit(2)*exp((1-fit(4))*xForce*deltaR/kBT)))./...
          ((fit(3)*fit(1)*exp(-fit(4)*xForce*deltaR/kBT)+fit(2)*fit(1)*exp((1-2*fit(4))*xForce*deltaR/kBT))...
          ./(fit(3)+fit(1)*exp(-fit(4)*xForce*deltaR/kBT)+fit(2)*exp((1-fit(4))*xForce*deltaR/kBT))) ;
hold on
box on

set (gca, 'fontsize',15)
xlabel("Force (pN)")
ylabel("l/m_1")
plot(xForce,(yForce),"LineWidth",6)

data = load("CondensinAverageBackstepData.dat");
scatter(x1,y2./y1,100,"filled")


%% Plot many dwell times together - Reverse Steps
clrs = linspecer(5);

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

forceArr = [0.1, 0.2, 0.3, 0.4, 0.5];
hold on
box on
set(gca,"linewidth",4)
set (gca, 'fontsize',15)
xlabel("time (s)")
ylim([0,0.25])
ylabel("P(\tau_R)")
i=1;
for xForce = forceArr
    k1plus = fit(1)*exp(-fit(4)*xForce*deltaR/kBT);
    k2plus = fit(3)*exp(0.58*(SampleLofF(forceLengthData,xForce,delF)-52));
    k1minus = fit(2)*exp((1-fit(4))*xForce*deltaR/kBT)
    t = 0:0.1:100;
    prob = (k1minus+k2plus).*exp(-(k1minus+k2plus).*t);
    
    plot(t,prob,"LineWidth",6,"Color",clrs(i,:))
    i = i+1;
end
legend("F=0.1 pN","F=0.2 pN","F=0.3 pN","F=0.4 pN","F=0.5 pN")


%% Plot many dwell times together - Forward Steps
clrs = linspecer(5);

x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

forceArr = [0.1, 0.2, 0.3, 0.4, 0.5];
hold on
box on
set(gca,"linewidth",4)
set (gca, 'fontsize',15)
xlabel("time (s)")
%ylim([0,0.25])
ylabel("P(\tau_F)")
i=1;
for xForce = forceArr
    k1plus = fit(1)*exp(-fit(4)*xForce*deltaR/kBT);
    k2plus = fit(3)*exp(0.58*(SampleLofF(forceLengthData,xForce,delF)-60));
    k1minus = fit(2)*exp((1-fit(4))*xForce*deltaR/kBT);
    t = 0:0.1:100;
    prob = (k1plus*(k1minus+k2plus)/(k1minus-k1plus+k2plus)*(exp(-k1plus*t)-exp(-(k1minus+k2plus)*t)))*(k2plus)/(k2plus+k1minus)+(k1minus)/(k2plus+k1minus)*k1plus*exp(-k1plus*t);

    plot(t,prob,"LineWidth",6,"Color",clrs(i,:))
    i = i+1;
end
legend("F=0.1 pN","F=0.2 pN","F=0.3 pN","F=0.4 pN","F=0.5 pN")