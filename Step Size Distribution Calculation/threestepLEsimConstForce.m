%% Load the cdf_inv of actual distribution

cdf_inv_pofl = load("pofL_cdf_inv.mat");
cdf_inv_pofl = cdf_inv_pofl.dummyVar;
delX = cdf_inv_pofl(2,1)-cdf_inv_pofl(1,1);

%% Calculate step size distributions

delT = 1; %s
tMax = 10000000; %s
nSteps = round(tMax/delT);

deltaR = 22;
kBT = 4.114;
%0.1636    0.0123    0.3290    0.4823
k1plus_0 = 0.16;
k1minus_0 = 0.012;
k2plus_0 = 0.329;
theta = 0.48;

currState = 1;
currRealLoop = 100; %setting it to 100 to keep things alright
currTransientLoop = 0;

forceLengthData = readmatrix('../ForceLength.dat');
delF = forceLengthData(2,1)-forceLengthData(1,1);
nForward = 0;
nTrueForward = 0; 
nReverse = 0;
currForce = 0.8; %pN
fileID = fopen('StepSizes/steps_'+string(currForce)+'.txt','w');

for i=1:nSteps

    %if currState == 1
    %randomly take k1plus
    %update currState = 2

    %if currState == 2
    %randomly take either k2plus or k1plus
    currT = i*delT;
    

    if currState == 1
        currK1plus = k1Plus(k1plus_0,theta,currForce,deltaR,kBT);
        pk1plus = exp(-currK1plus*delT);
        if rand>pk1plus %Transit to state 2
            currTransientLoop = SampleLofF(forceLengthData,currForce,delF);
            currState = 2;
            nForward = nForward+1;
            currStep = SamplePofL(cdf_inv_pofl,rand,delX)-14+normrnd(0,13); % /3.4*10/monomerSize
            fprintf(fileID,'%f\n',currStep);
        end
    elseif currState == 2
        currK1minus = k1Minus(k1minus_0,theta,currForce,deltaR,kBT);
        currK2plus = k2PlusFixedLength(k2plus_0,currTransientLoop);
        % disp([currK1minus,currK2plus])
        pLeave = exp(-(currK2plus+currK1minus)*delT);
        pk2plus = exp(-currK2plus*delT);
        pk1minus = exp(-currK1minus*delT);
        if rand>pLeave
           if rand<(currK2plus/(currK2plus+currK1minus)) %Transit to forward, true forward step
               currRealLoop = currRealLoop+currTransientLoop;
               currTransientLoop = 0;
               currState = 1;
               nTrueForward = nTrueForward+1;
           else
               currTransientLoop = 0;
               currState = 1;
               nReverse = nReverse+1;
               fprintf(fileID,'%f\n',-currStep);
           end
        end
    end
end
% disp([forceChangeFreq,nTrueForward/nForward,nReverse/nForward])
% nFreqArray(iFreq) = forceChangeFreq;
% nForwardArray(iFreq) = nForward;
% nTrueForwardArray(iFreq) = nTrueForward;
% nReverseArray(iFreq) = nReverse;

%% plot the data

hold on
box on
stepSizes = load('StepSizes/steps_'+string(currForce)+'.txt','-ascii');
set (gca, 'fontsize',15)
xlim([-200,200])
xlabel("Step sizes (nm)")
ylabel("Freq")

histogram(stepSizes(abs(stepSizes)<200),50,"LineWidth",4,"Normalization","probability",'DisplayStyle', 'stairs')


%% Plot all the data together
x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

hold on
box on
clrs = linspecer(4);
set(gca, 'linewidth', 4)
xlabel("S_t (nm)")
ylabel("P(S_t)")

currForce = 0.1;
stepSizes = load('StepSizes/steps_'+string(currForce)+'.txt','-ascii');
set (gca, 'fontsize',15)
xlim([-200,200])
histogram(stepSizes(abs(stepSizes)<200),100,"LineWidth",3,"Normalization","probability",'DisplayStyle', 'stairs',"EdgeColor",clrs(1,:));

currForce = 0.4;
stepSizes = load('StepSizes/steps_'+string(currForce)+'.txt','-ascii');
set (gca, 'fontsize',15)
xlim([-200,200])
hst = histogram(stepSizes(abs(stepSizes)<200),100,"LineWidth",3,"Normalization","probability",'DisplayStyle', 'stairs',"EdgeColor",clrs(2,:));

currForce = 0.8;
stepSizes = load('StepSizes/steps_'+string(currForce)+'.txt','-ascii');
set (gca, 'fontsize',15)
xlim([-200,200])
histogram(stepSizes(abs(stepSizes)<200),100,"LineWidth",3,"Normalization","probability",'DisplayStyle', 'stairs',"EdgeColor",clrs(3,:))

maxValBkwd = max(hst.Values(1:50));
maxValFwd = max(hst.Values(51:end));

frwdSteps=csvread("0p4FwdSteps.csv");
revSteps=csvread("0p4RevSteps.csv");
normFact = (maxValBkwd+maxValFwd)/(max(frwdSteps(:,2))+max(revSteps(:,2)));

histogram('BinEdges',frwdSteps(:,1)','BinCounts',normFact*frwdSteps(1:end-1,2),"LineWidth",4,"EdgeColor",clrs(4,:),'DisplayStyle','stairs')
histogram('BinEdges',revSteps(:,1)','BinCounts',normFact*revSteps(1:end-1,2),"LineWidth",4,"EdgeColor",clrs(4,:),'DisplayStyle','stairs')
histogram('BinEdges',frwdSteps(:,1)','BinCounts',normFact*frwdSteps(1:end-1,2),"LineWidth",4,"EdgeColor",clrs(4,:),"EdgeAlpha",0,"FaceColor",clrs(4,:),"FaceAlpha",0.5)
histogram('BinEdges',revSteps(:,1)','BinCounts',normFact*revSteps(1:end-1,2),"LineWidth",4,"EdgeColor",clrs(4,:),"EdgeAlpha",0,"FaceColor",clrs(4,:),"FaceAlpha",0.5)
% plot(-frwdSteps(:,1),frwdSteps(:,2)*normFact,"LineWidth",4,"Color",clrs(2,:),"LineStyle",":")
% plot(-revSteps(:,1),revSteps(:,2)*normFact,"LineWidth",4,"Color",clrs(2,:),"LineStyle",":")

legend("F = 0.1 pN","F = 0.4 pN","F = 0.8 pN","Experiment, F = 0.4 pN")


%% Force-velocity computation

nSteps = 100000;

deltaR = 22;
kBT = 4.114;

k1plus_0 = 0.16;
k1minus_0 = 0.012;
k2plus_0 = 0.329;
theta = 0.48;

currState = 1;
currRealLoop = 100; %setting it to 100 to keep things alright
currTransientLoop = 0;

forceLengthData = readmatrix('../ForceLength.dat');
delF = forceLengthData(2,1)-forceLengthData(1,1);
nForward = 0;
nTrueForward = 0; 
nReverse = 0;
forceArr = 0.1:0.1:2.0; %pN
velocityArr = forceArr;
for iForce=1:length(forceArr)
    currForce = forceArr(iForce)
    fileID = fopen('StepSizes/steps_w_time_'+string(currForce)+'.txt','w');
    currT = 0;
    timePassed = 0;
    totalLength = 0;
    for i=1:nSteps
    
        %if currState == 1
        %randomly take k1plus
        %update currState = 2
    
        %if currState == 2
        %randomly take either k2plus or k1plus
        
        
    
        if currState == 1
            currK1plus = k1Plus(k1plus_0,theta,currForce,deltaR,kBT);
            tau = exprnd(1/currK1plus); %fix: time till next step
            timePassed = timePassed+tau;
            %Transit to state 2
            currTransientLoop = SampleLofF(forceLengthData,currForce,delF);
            currState = 2;
            nForward = nForward+1;
            currStep = SamplePofL(cdf_inv_pofl,rand,delX)-14+normrnd(0,13); % /3.4*10/monomerSize
            totalLength = totalLength+currStep;
            fprintf(fileID,'%f\t%f\n',timePassed,currStep);
            
        elseif currState == 2
            currK1minus = k1Minus(k1minus_0,theta,currForce,deltaR,kBT);
            currK2plus = k2PlusFixedLength(k2plus_0,currTransientLoop);
            tau = exprnd(1/(currK2plus+currK1minus));
            timePassed = timePassed+tau;
            % pk2plus = exp(-currK2plus*delT);
            % pk1minus = exp(-currK1minus*delT);
            if rand<(currK2plus/(currK2plus+currK1minus)) %Transit to forward, true forward step
               currRealLoop = currRealLoop+currTransientLoop;
               currTransientLoop = 0;
               currState = 1;
               nTrueForward = nTrueForward+1;
            else
               currTransientLoop = 0;
               currState = 1;
               nReverse = nReverse+1;
               totalLength = totalLength-currStep;
               fprintf(fileID,'%f\t%f\n',timePassed,-currStep);
            end
            
        end
    end
    velocity = totalLength/timePassed;
    velocityArr(iForce) = velocity;
end
% disp([forceChangeFreq,nTrueForward/nForward,nReverse/nForward])
% nFreqArray(iFreq) = forceChangeFreq;
% nForwardArray(iFreq) = nForward;
% nTrueForwardArray(iFreq) = nTrueForward;
% nReverseArray(iFreq) = nReverse;
hold on
box on
set(gca,"linewidth",4)
x      = 0;   % Screen position
y      = 0;   % Screen position
width  = 900; % Width of figure
height = 700; % Height of figure (by default in pixels)

figure('Position', [x y width height]);

clrs = linspecer(4);
set (gca, 'fontsize',15)
xlabel("Force (pN)")
ylabel("Velocity (nm/s)")

plot((forceArr),(velocityArr),"LineWidth",6,"Color",clrs(1,:))
