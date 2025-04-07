% Skewness-kurtosis analysis of the DCM layer (L2) at Station ALOHA for chl-a.
% We import the data, process it, and calculate the skewness and kurtosis.

clear; clc; close all;
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

% Options and test cases.
threshold = 10;
limits = [-60 60];  % show results for this range in pressure values
                    % if the range is too large, the script will update to
                    % the nearest possible value.
testSel = 2;        % 2 = norm + logn; 4 = norm + logn + weib + gamm
season = 0;         % if zero, run default analysis. otherwise run seasonal
                    % analysis (1 = winter, 2 = spring, 3 = summer, 4 =
                    % autumn). NOTE that the seasonal analysis will skip
                    % certain depths due to insufficient measurements.
showYLabel = true;  % turn off for paper.
showLegend = true;  % turn off for paper.

% Data Import.
tmp = importdata("data/hot_chla.txt");
id = num2str(tmp.data(:,1));
p = tmp.data(:,4);
c = tmp.data(:,5);

%% Seasonal analysis.
if season == 0
    tmpT = ": L2";
elseif season ~= 0
    botidIn = tmp.data(:,2);
    n = length(p);
    botidIn(botidIn==-9) = nan;
    botId2 = num2str(botidIn);
    botMth = nan(n,1);
    for i = 1:n
        tmpX = str2num(botId2(i,1:end-4));
        if ~isnan(tmpX)
            botMth(i) = tmpX;
        end
    end
    winter = nan(n,1); spring = nan(n,1); summer = nan(n,1); autumn = nan(n,1);
    for i = 1:n
        tmpY = botMth(i);
        if (tmpY == 12) || (tmpY == 1) || (tmpY == 2)
            winter(i) = 1;
        end
        if (tmpY == 3) || (tmpY == 4) || (tmpY == 5)
            spring(i) = 1;
        end
        if (tmpY == 6) || (tmpY == 7) || (tmpY == 8)
            summer(i) = 1;
        end
        if (tmpY == 9) || (tmpY == 10) || (tmpY == 11)
            autumn(i) = 1;
        end
    end

    winIds = []; sprIds = []; sumIds = []; autIds = []; 
    for i = 1:n
        if winter(i) == 1
            winIds = [winIds i];
        end
        if spring(i) == 1
            sprIds = [sprIds i];
        end
        if summer(i) == 1
            sumIds = [sumIds i];
        end
        if autumn(i) == 1
            autIds = [autIds i];
        end
    end

    if season == 1
        c = c(winIds);
        p = p(winIds);
        id = id(winIds,:);
        tmpT = ": L2 Winter";
    elseif season == 2
        c = c(sprIds);
        p = p(sprIds);
        id = id(sprIds,:);
        tmpT = ": L2 Spring";
    elseif season == 3
        c = c(sumIds);
        p = p(sumIds);
        id = id(sumIds,:);
        tmpT = ": L2 Summer";
    elseif season == 4
        c = c(autIds);
        p = p(autIds);
        id = id(autIds,:);
        tmpT = ": L2 Autumn";
    end
end

%% Extract data beneath ML
c(c==-9) = nan;
id = id(~isnan(c),:);
p = p(~isnan(c));
X1 = c(~isnan(c));

% Cruise number 'crn'
crn = str2num(id(:,1:3));

% Stop evaluating after crn = 329
for i = 1:length(crn)
    if crn(i) == 330
        stop = i;
        break
    elseif crn(i) > 330
        stop = i;
        break
    else
        stop = length(p) + 1;
    end
end

% Extract measurements below mixed layer depth 'mld'
L = stop-1; % No. of casts with chl measurements in cruise 1 - 329
tmpP_subML = nan(L,1);
tmpCRN_subML = nan(L,1);
tmpX_subML = nan(L,1);
botID = [];

% To find the measurements taken beneath mld, we populate the arrays
% just defined for cases where p > mld...
for i = 1:L
    tmpMld = mld(crn(i));
    if p(i) > tmpMld
        tmpP_subML(i) = p(i);
        tmpCRN_subML(i) = crn(i);
        tmpX_subML(i) = X1(i);
        tmpStr = id(i,:);
        botID = [botID;tmpStr];
    end
end

% ...and remove nan values (which represent measurements above mld)
pSubml = tmpP_subML(~isnan(tmpP_subML));
cSubml = tmpX_subML(~isnan(tmpX_subML));
idSubml = botID;

% Extract cruise number 'crn' and 'cast'
crn = str2num(idSubml(:,1:3)); 
cast = str2num(idSubml(:,6:8));
cast(cast==100) = nan;

bottleArray = [crn cast pSubml];

% Create an array of all unique bottle cruise/cast combinations
botCrnCast = rmmissing(unique(bottleArray(:,1:2),"rows"));

% Find these unique cruise/cast combinations in the 'dcm' array
dcmCrnCast = [];
for i = 1:length(dcm(:,1))
    for x = 1:length(botCrnCast)
        if dcm(i,1:2) == botCrnCast(x,1:2) 
            dcmCrnCast = [dcmCrnCast i];
        end
    end
end

% Split bottle concentration by cruise & cast
tid = [];
for i = 2:length(pSubml)
    % check if CRN or CAST changes
    if bottleArray(i,1) > bottleArray(i-1,1) || bottleArray(i,2) > bottleArray(i-1,2)
        tid = [tid i];
    end
end

tPcm = nan(length(pSubml),1);
tPcm(1:tid(1)-1) = dcm(dcmCrnCast(1),3);
tPcm(tid(end):end) = dcm(dcmCrnCast(end),3);

Ltid = length(tid);

tmp = length(dcmCrnCast) - 1;
Ltid = tmp;
for i = 2:Ltid-2
    tPcm(tid(i):tid(i+1)-1) = dcm(dcmCrnCast(i),3);
end

bottleArray = [bottleArray tPcm];

% Convert to DCM-centred coordinates.
tPLagrangian = nan(length(pSubml),1);
tPLagrangian = bottleArray(:,3) - bottleArray(:,4);
bottleArray = [bottleArray tPLagrangian];

%% Binning.
% Bin pressures into 10 dbar bins.
pB10 = round(bottleArray(:,5),-1);
bottleArray = [bottleArray pB10];

bottleArray = [bottleArray cSubml];
tmin = min(bottleArray(:,6));
tmax = max(bottleArray(:,6));
tr = tmin:10:tmax;

X_out = bottleArray(:,7);
pB = bottleArray(:,6);

%% Calculate skewness and kurtosis.
obs = nan(1,length(tr));
sk = nan(1,length(tr));
ku = nan(1,length(tr));

for i = 1:length(tr)
    tmp = X_out(pB==tr(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
    end
end

% Set values that do not meet the threshold as NaN.
for i = 1:length(tr)
    if obs(i) < threshold
        sk(i) = nan;
        ku(i) = nan;
    end
end

% Set limits. If the limits specified above are outside of the range of the
% data, this code will set the limits to the nearest possible value.
[~,idA] = (min(abs(tr-limits(1))));
[~,idB] = (min(abs(tr-limits(2))));
a = tr(idA);
b = tr(idB);

% Remove NaN values
tmp = [];
for i = idA:idB
    if ~isnan(sk(i))
        tmp = [tmp i];
    end
end
tr2 = tr(tmp);
sk2 = sk(tmp);
ku2 = ku(tmp);
clear tmp;

%% Prepare figure.
% Draw theoretical curves for each distribution family tested.
% Lognormal
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
skLognN = -skLogn;
kuLognN = kuLogn;

if testSel == 4
    % Gamma
    kTh = linspace(0.04,3000,1500000);
    for i = 1:length(kTh)
        skGam(i) = 2/sqrt(kTh(i));
        kuGam(i) = 6/kTh(i) + 3;
    end
    skGamN= -skGam;
    kuGamN = kuGam;

    % Weibull
    kWbl = linspace(0.1,3.5,10000);
    for i = 1:length(kWbl)
        skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
            ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
        kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
           ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
    end
    skWblN = -skWbl;
    kuWblN = kuWbl;
end

% Automatically adjust limits.
kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
if max(ku2) > 10 & min(sk2) < 0
    kurtLimB = max(ku2) + 1;
    skewLimA = min(sk2) - 0.1;
    skewLimB = max(sk2) + 0.1;
elseif max(ku2) > 10
    kurtLimB = max(ku2) + 1;
    skewLimB = max(sk2) + 0.1;
elseif min(sk2) < 0 
    skewLimA = min(sk2) - 0.1;
elseif max(sk2) > 2.5
    kurtLimB = max(ku2) + 1;
    skewLimB = max(sk2) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

%% Create figure.
ax = figure;
scatter(0,3,72,[0.7725490196078432 0.10588235294117647 0.49019607843137253],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#4d9221',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#4d9221',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
if testSel == 4
    plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#d3d3d3',LineStyle='-',LineWidth=1.7);
    plot(skWblN,kuWblN,'Color','#d3d3d3',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
    plot(skGam,kuGam,'DisplayName','Gam.','Color','#000000',LineStyle='-',LineWidth=1.7);
    plot(skGamN,kuGamN,'Color','#000000',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
end
scatter(sk2,ku2,72,[0.8 0.8 0.8],HandleVisibility='off');
clr = 1:1:length(tr2);
scatter(sk2,ku2,54,clr,"filled","o",HandleVisibility="off");
colormap(gca,flipud(cbrewer2("PiYG")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(tr2);
cbar.TickLabels = tr2;
cbar.Label.String = "P [dbar]";
cbar.Label.Position = [0.7 1-0.7];
cbar.Label.Rotation = 0;
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness',FontSize=13,Interpreter='latex'); 
if showYLabel == true
    ylabel('Kurtosis',FontSize=13,Interpreter='latex');
else
    yticklabels({});
end
if showLegend == true
    legend(Location="best");
end
title("chl-$a$"+tmpT,"Interpreter","latex",FontSize=13);