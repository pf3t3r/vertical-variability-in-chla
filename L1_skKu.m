% Skewness-kurtosis analysis of the mixed layer (L1) at Station ALOHA for chl-a.
% We import the data, process it, and calculate the skewness and kurtosis.

clear; clc; close all;
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

% Options and Test Cases.
% NOTE that there is no option to set the y-limits here. To see results for
% a deeper mixed layer you can instead adjust the threshold (but do not
% trust values below 30).
threshold = 30; % default threshold for statistical tests
testSel = 2;    % 2 = norm + logn; 4 = norm + logn + weib + gamm
season = 0;     % if zero, run default analysis. otherwise run seasonal
                % analysis (1 = winter (DJF), 2 = spring, 3 = summer, 4 =
                % autumn)

% Load Data: chl-a and mixed layer depth 'mld'
tmp = importdata("data/hot_chla.txt");
mld = load("output/mldVals.mat").maxMld;

pIn = tmp.data(:,4);        % pressure 
cIn = tmp.data(:,5);        % variable to test (chl-a)
idIn = tmp.data(:,1);       % bottle ID, needed for seasonal analysis


%% Seasonal analysis code.
if season == 0
    tmpT = ": L1";
elseif season ~= 0
    botidIn = tmp.data(:,2);
    n = length(pIn);
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
        cIn = cIn(winIds);
        pIn = pIn(winIds);
        idIn = idIn(winIds);
        tmpT = ": L1 Winter";
    elseif season == 2
        cIn = cIn(sprIds);
        pIn = pIn(sprIds);
        idIn = idIn(sprIds);
        tmpT = ": L1 Spring";
    elseif season == 3
        cIn = cIn(sumIds);
        pIn = pIn(sumIds);
        idIn = idIn(sumIds);
        tmpT = ": L1 Summer";
    elseif season == 4
        cIn = cIn(autIds);
        pIn = pIn(autIds);
        idIn = idIn(autIds);
        tmpT = ": L1 Autumn";
    end
end

%% Extract data in mixed layer.

% Extract the cruise number 'crn'.
tmp = num2str(idIn);
crn = str2num(tmp(:,1:3)); clear tmp;

L = length(pIn);    % Length of dataset

% Initialise ...
tmpP = nan(1,L);
tmpCrn = nan(1,L);
tmpX = nan(1,L);
tmpId = nan(1,L);

% Stop at cruise number 329. This was in order to match the latest data
% available for CTD fluorometry but I could consider updating this.
for i = 1:L
    if crn(i) == 330
        stop = i;
        break
    elseif crn(i) > 330
        stop = i;
        break
    else
        stop = L+1;
    end
end

% Save values at pressures above (i.e. shallower than) the mixed layer depth
for i = 1:stop-1
    tmp = mld(crn(i));
    if pIn(i) < tmp
        tmpP(i) = pIn(i);
        tmpCrn(i) = crn(i);
        tmpX(i) = cIn(i);
        tmpId(i) = idIn(i);
    end
end

% Remove NaNs from new arrays of values in mixed layer.
pOut = tmpP(~isnan(tmpP));
crnOut = tmpCrn(~isnan(tmpCrn));
cOut = tmpX(~isnan(tmpX));
idOut = tmpId(~isnan(tmpId));

%% Data: Quality Control and Bin.

% Remove bottles that are too close to the surface (< 2.5 dbar)
idRm = pOut > 2.5;
pOut = pOut(idRm);
cOut = cOut(idRm);
botid = idOut(idRm)';

% Remove bottles where concentration of X = 0
idZero = cOut == 0;
pOut = pOut(~idZero);
cOut = cOut(~idZero);
botid = botid(~idZero);

% Save cruise number (CRN) of each bottle - needed below
tmp = num2str(botid);
crn = str2num(tmp(:,1:3)); clear tmp;

% Remove bottles from cruises 330 on (b/c fluorescence analysis not done)
for i = 1:length(pOut)
    if crn(i) > 329
        id329 = i - 1;
        break;
    else
        id329 = length(pOut);
    end
end

pOut = pOut(1:id329);
cOut = cOut(1:id329);
clear idRm idZero id329 i;

pb10 = discretize(pOut,0:10:200);
n10 = max(pb10);

cOutB = cOut;
pOutB = pb10;

%% Calculate skewness and kurtosis.
% Calculate the skewness and kurtosis of the cleaned data. NOTE that other
% parameters related to A-D calculation are still present. These are to be
% removed.

n = 15;
depth = 5:10:150;
obs = nan(n,1);
sk = nan(n,1);
ku = nan(n,1);
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = cOutB(pOutB==i);
    % apply statistical tests to the data   
    if length(X_i) > 3
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

% Remove values that do not meet the threshold
for i = 1:n
    if obs(i) < threshold
        sk(i) = nan;
        ku(i) = nan;
    end
end

% Remove these nan values
tmp = [];
for i = 1:n
    if ~isnan(sk(i))
        tmp = [tmp i];
    end
end
p = depth(tmp);
sk = sk(tmp);
ku = ku(tmp);

n = length(p);


%% Create figure.

% First draw theoretical curves for the various distribution families.

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
skLognN = -skLogn;
kuLognN = kuLogn;

if testSel == 4
    % Gamma family: generate theoretical skewness and kurtosis
    kTh = linspace(0.04,3000,1500000);
    for i = 1:length(kTh)
        skGam(i) = 2/sqrt(kTh(i));
        kuGam(i) = 6/kTh(i) + 3;
    end
    skGamN= -skGam;
    kuGamN = kuGam;
    
    % Weibull family: generate theoretical skewness and kurtosis
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
if max(ku) > 10 & min(sk) < 0
    kurtLimB = max(ku) + 1;
    skewLimA = min(sk) - 0.1;
    skewLimB = max(sk) + 0.1;
elseif max(ku) > 10
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
elseif min(sk) < 0 
    skewLimA = min(sk) - 0.1;
elseif max(sk) > 2.5
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

% Plot results.
ax = figure;
scatter(nan,nan,72,[0.8 0.8 0.8],DisplayName='Data');
hold on
scatter(0,3,72,[0.7725490196078432 0.10588235294117647 0.49019607843137253],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
if testSel == 2
    plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
    plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
elseif testSel == 4
    plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#4d9221',LineStyle='-',LineWidth=1.7);
    plot(skLognN,kuLognN,'Color','#4d9221',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
    plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#d3d3d3',LineStyle='-',LineWidth=1.7);
    plot(skWblN,kuWblN,'Color','#d3d3d3',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
    plot(skGam,kuGam,'DisplayName','Gam.','Color','#000000',LineStyle='-',LineWidth=1.7);
    plot(skGamN,kuGamN,'Color','#000000',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
end
scatter(sk,ku,72,[0.8 0.8 0.8],HandleVisibility="off");
if length(p) > 3
    clr = 1:1:length(p);
    scatter(sk,ku,54,clr,"filled","o",HandleVisibility="off");
    colormap(gca,cbrewer2("Greens"));
    cbar = colorbar;
    cbar.Ticks = 1:1:length(p);
    cbar.TickLabels = p;
    cbar.Label.Position = [0.7 1-0.35];
elseif length(p) == 3
    c = [0.9686    0.9882    0.9608; 0.4579    0.7700    0.4643;  0.0000    0.2667    0.1059];
    scatter(sk,ku,54,c,"filled","o",HandleVisibility="off");
    colormap(gca,c);
    cbar = colorbar;
    cbar.Ticks = [0.17 0.5 0.83];
    cbar.TickLabels = p;
    cbar.Label.Position = [0.7 -0.05];
elseif length(p) == 2
    c = [0.9686    0.9882    0.9608; 0.0000    0.2667    0.1059];
    scatter(sk,ku,54,c,"filled","o",HandleVisibility="off");
    colormap(gca,c);
    cbar = colorbar;
    cbar.Ticks = [0.25 0.75];
    cbar.TickLabels = p;
    cbar.Label.Position = [0.7 -0.05];
end
cbar.Direction = "reverse";
cbar.Label.String = "P [dbar]";
cbar.Label.Rotation = 0;
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness','FontSize',13,'Interpreter','latex'); ylabel('Kurtosis',FontSize=13,Interpreter='latex');
lgd = legend('Location','best');
title('chl-$a$'+tmpT,'Interpreter','latex','FontSize',13);