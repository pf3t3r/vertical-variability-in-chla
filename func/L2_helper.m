function [ax,X_out,pSubml,bA,ks,obs,sk,ku,rV,p,pV,ad,tr,vuongRes,cSubml,pB] = L2_helper(tmp,pMaxMld,dcm,threshold,testSel,hypTest,yLimits,yLimitsObs,season,suppressFig)
%%L2_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L2 (sub-mixed layer region that
% is centred on the DCM). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% pMaxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% dcm: pressure of deep chlorophyll maximum (by cruise)
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.

logAxis = true; % true => output p-values in log x-axis, otherwise no log plot.
showLegends  = true;

% 1. Assign default parameter values.
% Default threshold of 50 based on Mishra et al (2019), Ghasemi & Zahediasl 
% (2012), and Ahad et al (2011).
if nargin < 4
    threshold = 50;
end
% Default to test against four distributions in K-S or A-D.
if nargin < 5
    testSel = 4;
end
% Default = "ks" (alternative = "ad").
if nargin < 6
    hypTest = "ks";
end
% Default y-limits for K-S/A-D and Vuong plots.
if nargin < 7
    yLimits = [-60 60];
end
% Default y-limits for horizontal bar chart showing no. of observations per
% depth. This should vary in a similar way to yLimits above. Test run the
% code to output "tr" to check this.
if nargin < 8
    yLimitsObs = [7 19];
end
if nargin < 9
    season = 0;
    % 0 = no seasonal analysis, 1 = winter, 2 = spring, 3 = summer,
    % 4 = autumn
end
if nargin < 10
    suppressFig = false;
end

id = num2str(tmp.data(:,1));
p = tmp.data(:,4);
c = tmp.data(:,5);
% clear tmp;

idOld = id;
%%% seasonal analysis
if season ~= 0
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
    elseif season == 2
        c = c(sprIds);
        p = p(sprIds);
        id = id(sprIds,:);
    elseif season == 3
        c = c(sumIds);
        p = p(sumIds);
        id = id(sumIds,:);
    elseif season == 4
        c = c(autIds);
        p = p(autIds);
        id = id(autIds,:);
    end
end

% Extract data beneath ML
% Exclude data with no measurements (-9 here => NaN)
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

% Extract measurements below pMaxMld
L = stop-1; % No. of casts with chl measurements in cruise 1 - 329
tmpP_subML = nan(L,1);
tmpCRN_subML = nan(L,1);
tmpX_subML = nan(L,1);
botID = [];

% To find the measurements taken beneath pMaxMld, we populate the arrays
% just defined for cases where p > pMld...
for i = 1:L
    tmpMld = pMaxMld(crn(i));
    if p(i) > tmpMld
        tmpP_subML(i) = p(i);
        tmpCRN_subML(i) = crn(i);
        tmpX_subML(i) = X1(i);
        tmpStr = id(i,:);
        botID = [botID;tmpStr];
    end
end

% ...and remove nan values (which represent measurements above pMaxMld)
pSubml = tmpP_subML(~isnan(tmpP_subML));
cSubml = tmpX_subML(~isnan(tmpX_subML));
idSubml = botID;


% Calculate KS p-value, skewness, kurtosis
% ..., centre around DCM (?)
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

% THIS is where we convert to Lagrangian pressure coordinates!!!
tPLagrangian = nan(length(p),1);
tPLagrangian = bottleArray(:,3) - bottleArray(:,4);
bottleArray = [bottleArray tPLagrangian];

pB10 = round(bottleArray(:,5),-1);
bottleArray = [bottleArray pB10];

bottleArray = [bottleArray cSubml];
tmin = min(bottleArray(:,6));
tmax = max(bottleArray(:,6));
tr = tmin:10:tmax;

X_out = bottleArray(:,7);
pB = bottleArray(:,6);
ks = nan(5,length(tr));
ad = nan(4,length(tr));
rV = nan(10,length(tr));
pV = nan(10,length(tr));
obs = nan(1,length(tr));

%sk = nan(1,length(tr));

for i = 1:length(tr)
    tmp = X_out(pB==tr(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        gammaParams = mle(tmp,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
        [~,ad(3,i)] = adtest(tmp,"Distribution","weibull");
        [~,ad(4,i)] = adtest(tmp,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        %sk(i) = skewness(tmp);
        %ku(i) = kurtosis(tmp);
    end
end

for i = 1:length(tr)
    if obs(i) < threshold
        ks(:,i) = nan;
        ad(:,i) = nan;
        rV(:,i) = nan;
        %sk(i) = nan;
        %ku(i) = nan;
    end
end

bA = bottleArray;


% Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = zeros(1,length(tr));
rV(isnan(rV)) = 0;

if testSel==4
    for i = 1:length(tr)
        %disp(i);
        if rV(1,i) & rV(2,i) & rV(3,i) > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
            vuongRes(i) = 2;
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
            vuongRes(i) = 3;
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
            vuongRes(i) = 4;
        end
    end
    rV(rV==0) = nan;
elseif testSel == 2
    for i = 1:length(tr)
        if rV(1,i)  > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
        end
    end
    rV(rV==0) = nan;
end

limits = yLimits;
obsId = [yLimitsObs(1) yLimitsObs(2)];

% Plot results
if suppressFig == true
    ax = nan;
else
    tix = limits(1):10:limits(2);
    disp(tix);
    a = obsId(1); b = obsId(2);
    n = length(tr);
    
    alphaHy = 0.005;
    alphaLlr = 0.1;
    
    vuongRes2 = nan(length(vuongRes),1);
    if b > length(vuongRes)
        b = length(vuongRes);
    end
    vuongRes2(a:b) = vuongRes(a:b);
    vuongRes = vuongRes2;


    ax = figure;

    subplot(1,3,3)
    barh(obs(a:b),'FaceColor','#d3d3d3');
    hold on
    xline(threshold);
    hold off
    set(gca,'YDir','reverse');
    yticklabels({});
    xlabel('No. of Obs.',Interpreter='latex',FontSize=13);
    ylim([1 1+b-a]);
    
    subplot(1,3,[1 2])
    xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
    hold on
    if strcmp(hypTest,"ks")
        if testSel == 4
            plot(ks(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ks(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
            plot(ks(3,:),tr,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
            plot(ks(4,:),tr,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
        elseif testSel == 2
            plot(ks(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ks(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        end
        xlabel('K-S $p$-value',Interpreter='latex',FontSize=15);
    else
        if testSel == 4
            plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
            plot(ad(3,:),tr,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
            plot(ad(4,:),tr,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
        elseif testSel == 2
            for i = 1:n
                if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr && ad(2,i) > alphaHy
                    plot(ad(1,i),tr(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
                elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr && ad(2,i) > alphaHy
                    plot(ad(1,i),tr(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
                elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr && ad(1,i) > alphaHy
                    plot(ad(2,i),tr(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
                elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr && ad(1,i) > alphaHy
                    plot(ad(2,i),tr(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
                end
            end
            plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        end
        xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
    end
    if logAxis == true
        set(gca, 'XScale', 'log');
    end
    hold off
    grid minor;
    ylim([limits(1) limits(2)]); xlim([0.1*alphaHy 1]);
    set(gca,'YDir','reverse');
    set(gca,"YTick",limits(1):10:limits(2));
    ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
    if showLegends == true
        legend('Location','best',FontSize=11);
    end
end
end