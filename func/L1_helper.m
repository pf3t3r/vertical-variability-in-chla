function [ax,p,ks,obs,rV,pV,ad,cOut,pOut] = L1_helper(D,maxMld,threshold,testSel,hypTest,logAxis,season,suppressFig)
%L1_helper Calculate statistical tests for chl-a and other variables in the mixed layer (L1).
% Specifically, we calculate the hypothesis test
% results (A-D or K-S) here. Additionally, we also calculate the Vuong
% test, which tells us the best distribution when A-D or K-S predicts that
% multiple distributions fit the data at a given depth.

% INPUTS
% D: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% threshold: the statistical test (A-D) will only run above this threshold.
% testSel: choose between 2 (normal-lognormal) or 4
% (normal-lognormal-weibull-gamma).
% hypTest: choose "ad" to run the Anderson-Darling (A-D) test or "ks" to
% run the Kolmogorov-Smirnov (K-S) test.
% logAxis: if set to true, outputs the x-axis (p-values) in log scale.
% season: 0 = no seasonal analysis; otherwise 1-4 => run seasonal analysis
% on winter (1), spring (2), summer (3), or autumn (4).
% suppressFig: if set to true, the figure does not output.

% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.
% rV = Vuong Log-Likelihood Ratio. If positive, the first of two tested
% distributions is best; if negative, the second is best. 
% pV = p-value associated with the Vuong Log-Likelihood Ratio. If less than
% 0.1 (paper ref here), we deem this significant.
% ad = A-D test p-values at pressures p
% cOut = output cleaned variable
% pOut = output cleaned pressures

% Select expectation values as well as y-axis limits.
% NOTE: binLimA and binLimB -> limits of y-axis in terms of bins. 
% limits -> limits of y-axis in dbar.
alphaHy = 0.005;
alphaLlr = 0.10;
binLimA = 0.5;
binLimB = 10.5;
limits = [0 100];

if nargin < 8
    suppressFig = false;
end
if nargin < 7
    season = 0;
    % 0 = no seasonal analysis, 1 = winter, 2 = spring, 3 = summer,
    % 4 = autumn
end
if nargin < 6
    logAxis = true;
end
if nargin < 5
    hypTest = "ad";
end
if nargin < 4
    testSel = 2;
end
if nargin < 3
    threshold = 30;
end

pIn = D.data(:,4);
cIn = D.data(:,5);
idIn = D.data(:,1);

% Seasonal Analysis.
if season ~= 0
    botidIn = D.data(:,2);
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
    elseif season == 2
        cIn = cIn(sprIds);
        pIn = pIn(sprIds);
        idIn = idIn(sprIds);
    elseif season == 3
        cIn = cIn(sumIds);
        pIn = pIn(sumIds);
        idIn = idIn(sumIds);
    elseif season == 4
        cIn = cIn(autIds);
        pIn = pIn(autIds);
        idIn = idIn(autIds);
    end
end

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
    tmp = maxMld(crn(i));
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

% Apply statistical tests to the cleaned data.
% Specifically, we apply the Kolmogorov-Smirnov (K-S) and Anderson-Darling
% (A-D) tests to determine if the individual distributions tested are good
% fits to the data. Then, we apply the Vuong Test (bbvuong) to combinations
% of these distributions to determine which is best in the case that both
% are deemed good fits by the previous hypothesis test. Finally, the
% skewness and kurtosis of the data is recorded.
obs = nan(20,1); n = 20; depth = 5:10:200; ad = nan(4,20); ks = nan(5,20);
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = cOutB(pOutB==i);
    % apply statistical tests to the data   
    if length(X_i) > 3
        gammaParams = mle(X_i,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        [~,ad(2,i)] = adtest(X_i,"Distribution","logn","Alpha",0.005);
        [~,ad(1,i)] = adtest(X_i,"Distribution","norm","Alpha",0.005);
        [~,ad(3,i)] = adtest(X_i,"Distribution","weibull");
        [~,ad(4,i)] = adtest(X_i,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

% Remove values that do not meet the threshold
for i = 1:n
    if obs(i) < threshold
        ad(:,i) = nan;
        ks(:,i) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
    end
end

% Remove these nan values
tmp = [];
for i = 1:n
    if ~isnan(sum(ad(:,i)))
        tmp = [tmp i];
    end
end
p = depth(tmp);
rV = rV(:,tmp);
pV = pV(:,tmp);
ks = ks(:,~all(isnan(ks)));
ad = ad(:,~all(isnan(ad)));


% Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(p));
if testSel == 4
    for i = 1:length(p)
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
elseif testSel == 2         % normal-lognormal
    for i = 1:length(p)
        if rV(1,i) > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
        end
    end
end

% Plot the results.
if suppressFig == false
    ax = figure;
    n = length(p);
    
    subplot(1,3,3)
    barh(obs,'FaceColor','#d3d3d3');
    hold on
    xline(threshold);
    hold off
    set(gca,'YDir','reverse');
    ylim([binLimA binLimB]);
    yticks(0.5:1:19.5);
    yticklabels(0:10:200);
    yticklabels({});
    xlabel('No. of Obs.',Interpreter='latex',FontSize=13);
    
    subplot(1,3,[1 2])
    xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
    hold on
    if strcmp(hypTest,"ks")
        if testSel == 4
            plot(ks(1,:),p,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ks(2,:),p,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
            plot(ks(3,:),p,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
            plot(ks(4,:),p,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
        elseif testSel == 2
            plot(ks(1,:),p,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ks(2,:),p,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        end
        xlabel('K-S $p$-value',Interpreter='latex',FontSize=13);
    else
        if testSel == 4
            plot(ad(1,:),p,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ad(2,:),p,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
            plot(ad(3,:),p,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
            plot(ad(4,:),p,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
        elseif testSel == 2
            for i = 1:n
                if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr && ad(2,i) > alphaHy
                    plot(ad(1,i),p(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
                elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr && ad(2,i) > alphaHy
                    plot(ad(1,i),p(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
                elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr && ad(1,i) > alphaHy
                    plot(ad(2,i),p(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
                elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr && ad(1,i) > alphaHy
                    plot(ad(2,i),p(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
                end
            end
            plot(nan,nan,'square','Color','#808080','MarkerSize',15,'DisplayName','V-LLR best fit (p > 0.1)');        
            plot(nan,nan,'square','Color','#808080','MarkerSize',15,'LineWidth',4,'DisplayName','V-LLR best fit (p < 0.1)');        
            plot(ad(1,:),p,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
            plot(ad(2,:),p,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        end
        xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
    end
    hold off
    if logAxis == true
        set(gca, 'XScale', 'log');
    end
    grid minor;
    ylim(limits);
    xlim([0.1*alphaHy 1]); ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
    set(gca,'YDir','reverse');
    if season == 0
        legend('Position',[0.4 0.7 0.07 0.12],FontSize=11);
    else
        legend(Location="best",FontSize=11);
    end
    
    
else
        ax = nan;
end

end