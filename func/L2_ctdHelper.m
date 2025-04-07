function ax = L2_ctdHelper(X,pIn,maxMld,dcm,testSel,hypTest,limits,threshold,logAxis,legendOn)
%L2_ctdHelper Calculate statistical tests for CTD variables in the DCM layer (L2).
% Specifically, we calculate the hypothesis test results (A-D or K-S) here.
% Additionally, we also calculate the Vuong Log-Likelihood Ratio (V-LLR)
% test, which tells us the best distribution when A-D or K-S predicts that
% more than one distribution fits the data at a given depth.

% INPUTS
% X: the variable to be tested. 
% pIn: the pressure coordinates of the variable to be tested.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% dcm: the depth of the deep chlorophyll maximum (in dbar)
% testSel: choose between 2 (normal-lognormal) or 4
% (normal-lognormal-weibull-gamma).
% hypTest: choose "ad" to run the Anderson-Darling (A-D) test or "ks" to
% run the Kolmogorov-Smirnov (K-S) test.
% limits: y-limits (pressure coordinates) to output
% threshold: the statistical test (A-D) will only run above this threshold.
% logAxis: if set to true, outputs the x-axis (p-values) in log scale.
% season: 0 = no seasonal analysis; otherwise 1-4 => run seasonal analysis
% on winter (1), spring (2), summer (3), or autumn (4).
% suppressFig: if set to true, the figure does not output.

% OUTPUTS
% ax = figure identifier, needed for labelling and export,

% Set up default values.
if nargin < 10
    legendOn = true;
end
if nargin < 9
    logAxis = true; % true => output p-values in log x-axis, otherwise no log plot.
end
if nargin < 8
    threshold = 30;          
end                        
if nargin < 7
    limits = [-80 140];
end
if nargin < 6
    hypTest = "ad";
end
if nargin < 5
    testSel = 2;
end

% Options for the plot.
n = length(pIn);        % Depth range
nT = length(X(1,:));    % Time range
alphaHy = 0.005;        % Alpha for K-S/A-D p-value
alphaLlr = 0.1;         % Alpha for Vuong LLR p-value

% Extract data beneath the mixed layer.
pSubml = nan(n,nT);
xSubml = nan(n,nT);
for i = 1:n-1
    for j = 1:nT
        if pIn(i) >= maxMld(j)
            pSubml(i,j) = pIn(i+1);
            xSubml(i,j) = X(i+1,j);
        end
    end
end

% Convert vertical coordinates such that it centres on the DCM.
pL = nan(n,nT);
for i = 1:nT
    pL(:,i) = pSubml(:,i) - dcm(i);
end

% 3. Find KS p-value, Vuong p-value and LLR, skewness, and kurtosis for 
% each depth.

% Determine the maximum and minimum depths and the depth range.
l1 = min(min(pL));
l2 = max(max(pL));
range = l1:2:l2;
rangeLen = 1:1:length(range);
n2 = length(range);
bottom = floor((l1 - limits(1))./2);
top = floor((l2 - limits(2))./2); 

% Apply statistical analysis. Determine the A-D p-value 'ad', the K-S 
% p-value 'ks', the Vuong Log-Likelihood Ratio 'rV' (and associated
% p-value 'pV'), the number of observations per depth 'obs'.
ks = nan(5,n2);
ad = nan(4,n2);
rV = nan(10,n2); 
pV = nan(10,n2);
obs = nan(1,n2);

for i = rangeLen
    tmp = xSubml(pL==range(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        gammaParams = mle(tmp,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn");
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm");
        [~,ad(3,i)] = adtest(tmp,"Distribution","weibull");
        [~,ad(4,i)] = adtest(tmp,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
    end
end

% Set results = nan if they do not meet the threshold
for i = rangeLen
    if obs(i) < threshold
        ks(:,i) = nan;
        ad(:,i) = nan;
        pV(:,i) = nan;
        rV(:,i) = nan;
    end
end

% Vuong's Log-Likelihood Ratio (V-LLR) Test: Set up the visualisation.
vuongRes = zeros(1,n2);
annot = strings(1,n2);
anClr = strings(1,n2);
anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
tmpEmph = strings(1,n2); tmpEmph(cellfun(@isempty,tmpEmph)) = 'bold';
rV(isnan(rV)) = 0;

% NOT NEEDED? It is needed for vuongRes!!!
% V-LLR: Normal vs Lognormal vs Weibull vs Gamma
if testSel==4 && hypTest == "ks"
    for i = 1:n2
        if rV(1,i) > 0  & rV(2,i) > 0 & rV(3,i) > 0 & ks(1,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) > 1
                vuongRes(i) = 1;
            end
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0 & ks(2,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) > 1
                vuongRes(i) = 2;
            end
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0 & ks(3,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) > 1
                vuongRes(i) = 3;
            end
        end
    end
    rV(rV==0) = nan;
elseif testSel == 4 && hypTest == "ad"
    %%%
    for i = 1:n2
        if rV(1,i) > 0 & rV(2,i) > 0 & rV(3,i) > 0 & ad(1,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) > 1
                vuongRes(i) = 1;
            end
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0 & ad(2,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) > 1
                vuongRes(i) = 2;
            end
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0 & ad(3,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) > 1
                vuongRes(i) = 3;
            end
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0 & ad(4,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) > 1
                vuongRes(i) = 4;
            end
        end
    end
    rV(rV==0) = nan;
    %%%
elseif testSel == 2
    % 4.b. Vuong: Normal Vs. Lognormal Only
    for i = 1:n2
        if rV(1,i) > 0 
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
        end
    end
    rV(rV==0) = nan;
end

% Plot results.

ax = figure;

subplot(1,3,3)
barh(obs(rangeLen(1):rangeLen(end)),'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
xlabel('No. of Obs.','FontSize',13,Interpreter='latex');
set(gca,"YTick",2:5:n2,"YTickLabel",range(2):10:range(end));
ylim([rangeLen(1-bottom) rangeLen(end-top)]);
yticklabels({})

subplot(1,3,[1 2])
xline(alphaHy,'-','\alpha',LineWidth=1.5,Color="#808080",HandleVisibility="off");
hold on
if strcmp(hypTest,"ks")
    if testSel==4
        plot(ks(1,:),range,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),range,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(3,:),range,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ks(4,:),range,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel==2
        plot(ks(1,:),range,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),range,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('K-S $p$-value','Interpreter','latex',FontSize=15);
elseif strcmp(hypTest,"ad")
    if testSel==4
        plot(ad(1,:),range,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),range,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(3,:),range,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ad(4,:),range,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel==2
        for i = 1:n
            if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr && ad(2,i) > alphaHy
                plot(ad(1,i),range(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
            elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr && ad(2,i) > alphaHy
                plot(ad(1,i),range(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr && ad(1,i) > alphaHy
                plot(ad(2,i),range(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr && ad(1,i) > alphaHy
                plot(ad(2,i),range(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
            end
        end
        plot(ad(1,:),range,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),range,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('A-D $p$-value',FontSize=13,Interpreter='latex');
end
hold off
grid minor;
if logAxis == true
    set(gca, 'XScale', 'log');
end
ylim(limits); xlim([0.1*alphaHy 1]); 
set(gca,"YTick",limits(1):10:limits(2),"YTickLabel",limits(1):10:limits(2));
ylabel("Pressure [dbar]",Interpreter="latex",FontSize=13);
set(gca,'YDir','reverse');
if legendOn == true
    legend(Location="best");
end

end