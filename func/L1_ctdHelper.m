function ax = L1_ctdHelper(X,pIn,maxMld,threshold,testSel,hypTest,logAxis)
%L1_ctdHelper Calculate statistical tests for CTD variables in the mixed layer (L1).
% Specifically, we calculate the hypothesis test results (A-D or K-S) here.
% Additionally, we also calculate the Vuong Log-Likelihood Ratio (V-LLR)
% test, which tells us the best distribution when A-D or K-S predicts that
% more than one distribution fits the data at a given depth.

% INPUTS
% X: the variable to be tested. 
% pIn: the pressure coordinates of the variable to be tested.
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
% ax = figure identifier, needed for labelling and export.

if nargin < 7
    logAxis = true;
end
if nargin < 6
    hypTest = "ad";
end
if nargin < 5
    testSel = 2;
end
if nargin < 4
    threshold = 30;
end

% Options for the plot.
limits = [0 100];
alphaHy = 0.005;
alphaLlr = 0.10;

n = length(pIn);

% Determine how much of the variable is within the mixed layer.
mldCon = nan(size(X));
for k = 1:length(X(1,:))
    for j = 1:n
        if pIn(j) < maxMld(k)
            mldCon(j,k) = X(j,k);
        end
    end
end

% Apply statistical analysis. Determine the A-D p-value 'ad', the K-S 
% p-value 'ks', the Vuong Log-Likelihood Ratio 'rV' (and associated
% p-value 'pV'), the number of observations per depth 'obs'.
ad = nan(4,n); ks = nan(5,n); rV = nan(10,n); pV = nan(10,n);
obs = nan(1,n);

for i = 1:n
    tmp = mldCon(i,:);
    tmp(isnan(tmp) | tmp<0) = 0;
    tmp(tmp==0) = [];
    if length(tmp) > 3
        gammaParams = mle(tmp,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn");
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm");
        [~,ad(3,i)] = adtest(tmp,"Distribution","weibull");
        [~,ad(4,i)] = adtest(tmp,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        obs(i) = length(tmp);
    end
end

% Set values = NaN if they do not meet the threshold (i.e. there are
% insufficient values at that depth).
for i = 1:n
    if obs(i) < threshold
        ks(:,i) = nan;  
        ad(:,i) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
    end
end

% Remove those NaN values.
tmp = [];
for i = 1:n
    if ~isnan(sum(ad(:,i)))
        tmp = [tmp i];
    end
end
tr = pIn(tmp);
rV = rV(:,tmp);
pV = pV(:,tmp);
ad = ad(:,~all(isnan(ad)));
ks = ks(:,~all(isnan(ks)));

% Vuong's Log-Likelihood Ratio (V-LLR) Test: Set up the visualisation.
vuongRes = nan(1,length(tr));
if testSel == 4
    for i = 1:length(tr)
        if rV(1,i) > 0 & rV(2,i) > 0 & rV(3,i) > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
            vuongRes(i) = 2;
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
            vuongRes(i) = 3;
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
            vuongRes(i) = 4;
        end
    end
elseif testSel == 2
    % Normal vs Lognormal Case.
    for i = 1:length(tr)
        if rV(1,i) > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
        end
    end
end

n2 = length(tr);

% Plot the results.

ax = figure;

subplot(1,3,3)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
ylim([limits(1)+1 limits(2)/2 + 1]);
yticks(1:5:55);
yticklabels(0:10:150);
yticklabels({});
xlabel('No. of Obs.',Interpreter='latex',FontSize=13);

subplot(1,3,[1 2])
xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
hold on
if strcmp(hypTest,"ks")
    if testSel == 4
        plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),tr,'+-','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ks(4,:),tr,'.-','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel == 2
        plot(ks(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),tr,'+--','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('K-S $p$-value',Interpreter='latex',FontSize=13);
else
    if testSel == 4
        plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),tr,'+-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(3,:),tr,'x-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ad(4,:),tr,'.-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel == 2
        for i = 1:n2
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
        plot(nan,nan,'square','Color','#808080','MarkerSize',15,'DisplayName','V-LLR best fit (p > 0.1)');        
        plot(nan,nan,'square','Color','#808080','MarkerSize',15,'LineWidth',4,'DisplayName','V-LLR best fit (p < 0.1)');        
        plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
end
hold off
if logAxis == true
    set(gca, 'XScale', 'log');
end
grid minor;
ylim(limits); xlim([0.1*alphaHy 1]); ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
set(gca,'YDir','reverse');
legend(Location="best",FontSize=11);
end