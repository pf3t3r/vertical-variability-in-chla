function [ax,ks,obs,pB,X,ad] = L0_helper(tmp,threshold,hypTest,logAxis,season,figSuppress)
%L0_helper
%INPUT: tmp = text-file with pressure, bottle concentration, and bottle ID;
% threshold = threshold needed for test to be run (default for A-D = 30);
% hypTest = choose hypothesis test ("ad" for A-D or "ks" for K-S);
% logAxis = output x-axis of A-D test in log coordinates (true/false);
% season = do seasonal analysis or not. 0 = no seasonal analysis, while 1-4
% implies seasonal analysis (1 = winter, 2 = spring, 3 = summer, 4 =
% winter);
% figSuppress = suppress output of figure. 
%OUTPUT: ax = needed to save figure, ks = K-S p-value, obs = no. of 
% observations, pB = binned pressure, X = output (cleaned) bottle
% concentration, ad = A-D p-value. Any of these values can be left empty if
% you just want to show the plots.

meanDcm = load("data/hot_meanDcm.mat").timeSeriesMeanDcm;
meanPrc = load("data/hot_percDcm.mat").timeSeriesPrctl;

if nargin <6
    figSuppress = false;
end

if nargin < 5
    season = 0;
    % 0 = no seasonal analysis, 1 = winter, 2 = spring, 3 = summer,
    % 4 = autumn
end

if nargin < 4
    logAxis = true;
    % true => output p-values in log x-axis, otherwise no log plot.
end

if nargin < 3
    hypTest = 'ad';
end

if nargin <2
    threshold = 30;
end
tmpT = ""; alphaHy = 0.005;

if isstruct(tmp)
    if length(tmp.data(1,:))==6
        botId = tmp.data(:,2);
        pIn = tmp.data(:,5);
        X = tmp.data(:,6);
    else 
        botId = tmp.data(:,2);
        pIn = tmp.data(:,4);
        X = tmp.data(:,5);
    end
else
    botId = tmp(:,2);
    pIn = tmp(:,4);
    X = tmp(:,5);
end
n = length(pIn);
% nB = length(botId);
% n3 = length(X);

if season ~= 0
    botId(botId==-9) = nan;
    botId2 = num2str(botId);
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
        X = X(winIds);
        pIn = pIn(winIds);
        nSea = length(winIds);
    elseif season == 2
        X = X(sprIds);
        pIn = pIn(sprIds);
        nSea = length(sprIds);
    elseif season == 3
        X = X(sumIds);
        pIn = pIn(sumIds);
        nSea = length(sumIds);
    elseif season == 4
        X = X(autIds);
        pIn = pIn(autIds);
        nSea = length(autIds);
    end
end

%% Bin
pB = discretize(pIn,0:10:200);
pX = nan(n,1);

if season == 0
    for i = 1:n
        pX(i) = pB(i)*10 - 5;
    end
else
    for i = 1:nSea
        pX(i) = pB(i)*10 - 5;
    end
end

%% K-S: bottle chla

n2 = 20;
ks = nan(5,n2);
obs = nan(1,n2);
ad = nan(2,n2);

% Calculate KS p-value, skewness, kurtosis
for i = 1:n2
    % find concentration X_i at binned pressure i
    X_i = X(pB==i);
    % remove negative or null values
    X_i(X_i <= 0) = [];
    % apply KS test to X_i (only when at least 3 values at binned pressure)
    if length(X_i) > 3
        if strcmp(hypTest,'ks')
            [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        else
            [~,ad(1,i)] = adtest(X_i,'Distribution','logn','alpha',0.005);
            [~,ad(2,i)] = adtest(X_i,'Distribution','norm','alpha',0.005);
        end
    end
    obs(i) = length(X_i);
    clear X_i;
end

ks(:,obs<threshold) = nan;
ad(:,obs<threshold) = nan;

pXX = 5:10:195;

if figSuppress == false
    ax = figure;
    subplot(1,3,3)
    barh(obs,'FaceColor','#d3d3d3');
    xline(threshold);
    ylim([0.5 20.5]);
    yticks(0.5:2:20.5);
    yticklabels({});
    set(gca,"YDir","reverse");
    xlabel("No. of Obs.",Interpreter="latex",FontSize=13);
    
    subplot(1,3,[1 2])
    xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13); 
    hold on
    if strcmp(hypTest,'ks')
        plot(ks(2,:),pXX,'o-','Color','#4d9221',LineWidth=1.5,MarkerSize=5,DisplayName="Lognormal");
        xlabel('K-S $p$-value',Interpreter='latex',FontSize=13);
    else
        plot(ad(2,:),pXX,'o-','Color','#c51b7d',LineWidth=1.5,MarkerSize=5,DisplayName="Normal");
        hold on
        plot(ad(1,:),pXX,'o-','Color','#4d9221',LineWidth=1.5,MarkerSize=5,DisplayName="Lognormal");
        yline(meanDcm,DisplayName="$\bar{p}_{dcm} \pm$ 5/95");
        yline(meanPrc(1),LineWidth=0.25,HandleVisibility="off");
        yline(meanPrc(2),LineWidth=0.25,HandleVisibility="off");
        hold off
        xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
    end
    if logAxis == true
        set(gca, 'XScale', 'log');
    end
    ylim([0 200]); xlim([0.1*alphaHy 1]);
    set(gca,'YDir','reverse');
    ylabel("Pressure [dbar]",Interpreter="latex",FontSize=13);
    grid minor;
    legend(FontSize=13,Interpreter="latex");
else
    ax = nan;
end

end