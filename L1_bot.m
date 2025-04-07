% Statistical Analysis of the mixed layer (L1) at Station ALOHA for chl-a, other pigments, and BGC variables.
% We import the data and run a hypothesis test on it with the respective
% null hypotheses of normal and lognormal. We use the Anderson-Darling test
% since this is both more powerful than similar tests such as
% Kolmogorov-Smirnov and more flexible than tests such as Shapiro-Wilks
% which did not easily allow for testing of other distributions.

clear; clc; close all;
addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Options and Test Cases.
thresh = 30;                    % Threshold for A-D Test to be accepted
principle = true;               % Main analysis: A-D Test complete time series
seasonal = false;               % Seasonal analysis: A-D test applied to 
                                % four seasons
showHistograms = false;         % Show histogram analysis (true/false).
testSel = 2;                    % 2 = norm + logn (default), 
                                % 4 = norm + logn + weib + gamm
logAxes = true;                 % Output p-values as log values on x-axis
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Calculate Mixed Layer Depth
% Specifically, find the Maximum Mixed Layer Depth recorded for each 
% cruise "maxMld".

ctdData = importdata("data\hot_ctd.mat").ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;
save output/mldVals.mat maxMld;

%% Check histograms

if showHistograms == true

    % Select dataset and bin for analysis
    D = "data/hot_chla.txt";
    bin = 10;   % value here represents the midpoint of a bin,
                % e.g. 10 = bin from 5-15 dbar, 20 = bin from 15-25 dbar
    nameVar = "chl-$a$";

    % Import data; extract data within mixed layer.
    tmp = importdata(D);
    [~,~,~,~,~,~,~,X_ML,p_ML] = L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes,0,true);
    
    % Bin pressure
    pB = round(p_ML,-1);   
    
    % Plot 1.
    figure;
    histogram(X_ML(pB==bin)); title(""+nameVar+": "+bin+" dbar",Interpreter="latex");
    grid on

    % Plot 2. Histfit.
    figure
    histfit(X_ML(pB==bin),10,"lognormal"); title(""+nameVar+": "+bin+" dbar",Interpreter="latex");
    grid on

end

%% Principal Analysis: A-D

if principle == true

    tmpX = ": L1";
   
    % HPLC chl-a
    tmp = importdata("data/hot_chla.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    
    % chl-b
    tmp = importdata("data/hot_chlb.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    
    % Particulate Carbon
    tmp = importdata("data/hot_pc.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
end

%% Seasonal Analysis: A-D

if seasonal == true

    set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

    %%% WINTER
    tmpx = ": L1 Winter";

    % chla
    tmp = importdata("data/hot_chla.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("chl-$a$" + tmpx,"Interpreter","latex");
    clear tmp;
    
    % chlb (88-21)
    tmp = importdata("data/hot_chlb.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("chl-$b$"+tmpx,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data/hot_pc.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("Particulate Carbon"+tmpx,"Interpreter","latex");
    clear tmp;
    
    %%% SPRING
    tmpx = ": L1 Spring";

    % chla
    tmp = importdata("data/hot_chla.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("chl-$a$" + tmpx,"Interpreter","latex");
    clear tmp;
    
    % chlb (88-21)
    tmp = importdata("data/hot_chlb.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("chl-$b$"+tmpx,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data/hot_pc.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("Particulate Carbon"+tmpx,"Interpreter","latex");
    clear tmp;

    %%% SUMMER
    tmpx = ": L1 Summer";

    % chla
    tmp = importdata("data/hot_chla.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("chl-$a$" + tmpx,"Interpreter","latex");
    clear tmp;

    % chlb (88-21)
    tmp = importdata("data/hot_chlb.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("chl-$b$"+tmpx,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data/hot_pc.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("Particulate Carbon"+tmpx,"Interpreter","latex");
    clear tmp;
    
    %%% AUTUMN
    tmpx = ": L1 Autumn"; tmpT = "-ad-04";

    % chla
    tmp = importdata("data/hot_chla.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("chl-$a$" + tmpx,"Interpreter","latex");
    clear tmp;
    
    % chlb (88-21)
    tmp = importdata("data/hot_chlb.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("chl-$b$"+tmpx,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data/hot_pc.txt");
    L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("Particulate Carbon"+tmpx,"Interpreter","latex");
    clear tmp;
end