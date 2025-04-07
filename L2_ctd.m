% Script to output L2 ctd results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\"); addpath("output\");
set(groot,'defaultFigureUnits','centimeters','defaultFigurePosition',[3 3 15 15]);

% Load MLD and DCM
maxMld = load('mldVals.mat').maxMld;
dcm = load("dcm.mat").meanPcm;  

% Assign lower bound on pressure
lowerP = 129;
pIn = 0:2:2*(lowerP-1);

% Options and Test Cases.
thresh = 30;
principalAnalysisAd = false;
principalAnalysisKs = false;
setUpSeasonalAnalysis = true;
seasonalAnalysisAd = true;
seasonalAnalysisKs = false;
dists = 2;                      % 2 = norm + logn;
                                % 4 = norm + logn + weib + gamm
logAxes = true;                 % output p-values as log values (true)
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Load hydrographical variables

ctdData = load("data\hot_ctd.mat").iso;
timeData = load("data\hot_ctd.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);
clear msng cR;

% Temperature
Tin = nan(329,lowerP,31);
T = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).t(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            Tin(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    T(:,i) = mean(squeeze(Tin(i,:,:)),2,"omitnan");
end

% Salinity
SPin = nan(329,lowerP,31);
Sp = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).sp(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            SPin(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    Sp(:,i) = mean(squeeze(SPin(i,:,:)),2,"omitnan");
end

O2in = nan(329,lowerP,31);
o2 = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).o(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            O2in(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    o2(:,i) = mean(squeeze(O2in(i,:,:)),2,"omitnan");
end

clear i j cRm O2in SPin Tin tmp;

%% Additional Processing.
% Convert temperature and practical salinity into conservative temperature
% and absolute salinity, and then calculate potential density.

stnALOHA_lon = -158;
stnALOHA_lat = 22.75;
SA = gsw_SA_from_SP(Sp,pIn',stnALOHA_lon,stnALOHA_lat);
CT = gsw_CT_from_t(SA,T,pIn');
sigma0 = gsw_sigma0(SA,CT);
[N2,p_mid] = gsw_Nsquared(SA,CT,pIn',stnALOHA_lat);

%% Principal Analysis
if principalAnalysisAd == true
    
    tmpX = ": L2";

    % Chlorophyll a (88-21)
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,131:329);
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('chl-$a$'+tmpX,"Interpreter","latex");

    % Conservative Temperature (88-21)
    L2_ctdHelper(CT,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh,true,true);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");
    
    % Absolute Salinity (88-21)
    L2_ctdHelper(SA,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh,true,true);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % Potential Density (88-21)
    L2_ctdHelper(sigma0,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh,true,true);
    sgtitle('$\sigma_0$'+tmpX,"Interpreter","latex");
    
    % Stratification (N^2)
    L2_ctdHelper(N2,p_mid(:,1),maxMld,dcm,dists,"ad",[-60 60],thresh,true,true);
    sgtitle('$N^2$'+tmpX,"Interpreter","latex");
    
    % O2 (88-21)
    L2_ctdHelper(o2,pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");
    
end

%% Set up seasonal analysis

if setUpSeasonalAnalysis == true

    % Set-up seasonal test.

    % find "average" month of a cruise, then give an ID to that cruise 
    % saying which season it is (Spring, Summer, Autumn, or Winter).
    
    avgMonth = nan(329,1);
    winter = nan(329,1);
    spring = nan(329,1);
    summer = nan(329,1);
    autumn = nan(329,1);
    
    for i = 1:329
        tmp = round(mean(month(datetime(timeData(i).date,"ConvertFrom","datenum"))));
        
        if (tmp == 12) || (tmp == 1) || (tmp == 2)
            winter(i) = 1;
        end
        if (tmp == 3) || (tmp == 4) || (tmp == 5)
            spring(i) = 1;
        end
        if (tmp == 6) || (tmp == 7) || (tmp == 8)
            summer(i) = 1;
        end
        if (tmp == 9) || (tmp == 10) || (tmp == 11)
            autumn(i) = 1;
        end
    
        avgMonth(i) = tmp;
    end
    
    % ASTRO seasons -> winter = jfm; spring = amj; summer = jas; autumn = ond.
    
    winIds = []; sprIds = []; sumIds = []; autIds = []; 
    for i = 1:329
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
end
  

%% Seasonal Analysis.
if seasonalAnalysisAd == true
    % A-D
    % winter
    %tmpT = "-ad-01";
    tmpX = ": L2 Winter";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,winIds(36:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L2_ctdHelper(CT(:,winIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % SA
    L2_ctdHelper(SA(:,winIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,winIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

    % spring
    tmpT = "-ad-02"; tmpX = ": L2 Spring";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,sprIds(33:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L2_ctdHelper(CT(:,sprIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % Sp
    L2_ctdHelper(SA(:,sprIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,sprIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

    % SUMMER
    tmpX = ": L2 Summer";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,sumIds(33:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % T
    L2_ctdHelper(CT(:,sumIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % Sp
    L2_ctdHelper(SA(:,sumIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,sumIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

    % AUTUMN
    tmpX = ": L2 Autumn";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,autIds(30:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ad",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L2_ctdHelper(CT(:,autIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % Sp
    L2_ctdHelper(SA(:,autIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,autIds),pIn,maxMld,dcm,dists,"ad",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

end

%% Extra Analysis with Kolmogorov-Smirnov (K-S).
if principalAnalysisKs == true

    % K-S
    tmpT = ""; tmpX = ": L2";
    
    % Chlorophyll a (88-21)
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,131:329);
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('chl-$a$'+tmpX,"Interpreter","latex");
    
    % Conservative Temperature (88-21)
    L2_ctdHelper(CT,pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('\Theta'+tmpX);
       
    % Absolute Salinity (88-21)
    L2_ctdHelper(SA,pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('S_A'+tmpX);
        
    % O2 (88-21)
    L2_ctdHelper(o2,pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('O_2'+tmpX);
   
end

% Seasonal Analysis with K-S.
if seasonalAnalysisKs == true
    
    % winter
    tmpX = ": L2 Winter";
    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,winIds(36:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ks",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % T
    L2_ctdHelper(CT(:,winIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % Sp
    L2_ctdHelper(SA(:,winIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,winIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

    % spring
    tmpX = ": L2 Spring";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,sprIds(33:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ks",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % T
    L2_ctdHelper(CT(:,sprIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % Sp
    L2_ctdHelper(SA(:,sprIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,sprIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

    % summer
    tmpX = ": L2 Summer";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,sumIds(33:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ks",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % T
    L2_ctdHelper(CT(:,sumIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % Sp
    L2_ctdHelper(SA(:,sumIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,sumIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");

    % autumn
    tmpX = ": L2 Autumn";

    % chla
    chla = load("data\hot_fluo_L1_L2").meanLiN(1:lowerP,autIds(30:end));
    L2_ctdHelper(chla,pIn,maxMld,dcm,dists,"ks",[-60 60],thresh);
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L2_ctdHelper(CT(:,autIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$\Theta$'+tmpX,"Interpreter","latex");

    % SA
    L2_ctdHelper(SA(:,autIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$S_A$'+tmpX,"Interpreter","latex");

    % O2
    L2_ctdHelper(o2(:,autIds),pIn,maxMld,dcm,dists,"ks",[-60 60]);
    sgtitle('$O_2$'+tmpX,"Interpreter","latex");
end