% Statistical Analysis of the mixed layer (L1) at Station ALOHA for CTD
% data.
% We import the data and run a hypothesis test on it with the respective
% null hypotheses of normal and lognormal. We primarily use the 
% Anderson-Darling test since this is both more powerful than similar tests
% such as Kolmogorov-Smirnov and more flexible than tests such as 
% Shapiro-Wilks which did not easily allow for testing of other 
% distributions.

close all; clc; clear;
addpath("func\"); addpath("output\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Options and Test Cases.
principalAnalysisAd = true;     % main analysis with A-D test
seasonalAnalysisAd = false;     % seasonal analysis with A-D test
principleAnalysisKs = false;    % as above, but with K-S test
seasonalAnalysisKs = false;       % seasonality of statistics
noOfDists = 2;                  % 2 = norm + logn;
                                % 4 = norm + logn + weib + gamm
threshAd = 30;                  % only for chl-a in A-D
threshKs = 50;
logAxes = true;                 % output p-values as log values (true)

%% Load Data: processed fluorometry and mixed layer depth.

% The fluorometry data imported here has been previously saved as 
% cruise-means using only night-time casts to avoid the effect of
% non-photochemical quenching.

% Chloropigment fluorescence
% It was previously determined that fluorometry results from the initial
% years of Station ALOHA were different to those of later years. Therefore,
% we only use the results from the second, newer fluorometer.
F2 = [131 329];   
epN = load("data\hot_fluo_L1_L2.mat").meanEpN;
epN = epN(1:76,F2(1):F2(2));
pIn = 0:2:150;
maxMld = load("output/mldVals.mat").maxMld;

ctdData = load("data\hot_ctd.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);

% Temperature T
T = nan(329,76,31);
meanT = nan(76,329);
for i = cRm
    tmp = ctdData(i).t(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            T(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanT(:,i) = mean(squeeze(T(i,:,:)),2,"omitnan");
end

% Salinity SP
% S missing data same as for T? Yes!
SP = nan(329,76,31);
meanSp = nan(76,329);
for i = cRm
    tmp = ctdData(i).sp(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            SP(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanSp(:,i) = mean(squeeze(SP(i,:,:)),2,"omitnan");
end

% O2
O2 = nan(329,76,31);
meanO2 = nan(76,329);
for i = cRm
    tmp = ctdData(i).o(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            O2(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanO2(:,i) = mean(squeeze(O2(i,:,:)),2,"omitnan");
end

%% Further Processing.
% Convert temperature and practical salinity into conservative temperature
% and absolute salinity, and then calculate potential density sigma0 and
% stratification N2.

stnALOHA_lon = -158;
stnALOHA_lat = 22.75;
SA = gsw_SA_from_SP(meanSp,pIn',stnALOHA_lon,stnALOHA_lat);
CT = gsw_CT_from_t(SA,meanT,pIn');
sigma0 = gsw_sigma0(SA,CT);
[N2,p_mid] = gsw_Nsquared(SA,CT,pIn',stnALOHA_lat);

%% Visualise CT, SA, and sigma0.

% Conservative Temperature CT
figure;
scatter(CT,pIn,4,[0.6 0.6 0.6]); title("Conservative Temperature",Interpreter="latex");
xlabel("$\Theta$ [$^\circ$C]",Interpreter="latex"); ylabel("Pressure [dbar]");
set(gca,"YDir","reverse"); grid on;

% Absolute Salinity SA
figure;
scatter(SA,pIn,4,[0.6 0.6 0.6]); title("Absolute Salinity",Interpreter="latex");
xlabel("S_A [g/kg]"); ylabel("Pressure [dbar]");
set(gca,"YDir","reverse"); grid on;

% Potential Density Anomaly (wrt p=0) 'sigma0'
figure;
scatter(sigma0,pIn,5,[0.6 0.6 0.6]); title("Potential Density Anomaly",Interpreter="latex")
xlabel("$\sigma_0$ [kg/m$^3$]",Interpreter="latex"); ylabel("Pressure [dbar]");
set(gca,"YDir","reverse"); grid on;

% Stratification 'N^2'
figure;
scatter(N2,p_mid,5,[0.6 0.6 0.6]); title("Stratification",Interpreter="latex")
xlabel("$N^2$ [cycles $s^{-1}$]",Interpreter="latex"); ylabel("Pressure [dbar]");
set(gca,"YDir","reverse"); set(gca,'xscale','log'); grid on;

%% Principal Analysis: A-D.
if principalAnalysisAd == true

    tmpT = "";
    
    % CHL-A
    L1_ctdHelper(epN,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("chloropigment fluorescence");
    % save("output\L1_ctd\hot_fluo_L1_L2.mat","p","ad","ks","obs");

    % CT
    L1_ctdHelper(CT,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("\Theta");
    
    % SA
    L1_ctdHelper(SA,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("CTD S_A 88-21: L1"+tmpT);
   
    % sigma0
    L1_ctdHelper(sigma0,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("CTD \sigma_0 88-21: L1"+tmpT);
   
    % N^2 Stratification
    L1_ctdHelper(N2,p_mid(:,1),maxMld,threshAd,noOfDists,"ad");
    sgtitle("N^2 88-21: L1"+tmpT);
        
    % O2
    L1_ctdHelper(meanO2,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("CTD O2 88-21: L1"+tmpT);
end

%% Seasonal Analysis: A-D.

% find "average" month of a cruise, then give an ID to that cruise saying
% which season it is (Spring, Summer, Autumn, or Winter).

avgMonth = nan(329,1);
winter = nan(329,1);
spring = nan(329,1);
summer = nan(329,1);
autumn = nan(329,1);

for i = 1:329
    tmp = round(mean(month(datetime(ctdData(i).date,"ConvertFrom","datenum"))));
    
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

% Calculate average mld per season
mldW = mean(maxMld(winIds),"omitnan");
mldSp = mean(maxMld(sprIds),"omitnan");
mldSu = mean(maxMld(sumIds),"omitnan");
mldA = mean(maxMld(autIds),"omitnan");

% Average dcm per season
dcm = load("dcm.mat").meanPcm;  
dcmW = mean(dcm(winIds),"omitnan");
dcmSp = mean(dcm(sprIds),"omitnan");
dcmSu = mean(dcm(sumIds),"omitnan");
dcmA = mean(dcm(autIds),"omitnan");

% Naming scheme
% -01 = winter, -02 = spring, -03 = summer, -04 = autumn

if seasonalAnalysisAd == true

    % WINTER
    tmpX = ": L1 Winter";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,winIds(36:end));
    L1_ctdHelper(chla,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("chl-a" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,winIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$\Theta$ " + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,winIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$S_A$ " + tmpX,"Interpreter","latex");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,winIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$O_2$ " + tmpX,"Interpreter","latex");

    % SPRING
    tmpX = ": L1 Spring";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,sprIds(33:end));
    L1_ctdHelper(chla,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,sprIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,sprIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,sprIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$O_2$" + tmpX,"Interpreter","latex");

    % SUMMER
    tmpX = ": L1 Summer";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,sumIds(33:end));
    L1_ctdHelper(chla,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,sumIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,sumIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,sumIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$O_2$" + tmpX,"Interpreter","latex");

    % AUTUMN
    tmpX = ": L1 Autumn";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,sumIds(30:end));
    L1_ctdHelper(chla,pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,autIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,autIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,autIds),pIn,maxMld,threshAd,noOfDists,"ad");
    sgtitle("$O_2$ " + tmpX,"Interpreter","latex");
end

%% Extra Analysis with Kolmogorov-Smirnov (K-S).
if principleAnalysisKs == true

    tmpX = ": L1";
    
    % CHL-A
    L1_ctdHelper(epN,pIn,maxMld,threshKs,noOfDists);
    sgtitle("Fluorescence"+tmpX,"Interpreter","latex");
    
    % Ct
    L1_ctdHelper(CT,pIn,maxMld,threshKs,noOfDists);
    sgtitle("$\Theta$"+tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA,pIn,maxMld,threshKs,noOfDists);
    sgtitle("S_A"+tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2,pIn,maxMld,threshKs,noOfDists);
    sgtitle("O_2"+tmpX,"Interpreter","latex");
    
end

% Seasonal Analysis: K-S.
if seasonalAnalysisKs == true

    % WINTER
    tmpX = ": L1 Winter";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,winIds(36:end));
    L1_ctdHelper(chla,pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,winIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,winIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,winIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpX,"Interpreter","latex");

    % SPRING
    %tmpT = "-ks-02";
    tmpX = ": L1 Spring";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,sprIds(33:end));
    L1_ctdHelper(chla,pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,sprIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,sprIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,sprIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpX,"Interpreter","latex");

    % SUMMER
    tmpX = ": L1 Summer";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,sumIds(33:end));
    L1_ctdHelper(chla,pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,sumIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,sumIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,sumIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpX,"Interpreter","latex");

    % AUTUMN
    tmpX = ": L1 Autumn";

    chla = load("data\hot_fluo_L1_L2.mat").meanEpN(1:101,sumIds(30:end));
    L1_ctdHelper(chla,pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("chl-$a$" + tmpX,"Interpreter","latex");

    % CT
    L1_ctdHelper(CT(:,autIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$\Theta$" + tmpX,"Interpreter","latex");
    
    % SA
    L1_ctdHelper(SA(:,autIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$S_A$" + tmpX,"Interpreter","latex");
    
    % O2
    L1_ctdHelper(meanO2(:,autIds),pIn,maxMld,threshKs,noOfDists,"ks");
    sgtitle("$O_2$" + tmpX,"Interpreter","latex");
end
