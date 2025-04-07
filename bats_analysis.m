% Script to rerun the statistical analysis developed for Station ALOHA on BATS chl-a data.
% First, we visualise the depth- and time-series and mean profile of chl-a.
% Then, we apply the statistical analysis to the three "zones":
% i. L0: the upper 200 dbar,
% ii. L1: the mixed layer,
% iii. L2: the DCM layer.
% For each zone, we check the histograms, then apply the Anderson-Darling
% (A-D) test to determine the best-fitting distribution at each depth. In
% case two distributions fit the data well, we use an additional test
% (Vuong Log-Likelihood Ratio test, V-LLR) to determine which of the two is
% best.

clear;clc;close all; format longG;

% Options and Test Cases.
threshold = 30;
showTimeSeries = true;
showL0 = true;
showL1 = true;
showL2 = true;
showOtherTests = false;     % Evaluate data with Lilliefors and chi^2 tests.
                            % (in addition to Anderson-Darling)
showL0title = true;        % Leave as 'false' for paper figures.
showDcmMldComp = false;     % show depth comparison of DCM and MLD
%% Import DCM file for BATS
% This was calculated manually on Excel for BATS cruises 1-405.
A = importdata("data\bats_dcm.txt");
dcmBats = A.data;

%% Depth- and time-series. (Fig 1a)

% Import Data: ID, Time, Depth, Chl-a.
D = importdata('data\bats_pigments.txt').data;
id = D(:,1);                % Bottle ID with format !####$$$@@, where
                            % ! = cruise type (1 = core cruise),
                            % #### = cruise number, $$$ = cast number,
                            % @@ = Niskin number.
YMD = D(:,2);               % Year, Month, and Day.
hhmm = D(:,4);              % hours and minutes.
QF = D(:,7);                % Quality control factor.
depth = D(:,8);             % Depth (m).
chla = D(:,22);             % Chl-a concentration (ng/l).
chla_Turner = D(:,24);      % Chl-a concentration, Turner method (ng/l).

% Re-assign NaNs.
chla(chla==-999) = nan;
chla_Turner(chla_Turner==-999) = nan;

% Set up time vector.
YY = nan(length(YMD),1); MM = nan(length(YMD),1); DD = nan(length(YMD),1);
hh = nan(length(YMD),1); mm = nan(length(YMD),1);
for i = 1:length(YMD)
    tmp1 = num2str(YMD(i));
    YY(i) = str2num(tmp1(1:4));
    MM(i) = str2num(tmp1(5:6));
    DD(i) = str2num(tmp1(7:8));
    clear tmp1;

    tmpB = num2str(hhmm(i));
    if length(tmpB) == 3
        hh(i) = str2num(tmpB(1));
        mm(i) = str2num(tmpB(2:3));
    elseif length(tmpB) == 4
        hh(i) = str2num(tmpB(1:2));
        mm(i) = str2num(tmpB(3:4));
    else
        hh(i) = nan;
        mm(i) = nan;
    end
end
ss = zeros(length(YMD),1);
t = datetime(YY,MM,DD,hh,mm,ss);
clear YY MM DD hh mm ss;

% Bin by depth.
edges = 0:10:200;
depth_B = discretize(depth,edges);
clear edges;

% Extract cruise type 'cType' and cruise number 'CRN'.
tmp2 = num2str(id);
cType = str2num(tmp2(:,1));
CRN = str2num(tmp2(:,2:5));
cast = str2num(tmp2(:,6:8));
clear tmp2;

% Extract core cruises. These are labelled as cType = 1.
ids = find(cType==1);
CoreCRN = CRN(ids);
cast = cast(ids);
copyOfCoreCrn = CoreCRN;
t = t(ids);
chla = chla(ids);
depth = depth(ids);
depth_B = depth_B(ids);

% Save ID where new cruise starts.
id_nc = [1];
CRN_no = [1];
for i = 2:length(CoreCRN)
    if CoreCRN(i) > CoreCRN(i-1)
        disp(i)
        id_nc = [id_nc i];
        CRN_no = [CRN_no CoreCRN(i)];
    end
end

crnAndCast = [CoreCRN cast];

% Divide up the depth, chlorophyll, and time arrays by cruise.
% NOTE that there are 405 cruises but only 398 have chl-a data.
dgrid = nan(398,35);

% For cruise #1.
dgrid(1,1:(id_nc(2)-id_nc(1))) = depth_B(id_nc(1):id_nc(2)-1);
depthUnbinned(1,1:(id_nc(2)-id_nc(1))) = depth(id_nc(1):id_nc(2)-1);
tgrid(1,1:(id_nc(2)-id_nc(1))) = t(id_nc(1):id_nc(2)-1);
chlagrid(1,1:(id_nc(2)-id_nc(1))) = chla(id_nc(1):id_nc(2)-1);

% Loop for cruises #2-397.
for i = 2:397
    dgrid(i,1:(id_nc(i+1)-id_nc(i))) = depth_B(id_nc(i):id_nc(i+1)-1);
    depthUnbinned(i,1:(id_nc(i+1)-id_nc(i))) = depth(id_nc(i):id_nc(i+1)-1);
    tgrid(i,1:(id_nc(i+1)-id_nc(i))) = t(id_nc(i):id_nc(i+1)-1);
    chlagrid(i,1:(id_nc(i+1)-id_nc(i))) = chla(id_nc(i):id_nc(i+1)-1);
end

% Final cruise (#398).
dgrid(398,1:(length(depth_B)+1-id_nc(398))) = depth_B(id_nc(398):length(depth_B));
depthUnbinned(398,1:(length(depth)+1-id_nc(398))) = depth(id_nc(398):length(depth));
tgrid(398,1:(length(t)+1-id_nc(398))) = t(id_nc(398):length(t));
chlagrid(398,1:(length(chla)+1-id_nc(398))) = chla(id_nc(398):length(chla));

tgridDatenum = datenum(tgrid);      % Convert time vector to datenum format (for plotting).

if showTimeSeries == true
    % Figure 1a. chl-a(p,t).
    figure;
    contourf(tgridDatenum,dgrid,chlagrid,linspace(0,500,10),'LineColor','auto');
    set(gca,"YDir","reverse");
    datetickzoom('x','yyyy','keeplimits');
    colormap(flipud(cbrewer2("GnBu")));
    c = colorbar;
    c.Label.String = 'chl-a [ng/l]';
    c.FontSize = 13;
    ylim([1 18]);
    zlim([0 500]);
    yticks(1:1:18);
    yticklabels(5:10:175);
    ylabel("P [dbar]","FontSize",13); xlabel("Time",FontSize=13);
    ax = gca;
    ax.FontSize = 15;
    
    %% Mean profile of time-series. (Fig. 1b)
    
    % Calculate the mean chl-a at each depth bin as well as the 5th and 95th
    % percentile values.
    chlaProfile = nan(1,20); f5 = nan(1,20); f95 = nan(1,20);
    for i = 1:20
        chlaProfile(i) = mean(chlagrid(dgrid==i),"omitnan");
        f5(i) = prctile(chlagrid(dgrid==i),5);
        f95(i) = prctile(chlagrid(dgrid==i),95);
    end
    
    % Toggle show y-label and title (for the paper we don't need either)
    displayYLabelAndTitle = false;
    
    % Figure 1b. Mean chl-a profile.
    ax = figure;
    plot(chlagrid(:,1),dgrid(:,1),'.',Color=[0.8 0.8 0.8],DisplayName="raw data");
    hold on
    plot(chlagrid(:,2:20),dgrid(:,2:20),'.',Color=[0.8 0.8 0.8],HandleVisibility='off');
    plot(chlaProfile,1:1:20,'-',"Color",[0 0 0],DisplayName="mean");
    plot(f5,1:1:20,'-',"Color",[0.5 0.5 0.5],DisplayName="5%");
    plot(f95,1:1:20,'-',"Color",[0.5 0.5 0.5],DisplayName="95%");
    hold off
    set(gca,"YDir","reverse");
    legend();
    xlabel("chl-$a$ [ng/L]",Interpreter="latex");
    if displayYLabelAndTitle == true
        title("L0 Chl-$a$ 1988-2022",Interpreter="latex");
        ylabel("P [dbar]",Interpreter="latex");
        yticklabels({});
    end
    yticks(1:1:18); yticklabels(5:10:175);
    ylim([1 18]); xlim([0 600]);
    ax = gca;
    ax.FontSize = 15;
end


%% L0 Analysis.
if showL0 == true
    
    % HISTOGRAM.
    chla((chla==0)) = nan;
    figure
    histfit(chla(depth_B==2),10,"lognormal");
    
    % Calculate the mean depth of the Chlorophyll Maximum (CM) as well as the
    % 5th and 9th percentile interval values.
    meanCM = mean(dcmBats(:,2),"omitnan");
    CM_5pct = prctile(dcmBats(:,2),5);
    CM_95pct = prctile(dcmBats(:,2),95);
    
    % Use the Anderson-Darling (A-D) test to evaluate whether the data is
    % distributed normally or lognormally. The hypothesis test result 'h' will
    % return as h = 1 if the null hypothesis is rejected or h = 0 if there is a
    % failure to reject the null hypothesis.
    hN = nan(20,1); pN = nan(20,1); hL = nan(20,1); pL = nan(20,1);
    obs = nan(20,1);
    for i = 1:20
        tmp = chla(depth_B==i);
        if length(tmp) > 30
            obs(i) = length(tmp);
            tmp(tmp==0) = nan;
    
            [hN(i), pN(i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
            [hL(i), pL(i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);
    
            if showOtherTests == true
                pd0 = fitdist(tmp,'Normal');
                pd = fitdist(tmp,'Lognormal');
                [hN2(i), pN2(i)] = chi2gof(tmp,"CDF",pd0);
                [hN3(i), pN3(i)] = lillietest(tmp,"Distribution","norm");
                [hx1(i),px1(i)] = chi2gof(tmp,'CDF',pd);
                [hx2(i),px2(i)] = lillietest(log(tmp),"Distr","norm");
            end
        end
    end
    
    % Plot Figure 2. chl-a. Is the data normal or lognormal?
    figure;
    if showL0title == true
        sgtitle("chl-a (L0): " + "BATS "+num2str(YMD(1))+" - " + num2str(YMD(end))+"");
    end
    subplot(1,3,[1 2])
    yyaxis left
    semilogx(pN,0.5:1:19.5,'o-','Color','#c51b7d','DisplayName','Normal (A-D)','LineWidth',1.5,'MarkerSize',5);
    hold on
    semilogx(pL,0.5:1:19.5,'o-','Color','#4d9221','DisplayName','Lognormal (A-D)','LineWidth',1.5,'MarkerSize',5);
    if showOtherTests == true
        semilogx(pN2,1:1:20,'o-','Color','#c51b7d','DisplayName','Normal (chi^2)','LineWidth',1.5,'MarkerSize',5);
        semilogx(pN3,1:1:20,'o--','Color','#c51b7d','DisplayName','Normal (Lil.)','LineWidth',1.5,'MarkerSize',5);
        semilogx(px1,1:1:20,'o-','Color','#4d9221','DisplayName','Lognormal (chi^2)','LineWidth',1.5,'MarkerSize',5);
        semilogx(px2,1:1:20,'o--','Color','#4d9221','DisplayName','Lognormal (Lil.)','LineWidth',1.5,'MarkerSize',5);
    end
    set(gca,"YDir","reverse"); legend();
    yticklabels(0:20:200);
    ylim([0 20]);
    ylabel("Depth (m) (10-m bins)");
    yyaxis right
    yline(meanCM,DisplayName="p_{DCM} \pm 5/95",Interpreter="latex");
    yline(CM_5pct,HandleVisibility="off");
    yline(CM_95pct,HandleVisibility="off");
    xline(0.005,":","\alpha","DisplayName","\alpha = 0.005");
    hold off
    set(gca,"YDir","reverse"); legend();
    yticklabels({}); grid on;
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    ylim([0 200]);
    xlim([0.5e-3 1]);
    xlabel("p-value");
    
    subplot(1,3,3)
    barh(obs,'FaceColor','#d3d3d3');
    hold on
    xline(30);
    set(gca,"YDir","reverse"); xlabel("No. of Obs.");
    ylim([0.5 20.5]); yticklabels({});
end

%% Find Mixed Layer Depth
% Calculate the Mixed Layer Depth (MLD) from hydrographic data. This is
% needed for separating the L0 analysis into L1 and L2.

% Import data.
data = importdata("data\bats_bottle.txt").data;
id = data(:,1);     % Bottle ID with format !####$$$@@, where ! = cruise
                    % type (1 = core cruise), #### = cruise number, 
                    % $$$ = cast number, @@ = Niskin number.
YMD = data(:,2);                    % Year, Month, Day.
lat = data(:,5); long = data(:,6);  % Latitude and Longitude.
z = data(:,8);                      % Depth (m).
T = data(:,9);                      % Temperature (C).
Sp = data(:,11);                    % Practical Salinity.

% Setup time vector.
tmpB = num2str(YMD);
YY = str2num(tmpB(:,1:4));
MM = str2num(tmpB(:,5:6));
DD = str2num(tmpB(:,7:8));
t = datetime(YY,MM,DD);

% (re-) assign NaNs.
z(z==-999) = nan;

% Process Hydrography data according to TEOS-10 standard.
p = gsw_p_from_z(-z,lat);
SA = gsw_SA_from_SP(Sp,p,long,lat);
CT = gsw_CT_from_t(SA,T,p);

% Import cruise type "cType" and cruise number "CRN". 
tmp3 = num2str(id);
cType = str2num(tmp3(:,1));
CRN = str2num(tmp3(:,2:5));
clear tmp3;

% Extract hydrography from the core cruises (labelled as cType = 1).
ids = find(cType==1);
CoreCRN = CRN(ids);
SA = SA(ids);
CT = CT(ids);
p = p(ids);

% Save the point at which a new cruise starts to the array 'id_nc'.
% Additionally, save the actual cruise number of the new cruise to the
% array 'crn_mr'. This is necessary because not all cruises have MLD data.
id_nc = [1];    % ID (id_) of new cruise (nc).
crn_mr = [1];   % cruise number (crn_) where MLD was recorded (mr).
for i = 2:length(CoreCRN)
    if CoreCRN(i) > CoreCRN(i-1)
        id_nc = [id_nc i];
        crn_mr = [crn_mr CoreCRN(i)];
    end
end
t_nc = t(id_nc);    % Time at which new cruise starts.

% Calculate the Mixed Layer Depth per cruise 'MLD_pc'
mld_pc = 20*ones(405,1);        % Minimum possible MLD = 20 dbar.
mld_pc(1) = gsw_mlp(SA(1:24),CT(1:24),p(1:24));
for i = 2:159
    mld_pc(i) = gsw_mlp(SA(id_nc(i):id_nc(i+1)-1),CT(id_nc(i):id_nc(i+1)-1),p(id_nc(i):id_nc(i+1)-1));
end
mld_pc(161) = gsw_mlp(SA(14613:14697),CT(14613:14697),p(14613:14697));
for i = 161:403
    mld_pc(i+1) = gsw_mlp(SA(id_nc(i):id_nc(i+1)-1),CT(id_nc(i):id_nc(i+1)-1),p(id_nc(i):id_nc(i+1)-1));
end
mld_pc(405) = gsw_mlp(SA(64123:end),CT(64123:end),p(64123:end));

% Remove casts that are unrealistic.
mld_pc(mld_pc>500) = nan;
% mld_pc(isnan(mld_pc)) = 20; % reset removed cast as minimum MLD possible


%% Compare depths: DCM vs ML
% Compare the depth of the DCM and the ML for each cruise. We are
% interested in those cruises where the DCM is deeper than the ML.

% Convert DCM depth (m) to DCM pressure (dbar).
p_DCM = gsw_p_from_z(-dcmBats(:,2),lat(1));

% Which cruises to look at with Stn ALOHA Level Methodology. We need the
% DCM to be beneath the ML.
% Consolidate two vectors of MLD and DCM so they have the same size.
cruises = 1:1:405;                          % number per cruise
newMLDcrnVector = nan(length(cruises),1);
newMlp = nan(length(cruises),1);

for i = 1:405
    for k = 1:404
        if crn_mr(k) == i
            newMLDcrnVector(i) = crn_mr(k);
            newMlp(i) = mld_pc(k);
        end
    end
end

if showDcmMldComp == true
    figure;
    plot(dcmBats(:,1),dcmBats(:,2),DisplayName="Chl-a Maximum");
    hold on
    plot(newMLDcrnVector,newMlp,DisplayName="Mixed Layer Depth");
    hold off
    set(gca,"YDir","reverse");
    ylabel("Pressure [dbar]"); xlabel("Cruise No.");
    legend();
end

% Calculate where DCM is beneath the ML. This is important because it is on
% these cruises where we apply our L1/L2 analysis.

dcmBelow = dcmBats(:,2) - newMlp;
cruisesWhereDCMisBelowMLD = [];

for i = 1:405
    if dcmBelow(i) > 0
        cruisesWhereDCMisBelowMLD = [cruisesWhereDCMisBelowMLD i];
    end
end

if showDcmMldComp == true
    figure
    plot(1:1:405,dcmBelow);
end

% add column to start
newChla = [CRN_no' chlagrid];

test1 = 1:1:405;
test1 = test1';
tmpC = nan(length(test1),35); tmpD = nan(length(test1),35); tmpT = nan(length(test1),35);
for i = 1:405
    for j = 1:398
        if CRN_no(j) == i
            tmpC(i,:) = chlagrid(j,:);
            tmpD(i,:) = dgrid(j,:);
            tmpT(i,:) = datenum(tgrid(j,:));
        end
    end
end

test3 = [test1 tmpC]; % this has forced the chlagrid array onto a 405 cruise matrix.
test4 = [test1 tmpD]; % same for depth
test5 = [test1 tmpT]; % ... and time

% Now show only the chla arrays where DCM is beneath the MLD.
newChlaArray = test3(cruisesWhereDCMisBelowMLD,:);
newDepthArray = test4(cruisesWhereDCMisBelowMLD,:);
newTimeArray = test5(cruisesWhereDCMisBelowMLD,:);

% remember that the first entry is the CRUISE NO.

dcmsBelow = dcmBats(cruisesWhereDCMisBelowMLD,:);

% try to work out which bottles are from these cruises where the DCM is
% below the MLD.

test66 = [];
for i = 1:length(copyOfCoreCrn)
    for j = 1:length(cruisesWhereDCMisBelowMLD)
        if copyOfCoreCrn(i) == cruisesWhereDCMisBelowMLD(j)
            test66 = [test66 i];
        end
    end
end

chla_lowDCM = chla(test66);
dep_lowDCM = depth_B(test66);

% Calculate the mean depth of the Chlorophyll Maximum (CM) as well as the
% 5th and 9th percentile interval values.
meanCM = mean(dcmsBelow(:,2),"omitnan");
CM_5pct = prctile(dcmsBelow(:,2),5);
CM_95pct = prctile(dcmsBelow(:,2),95);


% Use the Anderson-Darling (A-D) test to evaluate whether the data is
% distributed normally or lognormally. The hypothesis test result 'h' will
% return as h = 1 if the null hypothesis is rejected or h = 0 if there is a
% failure to reject the null hypothesis.
hN = nan(20,1); pN = nan(20,1); hL = nan(20,1); pL = nan(20,1);
obs = nan(20,1);
for i = 1:20
    tmp = chla_lowDCM(dep_lowDCM==i);
    if length(tmp) > 30
        obs(i) = length(tmp);
        tmp(tmp==0) = nan;

        [hN(i), pN(i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
        [hL(i), pL(i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);

        if showOtherTests == true
            pd0 = fitdist(tmp,'Normal');
            pd = fitdist(tmp,'Lognormal');
            [hN2(i), pN2(i)] = chi2gof(tmp,"CDF",pd0);
            [hN3(i), pN3(i)] = lillietest(tmp,"Distribution","norm");
            [hx1(i),px1(i)] = chi2gof(tmp,'CDF',pd);
            [hx2(i),px2(i)] = lillietest(log(tmp),"Distr","norm");
        end
    end
end
%% L0 Analysis.
if showL0 == true
    % A-D test: L0. (Fig. 2) [DCM deeper than ML]
    % Figure 2X. chl-a. L0. Is the data normal or lognormal? This time we look
    % at only those cruises where the DCM was beneath the MLD.
    
    figure;
    if showL0title == true
        sgtitle("chl-a (L0): " + "BATS "+num2str(YMD(1))+" - " + num2str(YMD(end))+"");
    end
    subplot(1,3,[1 2])
    yyaxis left
    semilogx(pN,0.5:1:19.5,'o-','Color','#c51b7d','DisplayName','Normal (A-D)','LineWidth',1.5,'MarkerSize',5);
    hold on
    semilogx(pL,0.5:1:19.5,'o-','Color','#4d9221','DisplayName','Lognormal (A-D)','LineWidth',1.5,'MarkerSize',5);
    if showOtherTests == true
        semilogx(pN2,1:1:20,'o-','Color','#c51b7d','DisplayName','Normal (chi^2)','LineWidth',1.5,'MarkerSize',5);
        semilogx(pN3,1:1:20,'o--','Color','#c51b7d','DisplayName','Normal (Lil.)','LineWidth',1.5,'MarkerSize',5);
        semilogx(px1,1:1:20,'o-','Color','#4d9221','DisplayName','Lognormal (chi^2)','LineWidth',1.5,'MarkerSize',5);
        semilogx(px2,1:1:20,'o--','Color','#4d9221','DisplayName','Lognormal (Lil.)','LineWidth',1.5,'MarkerSize',5);
    end
    set(gca,"YDir","reverse"); legend();
    yticklabels(0:20:200);
    ylim([0 20]);
    ylabel("Pressure [dbar]",Interpreter='latex',FontSize=13);
    yyaxis right
    yline(meanCM,DisplayName="p_{DCM} \pm 5/95",LineWidth=1,Color="#808080",Interpreter="latex");
    yline(CM_5pct,LineWidth=1,Color="#808080",HandleVisibility="off");
    yline(CM_95pct,LineWidth=1,Color="#808080",HandleVisibility="off");
    xline(0.005,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
    hold off
    set(gca,"YDir","reverse"); legend();
    yticklabels({});
    ylim([0 200]);
    xlim([0.5e-3 1]);
    xlabel("A-D p-value",Interpreter='latex',FontSize=13);
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    legend('Position',[0.4 0.7 0.07 0.12],FontSize=11);
    grid on
    sgtitle("L0","Interpreter","latex");
    
    subplot(1,3,3)
    barh(obs,'FaceColor','#d3d3d3');
    hold on
    xline(30);
    set(gca,"YDir","reverse"); xlabel("No. of Obs.",Interpreter='latex',FontSize=13);
    ylim([0.5 20.5]); yticklabels({});

end

%% L1 A-D. (Fig. 3a)
if showL1 == true

% Set whether to analyse total data or...
% only data for cruises where the DCM is beneath the ML 
dcmBeneathML = true;

% Import data: bottle ID, depth, and chl-a concentration.
idIn = D(:,1);
depthIn = D(:,8);
cIn = D(:,22);

% Convert depth => pressure.
pIn = gsw_p_from_z(-depthIn,lat(1)); clear depthIn;

% Extract cruise no.
tmp = num2str(idIn);
cruiseNo = str2num(tmp(:,2:5)); clear tmp;

% Look only at cruises where DCM is beneath ML.
if dcmBeneathML == true
    disp("L1: examining only cruises where DCM is beneath ML")
    Lia = ismember(cruiseNo,cruisesWhereDCMisBelowMLD);
    cruiseNo(Lia==0) = [];
    pIn(Lia==0) = [];
    cIn(Lia==0) = [];
else
    disp("L1: examining total dataset...")
end

% Populate new arrays with only Mixed-Layer values.
L = length(pIn);
tmpP = nan(1,L);
tmpCrn = nan(1,L);
tmpC = nan(1,L);
tmpId = nan(1,L);

for i = 1:L
    disp(i)
    tmp = mld_pc(cruiseNo(i));
    if pIn(i) < tmp
        tmpP(i) = pIn(i);
        tmpCrn(i) = cruiseNo(i);
        tmpC(i) = cIn(i);
        tmpId(i) = idIn(i);
    end
end

% Remove nan values (i.e. values beneath Mixed Layer).
pOut = tmpP(~isnan(tmpP));
crnOut = tmpCrn(~isnan(tmpCrn));
COut = tmpC(~isnan(tmpC));
idOut = tmpId(~isnan(tmpId));

% Remove bottles that are too close to the surface (< 2.5 dbar)
idRm = pOut > 2.5;
pOut = pOut(idRm);
COut = COut(idRm);
idOut = idOut(idRm);

% Remove bottles where concentration of COut = 0
idZero = COut == 0;
pOut = pOut(~idZero);
COut = COut(~idZero);
idOut = idOut(~idZero);

% Bin pressures.
pb10 = discretize(pOut,0:10:200);
n10 = max(pb10);

%% L1 Histogram.

COut(COut<=0) =nan;
figure;
histogram(COut(pb10==1));

% histfit variation
figure
histfit(COut(pb10==4),10,"lognormal");

%% L1 A-D plot.
obs = nan(20,1);
n = n10;
depth = 5:10:200;
ad = nan(2,20);

% Apply A-D test to values of chl-a in the mixed layer.
for i = 1:n
    X_i = COut(pb10==i); % find concentration X_i at binned pressure i
    if length(X_i) > 3
        [~,ad(2,i)] = adtest(X_i,"Distribution","logn","Alpha",0.005);
        [~,ad(1,i)] = adtest(X_i,"Distribution","norm","Alpha",0.005);
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

for i = 1:n
    if obs(i) < threshold
        ad(:,i) = nan;
        sk(i) = nan;
        ku(i) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
    end
end

tmp = [];
for i = 1:n
    if ~isnan(sum(ad(:,i)))
        tmp = [tmp i];
    end
end
tr = depth(tmp);
sk = sk(tmp);
ku = ku(tmp);
rV = rV(:,tmp);
pV = pV(:,tmp);
ad = ad(:,~all(isnan(ad)));

% Plot results.
ax = figure;

vuongRes = nan(length(tr));
for i = 1:length(tr)
    if rV(1,i) > 0
        %disp('Normal');
        vuongRes(i) = 1;
    elseif rV(1,i) < 0
        %disp('Lognormal');
        vuongRes(i) = 2;
    end
end

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
skLognN = -skLogn;
kuLognN = kuLogn;

% Specify expectation value.
alphaHy = 0.005;
alphaLlr = 0.005;

subplot(1,3,3)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
yticks(0.5:1:9.5);
yticklabels(0:10:90);
ylim([0.5 9.5]);
yticklabels({});
xlabel('No. of Obs.',Interpreter='latex',FontSize=13);

subplot(1,3,[1 2])
xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
hold on
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
plot(nan,nan,'square','Color','#808080','MarkerSize',15,'DisplayName','V-LLR best fit (p > 0.1)');        
plot(nan,nan,'square','Color','#808080','MarkerSize',15,'LineWidth',4,'DisplayName','V-LLR best fit (p < 0.1)');        
plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
hold off
set(gca, 'XScale', 'log');
grid minor;
xlim([0.1*alphaHy 1]); 
ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
set(gca,'YDir','reverse');
legend('Position',[0.4 0.7 0.07 0.12],FontSize=11);
sgtitle("L1","Interpreter","latex");

%% L1 Skewness/Kurtosis. (Fig. 4a)

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
% Negative distributions. For plotting purposes only.
skLognN = -skLogn;
kuLognN = kuLogn;

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

ax = figure;
scatter(nan,nan,72,[0.8 0.8 0.8],DisplayName='Data');
hold on
scatter(0,3,72,[0.2 0.2 0.2],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(sk,ku,72,[0.8 0.8 0.8],HandleVisibility="off");
clr = 1:1:length(tr);
scatter(sk,ku,54,clr,"filled","o",HandleVisibility="off");
colormap(gca,cbrewer2("Greens"));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(tr);
cbar.TickLabels = tr;
cbar.Label.String = "P [dbar]";
cbar.Label.Position = [0.7 1-0.35];
cbar.Label.Rotation = 0;

hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness','FontSize',13,'Interpreter','latex'); ylabel('Kurtosis',FontSize=13,Interpreter='latex');
lgd = legend('Location','best');
title('L1','Interpreter','latex','FontSize',13);

end
%% L2 A-D. (Fig. 3b)
if showL2 == true
% Set whether to analyse total data or...
% only data for cruises where the DCM is beneath the ML 
dcmBeneathML = true;

% Import data.
idIn = num2str(D(:,1));
depthIn = D(:,8);
cIn = D(:,22);

% Convert depth => pressure.
pIn = gsw_p_from_z(-depthIn,lat(1)); clear depthIn;

% (included code from extractSMLC)
% Exclude NaN data (-9 = NaN)
cIn(cIn==-9) = nan;
idIn = idIn(~isnan(cIn),:);
pIn = pIn(~isnan(cIn));
cIn = cIn(~isnan(cIn));

% Cruise number 'crn'
crn = str2num(idIn(:,2:5));

% Look only at cruises where DCM is beneath ML.
if dcmBeneathML == true
    disp("L2: analysing only cruises where DCM is beneath ML...");
    Lia = ismember(crn,cruisesWhereDCMisBelowMLD);
    crn(Lia==0) = [];
    pIn(Lia==0) = [];
    cIn(Lia==0) = [];
else
    disp("L2: analysing complete dataset...");
end
% Extract measurements below pMaxMld
L = length(pIn);
tmpP_subML = nan(L,1);
tmpCrn_subML = nan(L,1);
tmpC_subML = nan(L,1);
botID = [];

% To find the measurements taken beneath pMaxMld, we populate the arrays
% just defined for cases where p > pMld...
for i = 1:L
    tmpMld = mld_pc(crn(i));
    if p(i) > tmpMld
        tmpP_subML(i) = pIn(i);
        tmpCrn_subML(i) = crn(i);
        tmpC_subML(i) = cIn(i);
        tmpStr = idIn(i,:);
        botID = [botID;tmpStr];
    end
end

% ...and remove nan values (which represent measurements above pMaxMld)
pOut = tmpP_subML(~isnan(tmpP_subML));
cOut = tmpC_subML(~isnan(tmpC_subML));
idOut = botID;

% Extract cruise number 'crn'
crn = str2num(idOut(:,2:5)); 

bottleArray = [crn pOut];

% Create an array of all unique bottle cruise/cast combinations
botCrnCast = rmmissing(unique(bottleArray(:,1),"rows"));

dcmCrnCast = [];
for i = 1:length(dcmBats(:,1))
    for x = 1:length(botCrnCast)
        if dcmBats(i,1) == botCrnCast(x,1) 
            dcmCrnCast = [dcmCrnCast i];
        end
    end
end

% Split bottle concentration by cruise & cast
tid = [];
for i = 2:length(pOut)
    % check if CRN changes
    if bottleArray(i,1) > bottleArray(i-1,1)
        tid = [tid i];
    end
end

tPcm = nan(length(pOut),1);
tPcm(1:tid(1)-1) = dcmBats(dcmCrnCast(1),2);
tPcm(tid(end):end) = dcmBats(dcmCrnCast(end),2);

Ltid = length(tid);

tmp = length(dcmCrnCast) - 1;
Ltid = tmp;
for i = 2:Ltid-2
    tPcm(tid(i):tid(i+1)-1) = dcmBats(dcmCrnCast(i),2);
end

bottleArray = [bottleArray tPcm];

% THIS is where we convert to Lagrangian pressure coordinates!!!
tPLagrangian = nan(length(p),1);
tPLagrangian = bottleArray(:,2) - bottleArray(:,3);
bottleArray = [bottleArray tPLagrangian];

pB10 = round(bottleArray(:,4),-1);
bottleArray = [bottleArray pB10];

bottleArray = [bottleArray cOut];
pmin = min(bottleArray(:,5));
pmax = max(bottleArray(:,5));
pr = pmin:10:pmax;

C_out = bottleArray(:,6);
pB = bottleArray(:,5);
ad = nan(2,length(pr));
rV = nan(10,length(pr));
pV = nan(10,length(pr));
obs = nan(1,length(pr));

sk = nan(1,length(pr));
ku = nan(1,length(pr));

for i = 1:length(pr)
    tmp = C_out(pB==pr(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
    end
end

for i = 1:length(pr)
    if obs(i) < threshold
        ad(:,i) = nan;
        rV(:,i) = nan;
        sk(i) = nan;
        ku(i) = nan;
    end
end

% 3.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = zeros(1,length(pr));
rV(isnan(rV)) = 0;

for i = 1:length(pr)
    if rV(1,i) > 0
        vuongRes(i) = 1;
    elseif rV(1,i) < 0
        vuongRes(i) = 2;
    end
end
rV(rV==0) = nan;
%% L2 Histograms.
% pB10, C_out

C_out(C_out < 0) = nan;

figure;
histogram(C_out(pB10==-100));

% histfit variation
% figure
% histfit(C_out(pB10==10),10,"lognormal");

%% L2 A-D. Plot.
% plotks2 code.
figure
n = length(pr);

alphaHy = 0.005;
alphaLlr = 0.1;

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
skLognN = -skLogn;
kuLognN = kuLogn;

subplot(1,3,3)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
yticklabels({});
xlabel('No. of Obs.',Interpreter='latex',FontSize=13);
% ylim([1 1+b-a]);
ylim([14 26]);

subplot(1,3,[1 2])
xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
hold on
for i = 1:n
    if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr && ad(2,i) > alphaHy
        plot(ad(1,i),pr(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
    elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr && ad(2,i) > alphaHy
        plot(ad(1,i),pr(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
    elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr && ad(1,i) > alphaHy
        plot(ad(2,i),pr(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
    elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr && ad(1,i) > alphaHy
        plot(ad(2,i),pr(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
    end
end
plot(ad(1,:),pr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5,'HandleVisibility','off');
plot(ad(2,:),pr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5,'HandleVisibility','off');
xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
set(gca, 'XScale', 'log');
hold off
grid minor;
ylim([-60 60]);
xlim([0.1*alphaHy 1]);
set(gca,'YDir','reverse');
ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
sgtitle("L2","Interpreter","latex");

%% L2 Skewness/Kurtosis. (Fig. 4b)

% Output only skewness/kurtosis values at depths defined by these limits.
% limits = [-60 60]
a = 14; b = 26;
tmp = [];
for i = a:b
    if ~isnan(sum(ad(:,i)))
        tmp = [tmp i];
    end
end
pr = pr(tmp);
sk = sk(tmp);
ku = ku(tmp);
clear tmp;

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

figure
scatter(0,3,72,[0.2 0.2 0.2],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(sk,ku,72,[0.8 0.8 0.8],HandleVisibility='off');
clr = 1:1:length(pr);
scatter(sk,ku,54,clr,"filled","o",HandleVisibility="off");
colormap(gca,flipud(cbrewer2("PiYG")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(pr);
cbar.TickLabels = pr;
cbar.Label.String = "P [dbar]";
cbar.Label.Position = [0.7 1-0.7];
cbar.Label.Rotation = 0;
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness',FontSize=13,Interpreter='latex'); 
ylabel('Kurtosis',FontSize=13,Interpreter='latex');
title('L2','Interpreter','latex','FontSize',13);

end