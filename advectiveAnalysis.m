% Advective Analysis for Station ALOHA.
% Idea -> Draw a cross around Station ALOHA. Calculate the gradients
% across Stn ALOHA in the horizontal (dC/dx) and vertical (dC/dy). Multiply
% these by the zonal (u) and meridional (v) components respectively of
% velocity at Stn ALOHA. This will give us a time rate change in
% concentration which we can directly compare to (Eulerian) static
% concentration change at Stn ALOHA.

%% Setup. Test Cases.
clear; clc; close all; addpath("output\");
set(groot,'defaultAxesXGrid','on'); set(groot,'defaultAxesYGrid','on');
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(0,'defaultAxesFontSize',10);

% Use subtightplot.
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.085 0.115], [0.1 0.05]);
if ~make_it_tight,  clear subplot;  end

% Test Cases.
readOscar = false;          % Set to TRUE if you need to read in the OSCAR 
                            % data. Otherwise leave as FALSE because this 
                            % will take some time.
regio = false;              % Set to TRUE to extract velocities surrounding
                            % STN ALOHA (instead of @Stn ALOHA).
fourTermComp = false;       % Set to TRUE to show a comparison of the 
                            % magnitude of the four terms in Case 2 (g2).
compareCaseEsts = false;    % Set to TRUE to compare the three different 
                            % methods (Cases) of determining advection. 
                            % NOTE that Case 3 is faulty...
compareTermsInG3 = false;   % Set to TRUE to show plot comparing size of 
                            % terms in Case 3 (g3).
showFft = false;             % Show FFT of CHL at Stn ALOHA.
applyCorrelations = false;  % Set to TRUE to calculate correlation between 
                            % total derivative, advection, and local
                            % derivative.
histfitOrGram = false;      % If TRUE, use histfit. Else use histogram.
                            % Histfit overlays a distribution whereas
                            % histogram can have uniform edges for easier
                            % comparison.
ptEast = false;              % Apply advective analysis to point 500km east 
                            % of Stn ALOHA.
%% View OSCAR.

info = ncinfo("..\..\Large Data Files\OSCAR2\oscar_currents_final_19970408.nc");
lon = ncread("../../Large Data Files/OSCAR2/oscar_currents_final_19970408.nc","lon");

%% Read in velocities from OSCAR Reanalysis.

if readOscar == true
    latAL = 451; lonAL = 829;       % Stn ALOHA location. 22.75N 158W. Lat = 451. Lon = 809.
    % 500km east of ALOHA: 153W. Lon= 829
    starT = datetime(1997,04,08);   % Start date (to coincide with 
                                    % earliest data from GLO_BGC).
    enD = 9249;                     % End date. 04/08/2022.
    endDD = starT + days(enD-1);
    suffix = "../../Large Data Files/OSCAR2/oscar_currents_final_";
    
    % Automate extraction.
    for i = 1:176
        FID = num2str(yyyymmdd(starT + days(i)));
        fname = strcat(suffix,FID,".nc");
        T1 = ncread(fname,"u");
        T2 = ncread(fname,"v");
        T3 = ncread(fname,"time");
        if regio == true
            u(i,:) = [T1(latAL+1,lonAL) T1(latAL,lonAL+1) T1(latAL-1,lonAL) T1(latAL,lonAL-1)];
            v(i,:) = [T2(latAL+1,lonAL) T2(latAL,lonAL+1) T2(latAL-1,lonAL) T2(latAL,lonAL-1)];
        else
            u(i) = T1(latAL,lonAL);
            v(i) = T2(latAL,lonAL);
        end
        t(i) = datetime(1990,01,01) + days(T3);
        clear fname T1 T2 T3;
    end
    
    % ERROR opening 19971002. Replace this value with NaNs.
    if regio == true
        u(177,:) = NaN; v(177,:) = NaN;
    else
        u(177) = NaN; v(177) = NaN;
    end
    t(177) = datetime(1990,01,01) + days(177);
    
    % Start again 19971003.
    for i = 178:255
        FID = num2str(yyyymmdd(datetime(1997,04,08) + days(i)));
        fname = strcat(suffix,FID,".nc");
        T1 = ncread(fname,"u");
        T2 = ncread(fname,"v");
        T3 = ncread(fname,"time");
        if regio == true
            u(i,:) = [T1(latAL+1,lonAL) T1(latAL,lonAL+1) T1(latAL-1,lonAL) T1(latAL,lonAL-1)];
            v(i,:) = [T2(latAL+1,lonAL) T2(latAL,lonAL+1) T2(latAL-1,lonAL) T2(latAL,lonAL-1)];
        else
            u(i) = T1(latAL,lonAL);
            v(i) = T2(latAL,lonAL);
        end
        t(i) = datetime(1990,01,01) + days(T3);
        clear fname T1 T2 T3;
    end
    
    % ERROR opening 19971220. Replace with NaNs.
    if regio == true
        u(256,:) = NaN; v(256,:) = NaN;
    else
        u(256) = NaN; v(256) = NaN;
    end
    t(256) = datetime(1990,01,01) + days(256);
    
    % Start again 19971221.
    for i = 257:4479
        FID = num2str(yyyymmdd(datetime(1997,04,08) + days(i)));
        fname = strcat(suffix,FID,".nc");
        T1 = ncread(fname,"u");
        T2 = ncread(fname,"v");
        T3 = ncread(fname,"time");
        if regio == true
            u(i,:) = [T1(latAL+1,lonAL) T1(latAL,lonAL+1) T1(latAL-1,lonAL) T1(latAL,lonAL-1)];
            v(i,:) = [T2(latAL+1,lonAL) T2(latAL,lonAL+1) T2(latAL-1,lonAL) T2(latAL,lonAL-1)];
        else
            u(i) = T1(latAL,lonAL);
            v(i) = T2(latAL,lonAL);
        end
        t(i) = datetime(1990,01,01) + days(T3);
        clear fname T1 T2 T3;
    end
    
    % ERROR opening 20090714. Set as NaN.
    if regio == true
        u(4480,:) = NaN; v(4480,:) = NaN;
    else
        u(4480) = NaN; v(4480) = NaN;
    end
    t(4480) = datetime(1990,01,01) + days(4480);
    
    % Start again 20090715.
    for i = 4481:6000
        FID = num2str(yyyymmdd(datetime(1997,04,08) + days(i)));
        fname = strcat(suffix,FID,".nc");
        T1 = ncread(fname,"u");
        T2 = ncread(fname,"v");
        T3 = ncread(fname,"time");
        if regio == true
            u(i,:) = [T1(latAL+1,lonAL) T1(latAL,lonAL+1) T1(latAL-1,lonAL) T1(latAL,lonAL-1)];
            v(i,:) = [T2(latAL+1,lonAL) T2(latAL,lonAL+1) T2(latAL-1,lonAL) T2(latAL,lonAL-1)];
        else
            u(i) = T1(latAL,lonAL);
            v(i) = T2(latAL,lonAL);
        end
        t(i) = datetime(1990,01,01) + days(T3);
        clear fname T1 T2 T3;
    end
    
    % ERROR opening 20130912. Set as NaN.
    if regio == true
        u(6001,:) = NaN; v(6001,:) = NaN;
    else
        u(6001) = NaN; v(6001) = NaN;
    end
    t(6001) = datetime(1990,01,01) + days(6001);
    
    % Start again 20130913.
    for i = 6002:enD
        FID = num2str(yyyymmdd(datetime(1997,04,08) + days(i)));
        fname = strcat(suffix,FID,".nc");
        T1 = ncread(fname,"u");
        T2 = ncread(fname,"v");
        T3 = ncread(fname,"time");
        if regio == true
            u(i,:) = [T1(latAL+1,lonAL) T1(latAL,lonAL+1) T1(latAL-1,lonAL) T1(latAL,lonAL-1)];
            v(i,:) = [T2(latAL+1,lonAL) T2(latAL,lonAL+1) T2(latAL-1,lonAL) T2(latAL,lonAL-1)];
        else
            u(i) = T1(latAL,lonAL);
            v(i) = T2(latAL,lonAL);
        end
        t(i) = datetime(1990,01,01) + days(T3);
        clear fname T1 T2 T3;
    end
    
    save output/alohaPlus500 t u v;
end

%% Load velocities + CHL.

% Load OSCAR Reanalysis velocities at Stn ALOHA, ...
t = load("output\oscar.mat","t").t;
t([177 256 4480 6001]) = NaT;
u = load("output\oscar.mat","u").u;
v = load("output\oscar.mat","v").v;
% ... and at four nearest points around Stn ALOHA.
uR = load("output\alohaRegVels.mat","u").u;
vR = load("output\alohaRegVels.mat","v").v;
% Load velocities at 500km east of Stn ALOHA.
uE = load("output\alohaPlus500.mat").u;
vE = load("output\alohaPlus500.mat").v;

% Load GLO_BGC satellite CHL data (L4 daily interpolated "finalised").
input3 = "../../Large Data Files/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_1759494520147.nc";

% Load GLO_BGC satellite CHL data (L4 daily) for pt 500 km E of Stn ALOHA.
input4 = "../../Large Data Files/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_1762786251875.nc";

% % UNCOMMENT to view variables, dimensions, etc.
% vf = ncinfo(input3);
% % Variables (dim.) = time (10249), latitude (96), longitude(100), CHL,
% % CHL_uncertainty

% Load variables from GLO_BGC L4 CHL satellite product.
tmp = ncread(input3,"time");
time = datetime(tmp,"ConvertFrom","posixtime"); % Units of days.
lat = ncread(input3,"latitude");  % 0.0417 deg. spacing.
lon = ncread(input3,"longitude"); % 0.0417 deg. spacing.
CHL = ncread(input3,"CHL");        % Format (LON,LAT,TIME). Units: mg.m-3.
A = 1; B = 9101;                    % Time limits. (1 = Sep 04 1997, 
                                    % 9101 = Aug 04 2022)

if ptEast == true
    % Load variables from GLO_BGC L4 CHL satellite product for POINT 500km east
    % of Stn ALOHA.
    tmp = ncread(input4,"time");
    time = datetime(tmp,"ConvertFrom","posixtime"); % Units of days.
    lat = ncread(input4,"latitude");  % 0.0417 deg. spacing.
    lon = ncread(input4,"longitude"); % 0.0417 deg. spacing.
    CHL = ncread(input4,"CHL");        % Format (LON,LAT,TIME). Units: mg.m-3.
    A = 1; B = 9101;                    % Time limits. (1 = Sep 04 1997, 
                                    % 9101 = Aug 04 2022)
end
%% Calculate DC/Dt + u.GRAD(c).

% Evaluate at point EAST of Stn ALOHA. Otherwise (if false) evaluate at Stn
% ALOHA.
if ptEast == true
    CHL_study = squeeze(CHL(24,24,A:B));
else
    % CHL at Stn ALOHA. Average of four nearest points. [Units = mg.m-3].
    CHL_study = squeeze(mean([CHL(50,49,A:B) CHL(51,49,A:B) CHL(50,48,A:B) CHL(51,48,A:B)],2,"omitnan"));

end

% edges = -3e-7:2.5e-8:3e-7;
figure;
histogram(CHL_study); %,BinEdges=edges);le
legend("CHL",Location="best");
if ptEast == true
    title("22.75N 153W: Satellite CHL 1997-2022","Fontsize",10,Interpreter="latex");
else
    title("Stn ALOHA: Satellite CHL 1997-2022","Fontsize",10,Interpreter="latex");
end
xlabel("Concentration (mg m$^{-3}$)",Interpreter="latex");
ylabel("No. of Counts",Interpreter="latex");

% dc/dt at Station ALOHA. 
% Divide by 86400 to convert from [day-1] to [s-1].
% [Units = mg.m-3.s-1].
% Label EUL for Eulerian (Local Change).
EUL = (1/86400) * diff(CHL_study); 

% CHL Gradient in X and Y. [Units of mg.m-3.m-1. = mg.m-4].
if ptEast == true
    gradCX = squeeze ( 0.5*( (CHL(25,24,A:B) - CHL(24,24,A:B)) + (CHL(25,25,A:B) - CHL(24,25,A:B)) ) ./ 4266);
    gradCY = squeeze ( 0.5*( (CHL(24,25,A:B) - CHL(24,24,A:B)) + (CHL(25,25,A:B) - CHL(25,24,A:B)) ) ./ 4266 );
else
    gradCX = squeeze ( 0.5*( (CHL(51,49,A:B) - CHL(50,49,A:B)) + (CHL(51,48,A:B) - CHL(50,48,A:B)) ) ./ 4266);
    gradCY = squeeze ( 0.5*( (CHL(51,49,A:B) - CHL(51,48,A:B)) + (CHL(50,49,A:B) - CHL(50,48,A:B)) ) ./ 4626 );
end
% Case 1. Advective Change in Concentration [Units of mg.m-3.s-1].
% c(du/dx + dv/dy) assumed to be zero. Label as ADV for ADVection.
ADV_1 = u(149:9249).*gradCX' + v(149:9249).*gradCY';

% Residual (TOTAL DC/Dt) for Case 1.
% Label as TD for Total Derivative. Units = mg.m-3.s-1.
TD = EUL + ADV_1(2:end)';

% Case 2. Advective Change in Concentration INCLUDING velocity gradients.
% <=> c(du/dx + dv/dy) NOT assumed to be zero
% [Units of mg.m-3.s-1]
% First, find the velocity gradient. [Units of m.s-1.m-1 = s-1].
% OSCAR resolution = 0.25 deg. These four points surround Stn ALOHA => they
% are all 0.5 degrees latitude/longitude apart. 
dxx = 51.27*1000; dyy = 55.60*1000; % distance [m].
gradU = (uR(149:9249,2) - uR(149:9249,4))./dxx;
gradV = (vR(149:9249,1) - vR(149:9249,3))./dyy;
% Calculate total advective component. Grad.(u.c).
% g2 = u(149:9249).*gradCX' + v(149:9249).*gradCY' + CHL_ALOHA'.*gradU' + CHL_ALOHA'.*gradV';
g2 = ADV_1 + CHL_study'.*gradU' + CHL_study'.*gradV';
T1 = u(149:9249).*gradCX';
T2 = v(149:9249).*gradCY';
T3 = CHL_study'.*gradU';
T4 = CHL_study'.*gradV';
% Residual for Case 2.
dcdt_2 = EUL - g2(2:end)';

% % Case 3. Calculate grad(cu) in one without using chain rule.
% % Find CHL at each of four points.... [mg.m-3]
% % Replace this analysis with gridded interpolant. And also expand.
% CHLr(:,1) = squeeze(0.25*(CHL(50,50,A:B)+CHL(51,50,A:B)+CHL(51,49,A:B)+CHL(50,49,A:B)));
% CHLr(:,2) = squeeze(0.25*(CHL(51,49,A:B)+CHL(52,49,A:B)+CHL(52,48,A:B)+CHL(51,48,A:B)));
% CHLr(:,3) = squeeze(0.25*(CHL(50,48,A:B)+CHL(51,48,A:B)+CHL(51,47,A:B)+CHL(50,47,A:B)));
% CHLr(:,4) = squeeze(0.25*(CHL(49,49,A:B)+CHL(50,49,A:B)+CHL(50,48,A:B)+CHL(49,48,A:B)));
% % dxx = 51.27*1000; dyy = 55.60*1000; % distance [m].
% % ... then find grad(uc). [mg.m-3.s-1]
% % g3 = (uR(149:9249,2).*CHLr(:,2) - uR(149:9249,4).*CHLr(:,4))./dxx + (vR(149:9249,1).*CHLr(:,1) - vR(149:9249,3).*CHLr(:,3))./dyy;
% g3_1 = (uR(149:9249,2).*CHLr(:,2) - uR(149:9249,4).*CHLr(:,4))./dxx;
% % g3_1 = (uR(149:9249,4).*CHLr(:,4) - uR(149:9249,2).*CHLr(:,2))./dxx;
% g3_2 = (vR(149:9249,1).*CHLr(:,1) - vR(149:9249,3).*CHLr(:,3))./dyy;
% % g3_2 = (vR(149:9249,3).*CHLr(:,3) - vR(149:9249,1).*CHLr(:,1))./dyy;
% g3 = g3_1 + g3_2;

% % TEMPORARY fix. Case 3 is too small. 
% fix = 1; % 5.5 works to make it similar to Cases 1 and 2 (???)
% g3 = (fix)*g3;

% % Residual for Case 3.
% dcdt_3 = EUL - g3(2:end);

% Fill missing values.
Af = fillmissing(ADV_1,"nearest");
g2f = fillmissing(g2,"nearest");
% g3f = fillmissing(g3,"nearest");
Ef = fillmissing(EUL,"nearest");
dcdt1f = fillmissing(TD,"nearest");
dcdt2f = fillmissing(dcdt_2,"nearest");
% dcdt3f = fillmissing(dcdt_3,"nearest");
CHLf = fillmissing(CHL_study,"nearest");

% if applyCorrelations == true
%     % CALCULATE the correlation of DC/Dt with each of the three advection
%     % estimates.
% 
%     rho_ = [corr(Ef,Af(2:end)') corr(Ef,g2f(2:end)') corr(Ef,g3f(2:end))];% corr(g1f(2:end),g2f(2:end)) corr(g1f(2:end),g3f(2:end)') corr(g2f(2:end),g3f(2:end)')] ;
%     rho__ = [corr(Af',g2f') corr(Af',g3f) corr(g2f',g3f)] ;
%     % Actually I need correlation of advection and local component.
%     rho2 = [corr(dcdt1f,Af(2:end)') corr(dcdt2f,g2f(2:end)')];
% end
%% 40-Day approach.

dCdt_40 = movmean(EUL,40,1,"omitnan");
g1_40 = movmean(ADV_1,40,2,"omitnan");

figure
plot(time(A+1:B),EUL); hold on
plot(time(A+1:B),dCdt_40); hold off

dcdt_40 = dCdt_40 - g1_40(2:end)';

if applyCorrelations == true, rho40 = corr(dcdt_40,g1_40(2:end)'); end

means_40m = [mean(abs(dCdt_40),1,"omitnan") mean(abs(dcdt_40),1,"omitnan") mean(abs(g1_40),2,"omitnan")];

%% Show figures incl. additional statistics calculations.

% SHOW surface CHL and DC/DT at Stn ALOHA. [mg.m-3]
figure; subplot(2,1,1); plot(time(A:B),CHL_study); ylabel("CHL [mg m$^{-3}$]","FontSize",12,Interpreter="latex");
title("Surface CHL at Stn ALOHA","FontSize",12,Interpreter="latex");
subplot(2,1,2); plot(time(A+1:B),EUL); ylabel("$\frac{D [CHL]}{Dt}$ [mg m$^{-3}$ day$^{-1}]$",Interpreter="latex",FontSize=12);
title("CHL Change over Time","FontSize",12,Interpreter="latex");

% SHOW Advective Change over Time. (Case One estimate)
figure;
plot(time(A:B),ADV_1,Color="b"); hold on
plot(time(A:B),Af,"LineStyle",":",Color="b"); 
hold off
ylabel("Advection [mg m$^{-3}$ s$^{-1}$]","FontSize",12,Interpreter="latex");
title(["Advection over Stn ALOHA","Case 1: $\vec{u}.\nabla c = u \frac{\partial c}{\partial x} + v \frac{\partial c}{\partial y}$ "],"FontSize",12,Interpreter="latex");

if compareCaseEsts == true
    % COMPARE different estimates of advection.
    % 1. Initial estimate. Advection = vec(u).GRAD(c) = u.dc/dx + v.dc/dy.
    % 2. Including gradient in velocity using chain rule. Advection = 
    % GRAD(u.c) = u.dc/dx + v.dc/dy + c(x).du/dx + c(y).dv/dy.
    % 3. Including gradient in velocity without chain rule. Advection =
    % (u2c4 - u1c3)/dx + (v2c2 - v1c1)/dy.
    figure; 
    subplot(2,1,1)
    plot(time(A:B),ADV_1,"LineWidth",1.5,"Color","k",DisplayName="Case 1"); hold on
    plot(time(A:B),Af,"LineWidth",1,"LineStyle",":","Color","k",HandleVisibility="off"); 
    plot(time(A:B),g2,"LineWidth",1.5,"Color","b",DisplayName="Case 2"); 
    plot(time(A:B),g2f,"LineWidth",1,"LineStyle",":","Color","b",HandleVisibility="off");
    semilogy(time(A:B),g3',"LineWidth",1.5,"Color","r",DisplayName="Case 3");
    plot(time(A:B),g3f,"LineWidth",1,"LineStyle",":","Color","r",HandleVisibility="off"); 
    hold off
    ylabel("Advection [mg m$^{-3}$ s$^{-1}$]",Interpreter="latex"); 
    title("Comparison of Advection Estimates"); legend();
    subplot(2,1,2)
    semilogy([1 2 3],[mean(abs(ADV_1),2,"omitnan") mean(abs(g2),2,"omitnan") mean(abs(g3),1,"omitnan")],"o-",DisplayName="Cases"); 
    title("Scale Analysis of Advection Estimates"); legend(); xticks([1 2 3]); xticklabels([1 2 3]);
end

if fourTermComp == true
    % SHOW Relative Importance of Each Term in Advection.
    figure;
    subplot(2,1,1)
    plot(time(A:B),T1); hold on; plot(time(A:B),T2);
    plot(time(A:B),T3); plot(time(A:B),T4); hold off
    legend("u.dc/dx","v.dc/dy","c.du/dx","c.dv/dy"); ylabel("Advection [mg m$^{-3}$ s$^{-1}$]",Interpreter="latex")
    title("Case 2: Four Terms Contributing Towards Advection");
    subplot(2,1,2)
    semilogy([1 2 3 4],[abs(mean(T1,2,"omitnan")) abs(mean(T2,2,"omitnan")) abs(mean(T3,2,"omitnan")) abs(mean(T4,2,"omitnan"))],'o:');% hold on; semilogy(2,abs(mean(T2,2,"omitnan"))); semilogy(3,abs(mean(T3,2,"omitnan"))); semilogy(4,abs(mean(T4,2,"omitnan"))); hold off;
    xticks([1 2 3 4]); ylabel("Mean Advection [mg m$^{-3}$ s$^{-1}$]",Interpreter="latex");
    xticklabels(["u.dc/dx" "v.dc/dy" "c.du/dx" "c.dv/dy"]);
    title("Scale Comparison of Four Terms");
end

%% Calculate Means (Scale Analysis).
% Comparison of the means and medians can tell us how comparable each
% dataset is.

% Convert E from day-1 to s-1. Already done.
% EUL = EUL;

% Advection, Eulerian
means = [mean(abs(ADV_1),2,"omitnan") mean(abs(EUL),1,"omitnan") mean(abs(TD),1,"omitnan")]; % mean(abs(g2),2,"omitnan") mean(abs(g3),1,"omitnan") mean(abs(dcdt_1),1,"omitnan") mean(abs(dcdt_2),1,"omitnan") ];
medians = [median(abs(ADV_1),2,"omitnan") median(abs(EUL),1,"omitnan") median(abs(TD),1,"omitnan")]; % median(abs(g2),2,"omitnan") median(abs(g3),1,"omitnan")];
% modes = [mode(abs(E),1)/86400 mode(abs(g1),2) mode(abs(g2),2) mode(abs(g3),1)];
% g3sss = [mean(abs(g3),1,"omitnan") mean(abs(g3_1),1,"omitnan") mean(abs(g3_2),1,"omitnan")];

stds = [std(abs(ADV_1),"omitnan") std(abs(EUL),"omitnan") std(abs(TD),"omitnan")];
% 40 Days done in previous section.


% 1/3 Year. 120 Days.
EUL_120 = movmean(EUL,120,1,"omitnan");
g1_120 = movmean(ADV_1,120,2,"omitnan");
% figure
% plot(time(A+1:B),EUL); hold on
% plot(time(A+1:B),DCDt_120); hold off
means_120 = [mean(abs(EUL_120),1,"omitnan") mean(abs(g1_120),2,"omitnan")];
std_120 = [std(abs(EUL_120),"omitnan") std(abs(g1_120),"omitnan")];

% 1/2 Year. 180 Days.
EUL_180 = movmean(EUL,180,1,"omitnan");
g1_180 = movmean(ADV_1,180,2,"omitnan");
means_180 = [mean(abs(EUL_180),1,"omitnan") mean(abs(g1_180),2,"omitnan")];
std_180 = [std(abs(EUL_180),"omitnan") std(abs(g1_180),"omitnan")];

% 1 Year. 365.25 Days.
EUL_365 = movmean(EUL,365.25,1,"omitnan");
g1_365 = movmean(ADV_1,365.25,2,"omitnan");
means_365 = [mean(abs(EUL_365),1,"omitnan") mean(abs(g1_365),2,"omitnan")];
std_365 = [std(abs(EUL_365),"omitnan") mean(std(g1_365),"omitnan")];

% 23 Days.
EUL_23 = movmean(EUL,23,1,"omitnan");
g1_23 = movmean(ADV_1,23,2,"omitnan");
means_23 = [mean(abs(EUL_23),1,"omitnan") mean(abs(g1_23),2,"omitnan")];
std_23 = [std(abs(EUL_23),"omitnan") std(abs(g1_23),"omitnan")];

% 12 Days.
EUL_12 = movmean(EUL,12,1,"omitnan");
g1_12 = movmean(ADV_1,12,2,"omitnan");
means_12 = [mean(abs(EUL_12),1,"omitnan") mean(abs(g1_12),2,"omitnan")];
std_12 = [std(abs(EUL_12),"omitnan") std(abs(g1_12),"omitnan")];

% 2 Days
EUL_2 = movmean(EUL,2,1,"omitnan");
g1_2 = movmean(ADV_1,2,2,"omitnan");
means_2 = [mean(abs(EUL_2),1,"omitnan") mean(abs(g1_2),2,"omitnan")];
std_2 = [std(abs(EUL_2),"omitnan") std(abs(g1_2),"omitnan")];


% TEST the distribution of 'advCon' in three steps...
% 1. Distribution.
% This can't be lognormal because of the negatives => test vs normal instead.
edges = -3e-7:2.5e-8:3e-7;

fig = figure; subplot(3,1,1);
if histfitOrGram == true, histfit(ADV_1,100,"normal"); legend("Data","Normal");
else, histogram(ADV_1,BinEdges=edges); legend("Data",Location="northwest"); end
title("Advected Concentration","Fontsize",10,Interpreter="latex");
%xlabel("Concentration (mg m$^{-3}$ s$^{-1}$)",Interpreter="latex");
xlim([-5e-7 5e-7]); ylim([0 3000]); xticklabels({});
subplot(3,1,2);
if histfitOrGram == true, histfit(EUL,100,"normal"); else,...
    histogram(EUL,BinEdges=edges); end
title("Eulerian/Local Concentration Change (dc/dt)","Fontsize",10,Interpreter="latex");
%xlabel("Concentration (mg m$^{-3}$ s$^{-1}$)",Interpreter="latex");
xlim([-5e-7 5e-7]); ylim([0 3000]); xticklabels({});
subplot(3,1,3); 
if histfitOrGram == true, histfit(TD,100,"normal"); else,...
        histogram(TD,BinEdges=edges); end
title("Total Concentration Change (Dc/Dt)","Fontsize",10,Interpreter="latex");
xlabel("Concentration (mg m$^{-3}$ s$^{-1}$)","Fontsize",8,Interpreter="latex");
xlim([-5e-7 5e-7]); ylim([0 3000]);
% ...
% 2. A-D Test.
[h_adv,p_adv] = adtest(ADV_1,"Distribution","norm");
if h_adv == 1, advN = "Not normal (p = "+num2str(p_adv)+")"; end
[h_dc,p_dc] = adtest(TD,Distribution="norm");
if h_adv == 1, dcN = "Not normal (p = "+num2str(p_adv)+")"; end
[h_Dc,p_Dc] = adtest(EUL,Distribution="norm");
if h_adv == 1, DcN = "Not normal (p = "+num2str(p_adv)+")"; end
% 3. Skewness-Kurtosis.
sk = [skewness(ADV_1) skewness(EUL) skewness(TD)];
ku = [kurtosis(ADV_1) kurtosis(EUL) kurtosis(TD)];
% ...
% Plot above results to same figure.
n(1,:) = {"Mean = "+means(1),"STD = "+stds(1)," Kurtosis = "+num2str(ku(1)), "Skewness = "+num2str(sk(1)), advN};
n(2,:) = {"Mean = "+means(2),"STD = "+stds(2),"Kurtosis = "+num2str(ku(2)), "Skewness = "+num2str(sk(2)), DcN};
n(3,:) = {"Mean = "+means(3),"STD = "+stds(3),"Kurtosis = "+num2str(ku(3)), "Skewness = "+num2str(sk(3)), dcN};
annotation('textbox',[.65 .77 .25 .16],'String',n(1,:),'FontSize',7); 
annotation('textbox',[.65 .43 .25 .16],'String',n(2,:),'FontSize',7);
annotation('textbox',[.65 .12 .25 .16],'String',n(3,:),'FontSize',7);
hold off
% ylabel("No. of counts","FontSize",12,Interpreter="latex");
% xlabel("Advection [mg m$^{-3}$ s$^{-1}$]","FontSize",12,Interpreter="latex")
% title("Distribution of Advected Concentration","FontSize",12,Interpreter="latex");

han=axes(fig,'visible','off'); 
han.YLabel.Visible='on';
yl = ylabel(han,'No. of Counts');
yl.Position(2) = 0.55;
yl.Position(1) = -0.12;


if compareTermsInG3 == true
    % COMPARE terms within case "g3". Trying to find a reason why g3 is so
    % anomalously large.
    figure;
    plot(time(A:B),g3); hold on
    plot(time(A:B),g3_1); plot(time(A:B),g3_2); hold off
    legend("Case 3","Term 1","Term 2");
    disp([mean(g3,1,"omitnan") mean(g3_1,1,"omitnan") mean(g3_2,1,"omitnan")]);
    title("Comparison of Terms in Case 3");
end

% Skewness-Kurtosis Plot.
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
skLognN = -skLogn;
kuLognN = kuLogn;

figure;
scatter(0,3,72,[0 0 0],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(sk(1),ku(1),"filled",DisplayName="Advected Concentration");
scatter(sk(2),ku(2),"filled",DisplayName="Local Change");
scatter(sk(3),ku(3),"filled",DisplayName="Total Change");
hold off;
xlabel("Skewness",Interpreter="latex"); ylabel("Kurtosis",Interpreter="latex"); 
ylim([0 40]); xlim([-0.5 2]);
legend(); title("Skewness vs. Kurtosis",Interpreter="latex");




%% FFT of CHL @Stn ALOHA.
if showFft == true

    % % --- CHL --- % %

    % Freq = day.-1
    CHLf = fillmissing(CHL_study,"nearest",1);
    
    Y = fft(CHLf);
    Fs = 1; % 1 per day
    L = length(CHL_study);
    
    smthi = movmean(abs(Y),5);
    
    figure;
    loglog(Fs/L*(0:L-1),abs(Y),"LineWidth",2,DisplayName="FFT"); hold on
    xline(0.0027,'-','1 Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0055,'-','1/2-Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0082,'-','1/3-Year','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    xline(0.0435,'-','23 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    loglog(Fs/L*(26:L-1),smthi(27:end),"LineWidth",2,DisplayName="Smoothed FFT");
    hold off;
    xlabel("Frequency (CPD)",Interpreter="latex"); ylabel("$\|$ FFT(CHL) $\|$",Interpreter="latex");
    legend("FFT","Mean FFT");
    title("FFT: CHL @Stn ALOHA (1997-2022)",Interpreter="latex");

    % % --- u --- % %
    
    uf = fillmissing(u,"nearest");
    Y1a = fft(uf);
    Fs = 1; % 1 per day
    L1a = length(uf);
    
    smthi = movmean(abs(Y1a),5);

    figure;
    loglog(Fs/L1a*(0:L1a-1),abs(Y1a),"LineWidth",2,DisplayName="FFT"); hold on
    xline(0.0055,'-','1/2 Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0119,'-','84 Days','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0278,'-','36 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    xline(0.0833,'-','12 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    loglog(Fs/L1a*(26:L1a-1),smthi(27:end),"LineWidth",2,DisplayName="Smoothed FFT");
    hold off;
    xlabel("Frequency (CPD)",Interpreter="latex"); ylabel("$\|$ FFT(u) $\|$",Interpreter="latex");
    legend("FFT","Mean FFT",Location="best");
    title("FFT: u @Stn ALOHA (1997-2002)",Interpreter="latex");

     % % --- v --- % %
    
    vf = fillmissing(v,"nearest");
    Y1b = fft(vf);
    Fs = 1; % 1 per day
    L1b = length(vf);
    
    smthi = movmean(abs(Y1b),5);

    figure;
    loglog(Fs/L1b*(0:L1b-1),abs(Y1b),"LineWidth",2,DisplayName="FFT"); hold on
    xline(0.0055,'-','1/2 Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0119,'-','84 Days','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0278,'-','36 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    xline(0.0833,'-','12 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    loglog(Fs/L1b*(26:L1b-1),smthi(27:end),"LineWidth",2,DisplayName="Smoothed FFT");
    hold off;
    xlabel("Frequency (CPD)",Interpreter="latex"); ylabel("$\|$ FFT(v) $\|$",Interpreter="latex");
    legend("FFT","Mean FFT",Location="best");
    title("FFT: v @Stn ALOHA (1997-2002)",Interpreter="latex");

    % % --- DC/Dt --- % %
    
    Y2 = fft(Ef);
    Fs = 1; % 1 per day
    L2 = length(Ef);
    
    smthi = movmean(abs(Y2),5);

    figure;
    loglog(Fs/L2*(0:L2-1),abs(Y2),"LineWidth",2,DisplayName="FFT"); hold on
    xline(0.0027,'-','1 Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0055,'-','1/2-Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0082,'-','1/3-Year','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    xline(0.0435,'-','23 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    xline(0.5000,'-','2 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    loglog(Fs/L2*(26:L2-1),smthi(27:end),"LineWidth",2,DisplayName="Smoothed FFT");
    hold off;
    xlabel("Frequency (CPD)"); ylabel("|fft(DC/Dt)|");
    legend("FFT","Mean FFT",Location="best");
    title("FFT: D[CHL]/Dt @Stn ALOHA (1997-2002)");

    % % --- u.GRAD(C) --- % %
    
    Y3 = fft(Af);
    Fs = 1; % 1 per day
    L3 = length(Af);
    
    smthi = movmean(abs(Y3),5);

    figure;
    loglog(Fs/L3*(0:L3-1),abs(Y3),"LineWidth",2,DisplayName="FFT"); hold on
    xline(0.0027,'-','1 Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0055,'-','1/2-Year','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center',HandleVisibility='off');
    xline(0.0119,'-','84 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    xline(0.0357,'-','28 Days','LabelHorizontalAlignment','center',LabelVerticalAlignment='bottom',HandleVisibility='off'); 
    loglog(Fs/L3*(26:L3-1),smthi(27:end),"LineWidth",2,DisplayName="Smoothed FFT");
    hold off;
    xlabel("Frequency (CPD)",Interpreter="latex"); ylabel("$\|$ FFT(u.$\nabla{c}$) $\|$",Interpreter="latex");
    legend("FFT","Mean FFT",Location="best");
    title("FFT: u.$\nabla{c}$ @Stn ALOHA (1997-2002)",Interpreter="latex");
end

%% Integrate terms.

% Time. [Units = s]
X = (1:1:9100)*86400;
X2 = (0:1:9100)*86400;

% Window length. [Units = s]
Window(1) = X(end)-X(1);
Window(2) = X2(end)-X2(1);

% Eulerian. [Units = mg.m-3.s-1]
tmpE = EUL;
tmpE(isnan(EUL)) = 0;

% Advective Change in Concentration [Units of mg.m-3.s-1].
tmpA = ADV_1;
tmpA(isnan(ADV_1)) = 0;

% Total derivative. [Units of mg.m-3.s-1].
tmpTD = TD;
tmpTD(isnan(TD)) = 0;

tmpCHL = CHL_study; tmpCHL(isnan(CHL_study)) = 0;

% Calculate integral of Eulerian term, Advective term, and Total Derivative
% between X_end and X_0.

% (this should equal CHL_ALOHA). Mean CHL_ALOHA is O(-2). Integral of 
% CHL_ALOHA is also O(-2). Integral of Eulerian is O(-11).
E_INT = trapz(X,abs(tmpE)); %/Window(1); 
tmp1 = cumsum(abs(tmpE),1,"omitnan");
E_INT2 = tmp1(9100)*86400;

CHL_INT = trapz(X2,tmpCHL)/Window(2);

A_INT = trapz(X2,abs(tmpA)); %/Window(2);
% tmp2 = cumsum(abs(ADV_1),2,"omitnan");
tmp2 = 86400*cumsum(ADV_1,2,"omitnan");
A_INT2 = tmp2(9101);

TD_INT = trapz(X,abs(tmpTD)); %/Window(1);
tmp3 = cumsum(abs(tmpTD),1,"omitnan");
TD_INT2 = tmp3(9100) * 86400;

ADV_timeseries = tmp2 + abs(min(tmp2));

figure;
subplot(1,4,[1 2])
plot(t(A:B),ADV_timeseries);
ylabel("Concentration (mg m^{-3})");
xlabel("Time"); title("Constructed Time Series");
subplot(1,4,3)
histogram(ADV_timeseries); xlabel("Concentration"); ylabel("No. of Counts");
title("Histogram");

subplot(1,4,4)
% Skewness-Kurtosis Plot.
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end
skLognN = -skLogn;
kuLognN = kuLogn;
scatter(0,3,72,[0 0 0],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(skewness(ADV_timeseries),kurtosis(ADV_timeseries),"filled",DisplayName="Advected Concentration");
% scatter(sk(2),ku(2),"filled",DisplayName="Local Change");
% scatter(sk(3),ku(3),"filled",DisplayName="Total Change");
hold off;
xlabel("Skewness",Interpreter="latex"); ylabel("Kurtosis",Interpreter="latex"); 
ylim([2 10]); xlim([-0.5 1.5]);
legend(); title("Skewness vs. Kurtosis",Interpreter="latex");
sgtitle("Advected Concentration");



% Actual time series.
figure; subplot(1,4,[1 2]);
plot(time(A:B),CHL_study); ylabel("CHL [mg m$^{-3}$]","FontSize",12,Interpreter="latex");
title("Surface CHL at Stn ALOHA","FontSize",12,Interpreter="latex");
subplot(1,4,3)
histogram(CHL_study); xlabel("CHL (mg m^{-3})"); ylabel("No. of Counts");
subplot(1,4,4)
% Skewness-Kurtosis Plot.
kuLognN = kuLogn;
scatter(0,3,72,[0 0 0],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(skewness(CHL_study),kurtosis(CHL_study),"filled",DisplayName="Eulerian Concentration");
% scatter(sk(2),ku(2),"filled",DisplayName="Local Change");
% scatter(sk(3),ku(3),"filled",DisplayName="Total Change");
hold off;
xlabel("Skewness",Interpreter="latex"); ylabel("Kurtosis",Interpreter="latex"); 
ylim([2 15]); xlim([-0.5 2]);
legend(); title("Skewness vs. Kurtosis",Interpreter="latex");
sgtitle("Eulerian Concentration");


% Lagrangian estimate.
LAG = CHL_study + ADV_timeseries';
figure;
subplot(1,4,[1 2]);
plot(time(A:B),LAG); ylabel("CHL [mg m$^{-3}$]","FontSize",12,Interpreter="latex");
title("Lagrangian CHL at Stn ALOHA","FontSize",12,Interpreter="latex");
subplot(1,4,3)
histogram(LAG); xlabel("CHL (mg m^{-3})"); ylabel("No. of Counts");
subplot(1,4,4)
% Skewness-Kurtosis Plot.
kuLognN = kuLogn;
scatter(0,3,72,[0 0 0],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(skewness(LAG),kurtosis(LAG),"filled",DisplayName="Lagrangian Concentration");
% scatter(sk(2),ku(2),"filled",DisplayName="Local Change");
% scatter(sk(3),ku(3),"filled",DisplayName="Total Change");
hold off;
xlabel("Skewness",Interpreter="latex"); ylabel("Kurtosis",Interpreter="latex"); 
ylim([2 15]); xlim([-0.5 2]);
legend(); title("Skewness vs. Kurtosis",Interpreter="latex");
sgtitle("Lagrangian Concentration");

figure;
subplot(1,4,[1 2]);
plot(time(A:B),ADV_timeseries); hold on
plot(time(A:B),CHL_study); 
plot(time(A:B),LAG); hold off
ylabel("CHL (mg m$^{-3}$)","FontSize",12,Interpreter="latex");
title("CHL at Stn ALOHA","FontSize",12,Interpreter="latex");
subplot(1,4,4)
% Skewness-Kurtosis Plot.
kuLognN = kuLogn;
scatter(0,3,72,[0 0 0],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
scatter(skewness(ADV_timeseries),kurtosis(ADV_timeseries),"filled",DisplayName="Advected");
scatter(skewness(CHL_study),kurtosis(CHL_study),"filled",DisplayName="Eulerian");
scatter(skewness(LAG),kurtosis(LAG),"filled",DisplayName="Lagrangian");
hold off;
xlabel("Skewness",Interpreter="latex"); ylabel("Kurtosis",Interpreter="latex"); 
ylim([2 15]); xlim([-0.2 2.5]);
legend(); title("Skewness vs. Kurtosis",Interpreter="latex");
sgtitle("Components of Concentration");


%% Correlate CHL and derivative.

TDf = fillmissing(TD,"nearest");
% TDf2 = TDf(2:end)
tmpA = corr(TDf,Ef);

tmpB = corr(CHLf(2:end),Ef);
tmpC = corr(CHLf,Af');
tmpD = corr(CHLf(2:end),TDf);

%% Clear unnecessary variables
clear advN applyCorrelations compareCaseEsts compareTermsInG3 dcN DcN ...
    fourTermComp i input3 make_it_tight readOscar regio showFft;