% Statistical Analysis of upper 200 dbar at Station ALOHA for chl-a, other pigments, and BGC variables.
% We import the data and run a hypothesis test on it with the respective
% null hypotheses of normal and lognormal. We use the Anderson-Darling test
% since this is both more powerful than similar tests such as
% Kolmogorov-Smirnov and more flexible than tests such as Shapiro-Wilks
% which did not easily allow for testing of other distributions.

clear; clc; close all;
addpath("func\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 15 15]);

% Options and test cases.
thresh = 30;
principle = true;       % main analysis
seasonal = false;       % seasonality of statistics: A-D
logAxes = true;         % output p-values as log values
showL0title = true;     % show title for chl-a plot. (off for paper)

if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Chl-a depth and time-series.
% Import Chl-a data and transform the data.
% Plot the depth and time-series of chla across the upper 200 dbar (approx
% 200 m) between 1988 and 2021.

% Import data file
tmp = importdata('data/hot_chla.txt');

% Import binned pressure and concentration of chl-a
[~,~,~,pB,chla,~] = L0_helper(tmp,thresh,'ad',true,0,true);

time = tmp.data(:,2);
hms = tmp.data(:,3);

YY1 = mod(time,100);

% Correct the years
for i=1:length(YY1)
    if YY1(i) < 85
        YY1(i) = YY1(i) + 2000;
    else
        YY1(i) = YY1(i) + 1900;
    end
end

YY2 = string(compose('%02d',YY1));
YY = str2double(YY2);

DD1 = mod(time,10000)-mod(time,100);
DD2 = DD1/100;
DD3 = string(compose('%02d',DD2));
DD = str2double(DD3);

MM1 = mod(time,1000000) - mod(time,10000);
MM2 = MM1/10000;
MM3 = string(compose('%02d',MM2));
MM = str2double(MM3);

hms(hms==-9) = nan;

ss1 = mod(hms,100);
ss2 = string(compose('%02d',ss1));
ss = str2double(ss2);

mm1 = mod(hms,10000)-mod(hms,100);
mm2 = mm1/100;
mm3 = string(compose('%02d',mm2));
mm = str2double(mm3);

hh1 = mod(hms,1000000) - mod(hms,10000);
hh2 = hh1/10000;
hh3 = string(compose('%02d',hh2));
hh = str2double(hh3);


time2 = datetime(YY,MM,DD,hh,mm,ss);

newCast = [1];

for i = 2:length(pB)
    if pB(i) < pB(i-1)
        disp(i);
        newCast = [newCast i];
    end
end

pgrid = nan(320,20);
tgrid = NaT(320,20);
chlagrid = NaN(320,20);

% First Case.
pgrid(1,1:(newCast(2)-newCast(1))) = pB(newCast(1):newCast(2)-1);
tgrid(1,1:(newCast(2)-newCast(1))) = time2(newCast(1):newCast(2)-1);
chlagrid(1,1:(newCast(2)-newCast(1))) = chla(newCast(1):newCast(2)-1);

% Then Loop.
for i = 2:319
    pgrid(i,1:(newCast(i+1)-newCast(i))) = pB(newCast(i):newCast(i+1)-1);
    tgrid(i,1:(newCast(i+1)-newCast(i))) = time2(newCast(i):newCast(i+1)-1);
    chlagrid(i,1:(newCast(i+1)-newCast(i))) = chla(newCast(i):newCast(i+1)-1);
end

% Final Case.
pgrid(320,1:(length(pB)+1-newCast(320))) = pB(newCast(320):length(pB));
tgrid(320,1:(length(time2)+1-newCast(320))) = time2(newCast(320):length(time2));
chlagrid(320,1:(length(chla)+1-newCast(320))) = chla(newCast(320):length(chla));

tgridDatenum = datenum(tgrid);

figure;
contourf(tgridDatenum,pgrid,chlagrid,'LineColor','auto');
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

%% Chl-a mean profile.

chlaProfile = nan(1,20);
f5 = nan(1,20);
f95 = nan(1,20);

for i = 1:20
    chlaProfile(i) = mean(chlagrid(pgrid==i),"omitnan");
    f5(i) = prctile(chlagrid(pgrid==i),5);
    f95(i) = prctile(chlagrid(pgrid==i),95);
end

% Toggle show y-label and title (for paper we don't need either)
displayYLabelAndTitle = false;

% Plot the mean profile of fluorescence with the mean and confidence
% interval.
figure;
plot(chlagrid(:,1),pgrid(:,1),'.',Color=[0.8 0.8 0.8],DisplayName="raw data");
hold on
plot(chlagrid(:,2:20),pgrid(:,2:20),'.',Color=[0.8 0.8 0.8],HandleVisibility='off');
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
ylim([1 18]);
ax = gca;
ax.FontSize = 15;


%% Chl-a. Parameters vs Depth (assuming a lognormal distribution).

XN = nan(20,500);
for i = 1:20
    loc = find(pB == i);
    XN(i,1:length(loc)) = chla(loc);
end

figure;
histogram(XN(15,:));
XN(XN<=0) = nan;

for i = 1:20
    phat(i,:) = mle(XN(i,:),distribution="Lognormal");
end

depths = 5:10:195;

figure;
subplot(1,2,1)
plot(exp(phat(:,1)),depths);
title("$\mu^*$",Interpreter="latex"); set(gca,"YDir","reverse"); grid on
subplot(1,2,2)
plot(exp(phat(:,2)),depths);
title("$\sigma^*$",Interpreter="latex"); set(gca,"YDir","reverse"); grid on
sgtitle("L0 chl-$a$ parameters (1988-2021)",Interpreter="latex");

%% A-D Test.
if principle == true    
    
    % A-D
    tmpT = "-ad";
    
    % chla
    tmpX ="";
    tmp = importdata('data/hot_chla.txt');
    L0_helper(tmp,thresh,'ad');
    if showL0title == true
        sgtitle("L0 chl-$a$"+tmpX,"Interpreter","latex");
    end

    % chlb
    tmp = importdata('data/hot_chlb.txt');
    L0_helper(tmp,thresh,'ad');
    sgtitle("L0: Chl $b$","Interpreter","latex");
   
    % pc
    tmp = importdata('data/hot_pc.txt');
    L0_helper(tmp,thresh,'ad');
    sgtitle("L0: Particulate Carbon","Interpreter","latex");
 
end
%% A-D Test: Seasonal.
if seasonal == true

    % WINTER
    tmpT = "-ad-01"; tmpX = ": L0 Winter";

    % chla
    tmp = importdata('data/hot_chla.txt');
    L0_helper(tmp,thresh,"ad",logAxes,1);
    sgtitle("Chl $a$ (1988-2021)"+tmpX,"Interpreter","latex");

    % pc
    tmp = importdata('data/hot_pc.txt');
    L0_helper(tmp,thresh,"ad",logAxes,1);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");


    % SPRING
    tmpX = ": L0 Spring";

    % chla
    tmp = importdata('data/hot_chla.txt');
    L0_helper(tmp,thresh,"ad",logAxes,2);
    sgtitle("Chl-$a$"+tmpX,"Interpreter","latex");

    % pc
    tmp = importdata('data/hot_pc.txt');
    L0_helper(tmp,thresh,"ad",logAxes,2);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");


    % SUMMER
    tmpX = ": L0 Summer";
    
    % chla
    tmp = importdata('data/hot_chla.txt');
    ax = L0_helper(tmp,thresh,"ad",logAxes,3);
    sgtitle("Chl-$a$"+tmpX,"Interpreter","latex");

    % pc
    tmp = importdata('data/hot_pc.txt');
    L0_helper(tmp,thresh,"ad",logAxes,3);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");


    % AUTUMN
    tmpX = ": L0 Autumn";
    
    % chla
    tmp = importdata('data/hot_chla.txt');
    L0_helper(tmp,thresh,"ad",logAxes,4);
    sgtitle("Chl-$a$"+tmpX,"Interpreter","latex");

    % pc
    tmp = importdata('data/hot_pc.txt');
    L0_helper(tmp,thresh,"ad",logAxes,4);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
end