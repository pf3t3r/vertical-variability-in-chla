function Output = calculateStatistics(Input)
%CALCULATESTATISTICS Summary of this function goes here
%   Detailed explanation goes here


if ~isfield(Input,"fourDist"), Input.fourDist = false; end
if ~isfield(Input,"useVuong"), Input.useVuong = true; end
if ~isfield(Input,"calculateMoments"), Input.calculateMoments = false; end
if ~isfield(Input,"threshold"), Input.threshold = 30; end
if ~isfield(Input,"hypTest"), Input.hypTest = 'ad'; end
if ~isfield(Input,"useL2"), Input.useL2 = false; end

N = Input.N;
TH = Input.threshold;
X = Input.X;
pB = Input.pB;
hypTest = Input.hypTest;
fourDist = Input.fourDist;
useVuong = Input.useVuong;
calculateMoments = Input.calculateMoments;
if Input.useL2 == true
    pr = Input.pr;
end


ks = nan(5,N);
obs = nan(1,N);
if fourDist == true
    ad = nan(4,N);
else
    ad = nan(2,N);
end

% Calculate KS p-value, skewness, kurtosis
for i = 1:N
    % find concentration X_i at binned pressure i
    if Input.useL2 == false
        X_i = X(pB==i);
    else
        X_i = X(pB==pr(i));
    end
    % remove negative or null values
    X_i(X_i <= 0) = [];
    % apply KS test to X_i (only when at least 30 values at binned pressure)
    if length(X_i) >= TH
        %disp(i);
        %disp(length(X_i));
        if strcmp(hypTest,'ks')
            [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        else
            [~,ad(1,i)] = adtest(X_i,'Distribution','logn','alpha',0.005);
            [~,ad(2,i)] = adtest(X_i,'Distribution','norm','alpha',0.005);
            if fourDist == true
                gammaParams = mle(X_i,"distribution","Gamma");
                pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
                [~,ad(3,i)] = adtest(X_i,"Distribution","weibull");
                [~,ad(4,i)] = adtest(X_i,Distribution=pdG,MCTol=0.05);
            end
            if useVuong == true
                [rV(:,i),pV(:,i)] = bbvuong(X_i);
            end
            if calculateMoments == true
                sk(i) = skewness(X_i);
                ku(i) = kurtosis(X_i);
            end
        end
        %disp(i);
    end
    obs(i) = length(X_i);
    clear X_i;

end

ks(:,obs<TH) = nan;
ad(:,obs<TH) = nan;
if useVuong == true
    rV(:,obs<TH) = nan;
    pV(:,obs<TH) = nan;
end
if calculateMoments == true
    sk(obs<TH) = nan;
    ku(obs<TH) = nan;
end


% % Remove these nan values
% tmp = [];
% for i = 1:N
%     if ~isnan(sum(ad(:,i)))
%         tmp = [tmp i];
%     end
% end
% p = pB(tmp);
% rV = rV(:,tmp);
% pV = pV(:,tmp);
% ks = ks(:,~all(isnan(ks)));
% ad = ad(:,~all(isnan(ad)));


if Input.useL2 == false
    pXX = 5:10:195;
else
    pXX = Input.pr;
end


Output.ks = ks;
Output.ad = ad;
Output.pXX = pXX;
Output.pB = pB;
Output.obs = obs;
if useVuong == true
    Output.rV = rV;
    Output.pV = pV;
end
if calculateMoments == true
    Output.sk = sk;
    Output.ku = ku;
end


end

