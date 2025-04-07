function [R,p2] = bbvuong(x)
%
%
%   function [R,p2] = bbvuong(x)
%
% Calculates Vuong's test Log-Likelihood Ratio (R) and p values (only two
%   tails) for some distribution:
%   Normal, Log-Normal ,Weibull, Gamma and Exponential.
%


MLEp = NaN*ones(5,2);
[MLEp(1,:),c95] = mle(x,'distribution','norm');
if sum(x<=0) == 0
    [MLEp(2,:),c95] = mle(x,'distribution','logn');
    [MLEp(3,:),c95] = mle(x,'distribution','wbl');
    [MLEp(4,:),c95] = mle(x,'distribution','gamma');
    [MLEp(5,1),c95] = mle(x,'distribution','exp');
else MLEp(2,:) = [NaN NaN];
    display('Can''t compute MLE parameters for Lognormal, Weibull, Gamma and Exponential distribution. Negative or zero values in input vector.')
end
l_norm = log(pdf('norm',x,MLEp(1,1),MLEp(1,2)));
l_logn = log(pdf('logn',x,MLEp(2,1),MLEp(2,2)));
l_wbl = log(pdf('wbl',x,MLEp(3,1),MLEp(3,2)));
l_gam = log(pdf('gamma',x,MLEp(4,1),MLEp(4,2)));
l_exp = log(pdf('exp',x,MLEp(5,1)));
n = length(x);

% Normal vs Log-normal
l_nl = l_norm - l_logn;
R_nl = sum(l_nl);
m_nl = mean(l_nl);
% disp(n);
% disp(m_nl);
s_nl = std(l_nl);
% disp(s_nl);
v_nl = sqrt(n)*m_nl/s_nl;
% disp(v_nl);
p1_nl = cdf('norm',v_nl,0,1);
% disp(p1_nl);
if p1_nl < 0.5
    p2_nl = 2*p1_nl;
else p2_nl = 2*(1-p1_nl);
end
% disp(p2_nl);
% Normal vs Weibull
l_nw = l_norm - l_wbl;
R_nw = sum(l_nw);
m_nw = mean(l_nw);
s_nw = std(l_nw);
v_nw = sqrt(n)*m_nw/s_nw;
p1_nw = cdf('norm',v_nw,0,1);
if p1_nw < 0.5
    p2_nw = 2*p1_nw;
else p2_nw = 2*(1-p1_nw);
end
% Normal vs Gamma
l_ng = l_norm - l_gam;
R_ng = sum(l_ng);
m_ng = mean(l_ng);
s_ng = std(l_ng);
v_ng = sqrt(n)*m_ng/s_ng;
p1_ng = cdf('norm',v_ng,0,1);
if p1_ng < 0.5
    p2_ng = 2*p1_ng;
else p2_ng = 2*(1-p1_ng);
end
% Normal vs Exponential
l_ne = l_norm - l_exp;
R_ne = sum(l_ne);
m_ne = mean(l_ne);
s_ne = std(l_ne);
v_ne = sqrt(n)*m_ne/s_ne;
p1_ne = cdf('norm',v_ne,0,1);
if p1_ne < 0.5
    p2_ne = 2*p1_ne;
else p2_ne = 2*(1-p1_ne);
end
% Log-normal vs Weibull
l_lw = l_logn - l_wbl;
R_lw = sum(l_lw);
m_lw = mean(l_lw);
s_lw = std(l_lw);
v_lw = sqrt(n)*m_lw/s_lw;
p1_lw = cdf('norm',v_lw,0,1);
if p1_lw < 0.5
    p2_lw = 2*p1_lw;
else p2_lw = 2*(1-p1_lw);
end
% Log-normal vs Gamma
l_lg = l_logn - l_gam;
R_lg = sum(l_lg);
m_lg = mean(l_lg);
s_lg = std(l_lg);
v_lg = sqrt(n)*m_lg/s_lg;
p1_lg = cdf('norm',v_lg,0,1);
if p1_lg < 0.5
    p2_lg = 2*p1_lg;
else p2_lg = 2*(1-p1_lg);
end
% Log-normal vs Exponential
l_le = l_logn - l_exp;
R_le = sum(l_le);
m_le = mean(l_le);
s_le = std(l_le);
v_le = sqrt(n)*m_le/s_le;
p1_le = cdf('norm',v_le,0,1);
if p1_le < 0.5
    p2_le = 2*p1_le;
else p2_le = 2*(1-p1_le);
end
% Weibull vs Gamma
l_wg = l_wbl - l_gam;
R_wg = sum(l_wg);
m_wg = mean(l_wg);
s_wg = std(l_wg);
v_wg = sqrt(n)*m_wg/s_wg;
p1_wg = cdf('norm',v_wg,0,1);
if p1_wg < 0.5
    p2_wg = 2*p1_wg;
else p2_wg = 2*(1-p1_wg);
end
% Weibull vs Exponential
l_we = l_wbl - l_exp;
R_we = sum(l_we);
m_we = mean(l_we);
s_we = std(l_we);
v_we = sqrt(n)*m_we/s_we;
p1_we = cdf('norm',v_we,0,1);
if p1_we < 0.5
    p2_we = 2*p1_we;
else p2_we = 2*(1-p1_we);
end
% Log-normal vs Exponential
l_ge = l_gam - l_exp;
R_ge = sum(l_ge);
m_ge = mean(l_ge);
s_ge = std(l_ge);
v_ge = sqrt(n)*m_ge/s_ge;
p1_ge = cdf('norm',v_ge,0,1);
if p1_ge < 0.5
    p2_ge = 2*p1_ge;
else p2_ge = 2*(1-p1_ge);
end

rows_outputs = char('1. Normal-Log-normal','2. Normal-Weibull','3. Normal-Gamma','4. Normal-Exponential','5. Log-normal-Weibull','6. Log-normal-Gamma','7. Log-normal-Exponential','8. Weibull-Gamma','9. Weibull-Exponential','10. Gamma-Exponential');
% l = [l_nl;l_nw;l_ng;l_ne;l_lw;l_lg;l_le;l_wg;l_we;l_ge]; % EDIT PETER: compare with R
R = [R_nl;R_nw;R_ng;R_ne;R_lw;R_lg;R_le;R_wg;R_we;R_ge];
p1 = [p1_nl;p1_nw;p1_ng;p1_ne;p1_lw;p1_lg;p1_le;p1_wg;p1_we;p1_ge];
p2 = [p2_nl;p2_nw;p2_ng;p2_ne;p2_lw;p2_lg;p2_le;p2_wg;p2_we;p2_ge];