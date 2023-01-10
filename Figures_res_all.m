% FADE-pqrs: Figures 5, 6, 7, 8, 9, 10
% written   by Joram Soch <Joram.Soch@DZNE.de>, 19/11/2021, 15:35 (V1)
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 28/07/2022, 16:27 (V2)
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 17/11/2022, 21:57 (V3)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/01/2023, 14:45

% clear
% close all

%%% Step 0: load the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add toolbox to path
pqrs_path = 'C:\Users\sochj\ownCloud\Tools\pqrs\pqrs\';
addpath(pqrs_path);

% load subject info
subj_file = 'subjects/subjects_FADE.xls';
[num, txt, raw] = xlsread(subj_file);
tab_hdr   = raw(1,:);
subj_info = raw(2:end,:);
clear num txt raw

% load behavioral data
bhvr_file = 'data/logfiles_FADE.mat';
load(bhvr_file);

% number of observations
n = 2*44 + 2*22;                % number of trials
N = numel(subj_ids);            % number of subjects

% define subject groups
col    =  'rb';
lab    = {'young adults', 'older adults'};
age    = cell2mat(subj_info(:,4));
ind{1} = [age <  50];
ind{2} = [age >= 50];

% define model of interest
moi      = 'pqqrrss';           % model of interest
ab_prior = [1,1];               % prior distribution
x_bins   = [0:0.01:1];          % density plot bins
p_thr    = [0.05, 0.01, 0.001]; % p-value thresholds
conf_lvl = 0.95;                % confidence level


%%% Figure 5: Bayesian model selection / family comparison %%%%%%%%%%%%%%%%

% define models and families
m = {'pqr', 'pqrs', 'pqrr', 'pqrrss', ...
     'pqqr', 'pqqrs', 'pqqrr', 'pqqrrss', ...
     'ppqr', 'ppqrs', 'ppqrr', 'ppqrrss', ...
     'ppqqr', 'ppqqrs', 'ppqqrr', 'ppqqrrss'};
f = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2;
     1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2;
     1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2;
     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2];
l = {'p', 'pp';
     'q', 'qq';
     'r', 'rr';
     'r', 'rs'};

% preallocate results
M = numel(m);
D = size(f,1);
ab_post = cell(N,M);
p_lab   = cell(1,M);
LME     = zeros(N,M);
k       = zeros(1,M);

% estimate models
fprintf('\n-> Bayesian estimation of pqrs models:\n');
for i = 1:N
    fprintf('   - Subject %s (%d out of %d): Model ', subj_ids{i}, i, N)
    x = Subj(i).rTrials( Subj(i).rTrials(:,2)==1, 4 );
    y = Subj(i).rTrials( Subj(i).rTrials(:,2)==1, 6 );
    for j = 1:M
        fprintf('%s, ', m{j});
        [ab_post{i,j}, p_lab{j}, LME(i,j), k(j)] = ME_pqrs_Bayes(y, x, m{j}, ab_prior);
    end;
    fprintf('done.\n');
end;

% calculate LFEs
LFE = zeros(N,2,D);
for d = 1:D
    LFE(:,1,d) = ME_MF_LFE(LME(:,f(d,:)==1)')';
    LFE(:,2,d) = ME_MF_LFE(LME(:,f(d,:)==2)')';
end;

% calculate PPs
PP_fam = zeros(N,2,D);
EP_fam = zeros(numel(ind),D);
for d = 1:D
    LBF = LFE(:,1,d)-LFE(:,2,d);
    PP_fam(:,1,d) = exp(LBF)./(exp(LBF)+1);
    PP_fam(:,2,d) = 1-PP_fam(:,1,d);
    for h = 1:numel(ind)
        EP_fam(h,d) = mean( PP_fam(ind{h},2,d) > PP_fam(ind{h},1,d) );
    end;
    l{d,3} = sprintf('%s vs. %s', l{d,2}, l{d,1});
end;

% estimate kernel densities
PP_med = zeros(2,D);
PP_ksd = cell(2,D);
for d = 1:D
    pp1 = PP_fam(ind{1},2,d);
    pp2 = PP_fam(ind{2},2,d);
    PP_ksd{1,d} = ksdensity(pp1, x_bins);
    PP_ksd{2,d} = ksdensity(pp2, x_bins);
    PP_med(:,d) = [median(pp1), median(pp2)]';
end;

% visualize posterior probabilities
figure('Name', 'Analysis 1: complete model space', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
x_lim = [min(x_bins), max(x_bins)];

for d = 1:D
    subplot(2,2,d);
    hold on;
    if d == 1
        area(-1, -1, 'EdgeColor', 'r', 'FaceColor', [1, 2/3, 2/3], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
        area(-1, -1, 'EdgeColor', 'b', 'FaceColor', [2/3, 2/3, 1], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
    end;
    y_lim = [0, (11/10)*max(max([PP_ksd{1,d}; PP_ksd{2,d}]))];
    if d == 2, y_lim = [0, 10]; end;
    area(x_bins, PP_ksd{1,d}, 'EdgeColor', 'r', 'FaceColor', [1, 2/3, 2/3], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
    plot([PP_med(1,d), PP_med(1,d)], [0, PP_ksd{1,d}(round(x_bins,2)==round(PP_med(1,d),2))], 'r', 'LineWidth', 1);
    area(x_bins, PP_ksd{2,d}, 'EdgeColor', 'b', 'FaceColor', [2/3, 2/3, 1], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
    plot([PP_med(2,d), PP_med(2,d)], [0, PP_ksd{2,d}(round(x_bins,2)==round(PP_med(2,d),2))], 'b', 'LineWidth', 1);
    if d == 2
        pp1 = PP_fam(ind{1} & round(PP_fam(:,2,d),2)<1, 2, d);
        pp2 = PP_fam(ind{2} & round(PP_fam(:,2,d),2)<1, 2, d);
        plot(pp1, 2/3, '+r', 'LineWidth', 1);
        plot(pp2, 1/3, '+b', 'LineWidth', 1);
    end;
    axis([x_lim, y_lim]);
    set(gca,'Box','On');
    set(gca,'FontSize',12);
    if d == 1, legend(lab, 'Location', 'NorthEast'); end;
    xlabel('posterior probability', 'FontSize', 14);
    ylabel('probability density', 'FontSize', 14);
    title(l{d,3}, 'FontSize', 16);
    text(x_lim(1), max(y_lim)+1/100*y_lim(2), sprintf('EP = %2.2f %%', EP_fam(1,d)*100), ...
         'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'r', ...
         'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');
    text(x_lim(2), max(y_lim)+1/100*y_lim(2), sprintf('EP = %2.2f %%', EP_fam(2,d)*100), ...
         'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'b', ...
         'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
end;


%%% Figure 6: model parameter estimates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract parameter estimates
koi   = k(strcmp(m,moi));
p_lab = p_lab{strcmp(m,moi)};
p_MAP = zeros(N,koi);
for i = 1:N
    ab_est= ab_post{i,strcmp(m,moi)};
    for j = 1:koi
        p_MAP(i,j) = (ab_est(j,1)-1)/(ab_est(j,1)+ab_est(j,2)-2);
    end;        
end;
clear ab_est

% clean parameters estimates
p_ok  = sum(isnan(p_MAP),2)==0;
p_MAP = p_MAP(p_ok,:);
i1    = ind{1}(p_ok);
i2    = ind{2}(p_ok);

% estimate kernel densities
p_med = zeros(2,koi);
p_ksd = cell(2,koi);
for j = 1:koi
    y1 = p_MAP(i1,j);
    y2 = p_MAP(i2,j);
    p_ksd{1,j} = ksdensity(y1, x_bins);
    p_ksd{2,j} = ksdensity(y2, x_bins);
    p_med(:,j) = [median(y1), median(y2)]';
end;

% visualize model parameters
figure('Name', 'Analysis 3: model parameter estimates', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
x_lim = [min(x_bins), max(x_bins)];
sp = [1, 5, 7, 2, 4, 6, 8];
at = {'non-neutral (old & new)', 'affirmative (old)', 'affirmative (new)', ...
      'confident (old, affirmative)', 'confident (new, affirmative)', ...
      'confident (old, non-affirmative)', 'confident (new, non-affirmative)'};

for j = 1:koi
    if j == 1
        subplot(8,2,[1 3 5]);
    else
        subplot(4,2,sp(j));
    end;
    hold on;
    if j == 1
        area(-1, -1, 'EdgeColor', 'r', 'FaceColor', [1, 2/3, 2/3], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
        area(-1, -1, 'EdgeColor', 'b', 'FaceColor', [2/3, 2/3, 1], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
    end;
    y_lim = [0, (11/10)*max(max([p_ksd{1,j}; p_ksd{2,j}]))];
    area(x_bins, p_ksd{1,j}, 'EdgeColor', 'r', 'FaceColor', [1, 2/3, 2/3], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
    plot([p_med(1,j), p_med(1,j)], [0, p_ksd{1,j}(round(x_bins,2)==round(p_med(1,j),2))], 'r', 'LineWidth', 1);
    area(x_bins, p_ksd{2,j}, 'EdgeColor', 'b', 'FaceColor', [2/3, 2/3, 1], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
    plot([p_med(2,j), p_med(2,j)], [0, p_ksd{2,j}(round(x_bins,2)==round(p_med(2,j),2))], 'b', 'LineWidth', 1);
    axis([x_lim, y_lim]);
    set(gca,'Box','On');
    set(gca,'FontSize',12);
    if j == 1, legend(lab, 'Location', 'NorthWest'); end;
    xlabel(p_lab{j}, 'FontSize', 14);
    if j == 1, ylabel('probability density', 'FontSize', 14); end;
    if j  > 1, ylabel(sprintf('p(%s)', p_lab{j}), 'FontSize', 14); end;
    title(at{j}, 'FontSize', 16);
end;


%%% Figures 7/8: group comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble parameter estimates
i1 = ind{1}(p_ok);
i2 = ind{2}(p_ok);
Y1 = p_MAP;
l1 = p_lab'; % = {'p', 'q_o', 'q_n', 'r_o', 'r_n', 's_o', 's_n'};

% extract parameter estimates
p_a = Y1(:,1);
q_o = Y1(:,2);
q_n = Y1(:,3);
r_o = Y1(:,4);
r_n = Y1(:,5);
s_o = Y1(:,6);
s_n = Y1(:,7);

% calculate performance measures
Y1  = [Y1, p_a];                            % level of decidedness
Y1  = [Y1, q_o, ...                         % memory performance
           1-q_n, ...
           (q_o + (1-q_n))./2];
Y1  = [Y1, (r_o + s_o)./2, ...              % confidence level
           (r_n + s_n)./2, ...
           ((r_o + s_o)./2 + (r_n + s_n)./2)./2];
Y1  = [Y1, (r_o + (1-s_o))./2, ...          % meta-cognitive accuracy
           ((1-r_n) + s_n)./2, ...
           ((r_o + (1-s_o))./2 + ((1-r_n) + s_n)./2)./2];
l1  = [l1, {'LD = Pr("dec")_{ }', ...
            'MP_o',  'MP_n',  'MP = Pr("corr")', ...
            'LC_o',  'LC_n',  'LC = Pr("conf")', ...
            'MCA_o', 'MCA_n', '        MCA = Pr("conf=corr")'}];
clear p_a q_o q_n r_o r_n s_o s_n

% test parameter estimates
k1  = size(Y1,2);
p1  = zeros(1,k1);
p1c = zeros(1,k1);
z1  = zeros(1,k1);
m1  = zeros(2,k1);
ci1 = zeros(2,k1);
med1= zeros(2,k1);
for j = 1:k1
    y1 = Y1(i1,j);
    y2 = Y1(i2,j);
   [p1(j), h, stats] = ranksum(y1, y2);
    z1(j)     = stats.zval;
    m1(:,j)   = [mean(y1), mean(y2)]';
    ci1(:,j)  = [std(y1)/sqrt(numel(y1)), std(y2)/sqrt(numel(y2))]' * norminv(mean([conf_lvl 1]), 0, 1);
    med1(:,j) = [median(y1), median(y2)]';
end;
p0 = p1;
clear h stats

% visualize parameter estimates
figure('Name', sprintf('Analysis 4: significance tests (%d)', 1), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
ii = [1:7];
p1c(ii) = bonf_holm(p1(ii), 0.05);
hold on;
bp = bar(m1(:,ii)', 'grouped');
xo = 1/7; % abs(get(bp(1),'XOffset'));
set(bp(1), 'FaceColor', col(1));
set(bp(2), 'FaceColor', col(2));
errorbar([1:numel(ii)]-xo, m1(1,ii), ci1(1,ii), '.k', 'LineWidth', 2, 'CapSize', 15);
errorbar([1:numel(ii)]+xo, m1(2,ii), ci1(2,ii), '.k', 'LineWidth', 2, 'CapSize', 15);
axis([(1-0.5), (numel(ii)+0.5), (0-0.05), (1+0.05)]);
set(gca, 'Box', 'On');
set(gca, 'FontSize', 12);
set(gca, 'XTick', [1:numel(l1(ii))], 'XTickLabel', l1(ii), 'XTickLabelRotation', 0);
legend(lab, 'Location', 'SouthWest');
xlabel('model parameter', 'FontSize', 14);
ylabel(sprintf('mean and %d%% CI', round(conf_lvl*100)), 'FontSize', 14);
for j = 1:numel(ii)
    text(j, 0.95, repmat('*',[1 sum(p1c(ii(j))<p_thr)]), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    if p1c(ii(j)) < p_thr(end)
        text(j, 1, sprintf('z = %0.2f\np < 0.001', z1(ii(j))), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    else
        text(j, 1, sprintf('z = %0.2f\np = %0.3f', z1(ii(j)), p1c(ii(j))), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
end;

% visualize behavioral performance
figure('Name', sprintf('Analysis 4: significance tests (%d)', 2), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
ii = [8,11,14,17,9,10,12,13,15,16];
p1c(ii) = bonf_holm(p1(ii), 0.05);

for i = 1:2
    if i == 1, ii = [8,11,14,17];       end;
    if i == 2, ii = [9,10,12,13,15,16]; end;
    subplot(1,2,i); hold on;
    bp = bar(m1(:,ii)', 'grouped');
    xo = 1/7; % abs(get(bp(1),'XOffset'));
    set(bp(1), 'FaceColor', col(1));
    set(bp(2), 'FaceColor', col(2));
    errorbar([1:numel(ii)]-xo, m1(1,ii), ci1(1,ii), '.k', 'LineWidth', 2, 'CapSize', 15-(i-1)*5);
    errorbar([1:numel(ii)]+xo, m1(2,ii), ci1(2,ii), '.k', 'LineWidth', 2, 'CapSize', 15-(i-1)*5);
    axis([(1-0.5), (numel(ii)+0.5), (0-0.05), (1+0.05)]);
    set(gca, 'Box', 'On');
    set(gca, 'FontSize', 12);
    set(gca, 'XTick', [1:numel(l1(ii))], 'XTickLabel', l1(ii), 'XTickLabelRotation', 0);
    if i == 1, legend(lab, 'Location', 'SouthWest'); end;
    if i == 1, xlabel('behavioral performance measure', 'FontSize', 14); end;
    if i == 2, xlabel('behavioral performance, separately for old and new items', 'FontSize', 14); end;
    ylabel(sprintf('mean and %d%% CI', round(conf_lvl*100)), 'FontSize', 14);
    for j = 1:numel(ii)
        text(j, 0.95, repmat('*',[1 sum(p1c(ii(j))<p_thr)]), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        if p1c(ii(j)) < p_thr(end)
            text(j, 1, sprintf('z = %0.2f\np < 0.001', z1(ii(j))), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        else
            text(j, 1, sprintf('z = %0.2f\np = %0.3f', z1(ii(j)), p1c(ii(j))), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
end;


%%% Figure 9: correlation & regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare parameter extraction
ii = [8, 11, 14, 17];
l2 = {'LD', 'MP', 'LC', 'MCA'};
v2 = numel(ii);

% prepare regression lines
x1 = cell(v2,v2);               % parameter estimates
x2 = cell(v2,v2);
y1 = cell(v2,v2);               % parameter estimates
y2 = cell(v2,v2);
mn1= cell(v2,v2);               % regression lines
mn2= cell(v2,v2);
r1 = zeros(v2,v2);              % correlation coefficients
r2 = zeros(v2,v2);
p1 = zeros(v2,v2);              % correlation p-values
p2 = zeros(v2,v2);
d1 = cell(1,v2);                % kernel density estimates
d2 = cell(1,v2);

% estimate regression lines
for j1 = 1:v2
    for j2 = 1:(j1-1)
        % correlation within young subjects
        y1{j1,j2} = Y1(i1,ii(j1));
        x1{j1,j2} = Y1(i1,ii(j2));
        [r1(j1,j2), p1(j1,j2)] = corr(y1{j1,j2}, x1{j1,j2}, 'type', 'Spearman');
        mn1{j1,j2} = polyfit(x1{j1,j2}, y1{j1,j2}, 1);
        % correlation within older subjects
        y2{j1,j2} = Y1(i2,ii(j1));
        x2{j1,j2} = Y1(i2,ii(j2));
        [r2(j1,j2), p2(j1,j2)] = corr(y2{j1,j2}, x2{j1,j2}, 'type', 'Spearman');
        mn2{j1,j2} = polyfit(x2{j1,j2}, y2{j1,j2}, 1);
    end;
end;

% estimate kernel densities
for j = 1:v2
    d1{j} = ksdensity(Y1(i1,ii(j)), x_bins);
    d2{j} = ksdensity(Y1(i2,ii(j)), x_bins);
    d1{j} = d1{j}./trapz(x_bins, d1{j});
    d2{j} = d2{j}./trapz(x_bins, d2{j});
end;

% visualize performance measures
figure('Name', 'Analysis 6: regression lines', 'Color', [1 1 1], 'Position', [50 50 1500 900]);
at    = {'level of decidedness', 'memory performance', 'level of confidence', 'meta-cognitive accuracy'};
x_lim = [0,1];
y_lim = [0,1];

for j1 = 1:v2
    for j2 = 1:j1
        subplot(v2,v2,(j1-1)*v2+j2); hold on;
        % scatter plot
        if j2 < j1
            plot(x1{j1,j2}, y1{j1,j2}, '.r', 'Color', [1, (p1(j1,j2)>p_thr(1))*2/3, (p1(j1,j2)>p_thr(1))*2/3]);
            plot(x2{j1,j2}, y2{j1,j2}, '.b', 'Color', [(p2(j1,j2)>p_thr(1))*2/3, (p2(j1,j2)>p_thr(1))*2/3, 1]);
            plot(x_lim, mn1{j1,j2}(1)*x_lim + mn1{j1,j2}(2), '-r', 'Color', [1, (p1(j1,j2)>p_thr(1))*2/3, (p1(j1,j2)>p_thr(1))*2/3]);
            plot(x_lim, mn2{j1,j2}(1)*x_lim + mn2{j1,j2}(2), '-b', 'Color', [(p2(j1,j2)>p_thr(1))*2/3, (p2(j1,j2)>p_thr(1))*2/3, 1]);
            axis([x_lim, y_lim]);
            set(gca,'Box','On');
          % set(gca,'FontSize',12);
            xlabel(l2{j2}, 'FontSize', 14);
            ylabel(l2{j1}, 'FontSize', 14);
            if p1(j1,j2) < p_thr(end)
                txt1 = sprintf('\\rho = %0.2f, p < 0.001', r1(j1,j2));
            else
                txt1 = sprintf('\\rho = %0.2f, p = %0.3f', r1(j1,j2), p1(j1,j2));
            end;
            if p2(j1,j2) < p_thr(end)
                txt2 = sprintf('\\rho = %0.2f, p < 0.001', r2(j1,j2));
            else
                txt2 = sprintf('\\rho = %0.2f, p = %0.3f', r2(j1,j2), p2(j1,j2));
            end;
            text(x_lim(1)-0.01, max(y_lim)+1/25*range(y_lim), txt1, ...
                 'FontSize', 10, 'FontWeight', 'Bold', 'Color', [1, (p1(j1,j2)>p_thr(1))*2/3, (p1(j1,j2)>p_thr(1))*2/3], ...
                 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');
            text(x_lim(2)+0.01, max(y_lim)+1/25*range(y_lim), txt2, ...
                 'FontSize', 10, 'FontWeight', 'Bold', 'Color', [(p2(j1,j2)>p_thr(1))*2/3, (p2(j1,j2)>p_thr(1))*2/3, 1], ...
                 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
        % density plot
        else
            area(x_bins, d1{j2}, 'EdgeColor', 'r', 'FaceColor', [1, 2/3, 2/3], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
            area(x_bins, d2{j2}, 'EdgeColor', 'b', 'FaceColor', [2/3, 2/3, 1], 'FaceAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1);
            axis([x_lim, 0, (11/10)*max(max([d1{j1}; d2{j1}]))]);
            set(gca,'Box','On');
          % set(gca,'FontSize',12);
            if j2 == 1, legend(lab, 'FontSize', 12, 'Location', 'NorthWest'); end;
            xlabel(l2{j2}, 'FontSize', 14);
            ylabel(sprintf('p(%s)', l2{j1}), 'FontSize', 14);
            title(at{j2}, 'FontSize', 16);
        end;
    end;
end;


%%% Figure 10: results summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% behavioral performance measures
j_ind = [8, 11, 14, 17];
p_thr = [0.05, 0.01, 0.001];

% report median differences
fprintf('\n-> Age-related differences:\n');
for j = 1:v2
    p_str = repmat('*',[1 sum(p0(j_ind(j))<p_thr)]);
    fprintf('   - %s: delta = %0.2f%s.\n', l2{j}, med1(1,j_ind(j))-med1(2,j_ind(j)), p_str);
end;

% report correlation coefficients
fprintf('\n-> Inter-measure correlations:\n');
for j1 = 1:v2
    for j2 = 1:(j1-1)
        if p1(j1,j2) < p_thr(1) || p2(j1,j2) < p_thr(1)
            p1_str = repmat('*',[1 sum(p1(j1,j2)<p_thr)]);
            p2_str = repmat('*',[1 sum(p2(j1,j2)<p_thr)]);
            fprintf('   - %s vs. %s: rho = %0.2f%s, %0.2f%s.\n', l2{j1}, l2{j2}, r1(j1,j2), p1_str, r2(j1,j2), p2_str);
        end;
    end;
end;
fprintf('\n');
clear j_ind p_str p1_str p2_str