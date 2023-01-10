% FADE-pqrs: Figure 2B
% written   by Joram Soch <Joram.Soch@DZNE.de>, 19/11/2021, 12:02
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/01/2023, 14:24

% clear
% close all

% define parameters
N = 15;                         % number of subjects
n = 100;                        % number of trials
o = 1;                          % fraction of old items
p = 0.9;
q = 2/3;
r = 2/3;
p     = [p, q, r]';
p_lab = {'p', 'q', 'r'}';

% sample data
rng(1);
Y = zeros(n,N);
C = zeros(5,N);
for i = 1:N
    [Y(:,i), x, m] = ME_pqrs_Sim(p, p_lab, n, o);
    C(:,i) = hist(Y(:,i), [1:5])';
end;

% visualize counts
figure('Name', 'pqrs: Figure 2B', 'Color', [1 1 1], 'Position', [50 50 600 900]);
for i = 1:N
    subplot(N/3,3,i);
    bar([1:5], C(:,i)', 'FaceColor', 2/3*[1 1 1]);
    axis([(1-1), (5+1), 0, max(max(C))+1]);
    if i > N-3,     xlabel('response', 'FontSize', 12); end;
    if mod(i,3)==1, ylabel('count', 'FontSize', 12); end;
  % title(sprintf('Subject %d', i), 'FontSize', 12);
end;