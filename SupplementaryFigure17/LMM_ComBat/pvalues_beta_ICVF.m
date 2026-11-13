function pvalues_beta_ICVF

% This script aims to replicate a panel in Supplementary Figure 17 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

baseDirs = {'../../TractometryData/Protocol1', ...
            '../../TractometryData/Protocol2'};
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};
sessions_p = {'ses-VerioA1','ses-VerioA2','ses-VerioB1', ...
              'ses-VerioB2','ses-PrismaFitA1','ses-PrismaFitA2'};
prop = 'icvf.csv';
metricLabel = 'ICVF';
nSub = numel(subjects);
nP = numel(sessions_p);
sessions = [strcat('Protocol1-',sessions_p), ...
            strcat('Protocol2-',sessions_p)];
nSes = numel(sessions);

addpath(genpath('../../ExternalTools'))

% Read tract-wise means: subject x session x tract
data = nan(nSub,nSes,20);
for p = 1:2
    for k = 1:nP
        c = (p-1)*nP + k;
        for s = 1:nSub
            A = readmatrix(fullfile(baseDirs{p},subjects{s}, ...
                                    sessions_p{k},prop));
            data(s,c,:) = mean(A(1:100,1:20),1);
        end
    end
end

% Crossed LMM: change in absolute discrepancy after pairwise ComBat
pValues = nan(nSes);
[subIdx,tractIdx] = ndgrid(1:nSub,1:20);
Subject = categorical(subIdx(:));
Tract = categorical(tractIdx(:));

for i = 1:nSes
    Xi = reshape(data(:,i,:),nSub,20);
    for j = i+1:nSes
        % Skip repeated scans with the same protocol and scanner.
        if ceil(i/nP) == ceil(j/nP) && ...
                strcmp(regexprep(sessions_p{mod(i-1,nP)+1},'\d+$',''), ...
                       regexprep(sessions_p{mod(j-1,nP)+1},'\d+$',''))
            continue
        end

        Xj = reshape(data(:,j,:),nSub,20);
        X = [Xi;Xj];
        Xcb = combat(X',[ones(nSub,1);2*ones(nSub,1)],[],1)';
        D = abs(Xcb(1:nSub,:) - Xcb(nSub+1:end,:)) - abs(Xi-Xj);

        tbl = table(D(:),Subject,Tract, ...
                    'VariableNames',{'Delta','Subject','Tract'});
        lme = fitlme(tbl,'Delta ~ 1 + (1|Subject) + (1|Tract)', ...
                     'FitMethod','REML');
        [~,~,stats] = fixedEffects(lme,'DFMethod','satterthwaite');
        pValues(i,j) = 2*tcdf(-abs(stats.tStat(1)),stats.DF(1));
    end
end
pValues = triu(pValues,1);
pValues = pValues + pValues';
pValues(1:nSes+1:end) = NaN;

% Two-sided p-value heatmap
fig = figure('Position',[100 100 900 760],'Color','w');
ax = axes(fig);
h = imagesc(ax,pValues);
h.AlphaData = isfinite(pValues);
axis(ax,'square');
ax.Color = 'w';
colormap(ax,parula);
caxis(ax,[0 0.2]);
cb = colorbar(ax);
cb.Label.String = 'Two-sided p-value';
cb.Ticks = 0:0.05:0.2;
cb.TickLabels = {'0','.05','.10','.15','.20'};
title(ax,metricLabel,'FontSize',18,'FontWeight','bold');
xlabel(ax,'Session'); ylabel(ax,'Session');
ax.XTick = 1:nSes; ax.YTick = 1:nSes;
ax.XTickLabel = sessions; ax.YTickLabel = sessions;
ax.TickLabelInterpreter = 'none';
ax.XTickLabelRotation = 45;
ax.TickLength = [0 0];
ax.FontSize = 9;

for i = 1:nSes
    for j = 1:nSes
        p = pValues(i,j);
        if ~isfinite(p), continue; end
        if p <= 0.001
            label = '<.001';
        else
            label = regexprep(sprintf('%.3f',p),'^0','');
        end
        color = 'k';
        if p <= 0.05, color = 'w'; end
        text(ax,j,i,label,'HorizontalAlignment','center', ...
             'FontSize',10,'Color',color);
    end
end

exportgraphics(fig,sprintf('pvalues_beta_%s.eps', ...
                           metricLabel),'ContentType','vector');
end