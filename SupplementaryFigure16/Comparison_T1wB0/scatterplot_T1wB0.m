function scatterplot_T1wB0

% This script aims to replicate a panel in Supplementary Figure 16 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Define subject information
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};

%% Define directory/session pairs and property file
baseDirPairs = {
    {'../../TractometryData/Protocol1/','../../TractometryData/Protocol1/'}
    {'../../TractometryData/Protocol1/','../../TractometryData/Protocol1/'}
    {'../../TractometryData/Protocol1/','../../TractometryData/Protocol2/'}
};

sessionPairs = {
    {'ses-VerioA1','ses-VerioA2'}
    {'ses-VerioA1','ses-PrismaFitA1'}
    {'ses-VerioA1','ses-PrismaFitA1'}
};

prop = 't1w_over_b0.csv';
propCols = 1:20;

%% Add path to external tools
addpath(genpath('../../ExternalTools'))

%% Count the number
nSub   = numel(subjects);
nPairs = numel(sessionPairs);
nFeat  = numel(propCols);

%% Figure
fig = figure('Position',[100,100,1200,750], 'Renderer','opengl');
tiledlayout(2,3,'Padding','loose','TileSpacing','compact');

% Axis range
M = 2.2;
m = 0;

% Subject markers
subMarkers = {'o','s','^','d','p'};

%% Main loop
for pairIdx = 1:nPairs
    baseDirA = baseDirPairs{pairIdx}{1};
    baseDirB = baseDirPairs{pairIdx}{2};
    sesa = sessionPairs{pairIdx}{1};
    sesb = sessionPairs{pairIdx}{2};

    % =======================
    % Top row: Original
    % =======================
    ax1 = nexttile(pairIdx);
    hold(ax1,'on');

    R_orig = nan(nSub,1);
    Xa_all = nan(nFeat, nSub);
    Xb_all = nan(nFeat, nSub);

    for subIdx = 1:nSub
        subjectID = subjects{subIdx};

        fileA = fullfile(baseDirA, subjectID, sesa, prop);
        fileB = fullfile(baseDirB, subjectID, sesb, prop);

        A = readmatrix(fileA);
        B = readmatrix(fileB);

        Xa = mean(A(1:100,propCols),1).';
        Xb = mean(B(1:100,propCols),1).';

        Xa_all(:,subIdx) = Xa;
        Xb_all(:,subIdx) = Xb;

        mk = subMarkers{subIdx};
        scatter(ax1, Xa, Xb, 20, 'Marker', mk);

        T = [Xa, Xb];
        [r, ~, ~, ~, ~, ~, ~] = ICC(T, 'A-1', 0.05, 0);
        R_orig(subIdx) = r;
    end

    axis(ax1,'square');
    set(ax1,'tickdir','out','box','off', ...
        'ylim',[m M],'ytick',[m (m+M)/2 M], ...
        'xlim',[m M],'xtick',[m (m+M)/2 M], ...
        'fontsize',16);
    plot(ax1,[m M],[m M],'--k','LineWidth',1.5);

    ttlStr = panel_title(pairIdx);
    title(ax1, ttlStr, ...
        'FontSize', 14, ...
        'Units', 'normalized', ...
        'Position', [0.52 0.91 0], ...
        'HorizontalAlignment', 'center');

    [xlbl, ylbl] = axis_labels(pairIdx);
    xlabel(ax1, xlbl, 'FontSize',16);
    ylabel(ax1, ylbl, 'FontSize',16);

    meanICC = mean(R_orig,'omitnan');
    text(ax1, 1.05, 0.02, sprintf('Mean ICC = %.3f', meanICC), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontSize',15,'Clipping','off');

    % Regression line
    x = Xa_all(:);
    y = Xb_all(:);
    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 200);
    yy = polyval(p, xx);
    plot(ax1, xx, yy, '-', 'Color', [0 0 0 0.6], 'LineWidth', 1.5);

    % =======================
    % Bottom row: ComBat-harmonized
    % =======================
    ax2 = nexttile(pairIdx+3);

    % Do not apply ComBat for the same-condition comparison
    if pairIdx == 1
        axis(ax2, 'off');
        continue
    end
    
    hold(ax2,'on');

    all_data = [Xa_all, Xb_all];
    batch_labels = [ones(nSub,1); 2*ones(nSub,1)];

    data_harm = combat(all_data, batch_labels, [], 1)';

    R_cb = nan(nSub,1);
    Xh1_all = nan(nFeat, nSub);
    Xh2_all = nan(nFeat, nSub);

    for subIdx = 1:nSub
        Xh1 = data_harm(subIdx,      :)';
        Xh2 = data_harm(subIdx+nSub, :)';

        Xh1_all(:,subIdx) = Xh1;
        Xh2_all(:,subIdx) = Xh2;

        mk = subMarkers{subIdx};
        scatter(ax2, Xh1, Xh2, 20, 'Marker', mk);

        T = [Xh1, Xh2];
        [r, ~, ~, ~, ~, ~, ~] = ICC(T, 'A-1', 0.05, 0);
        R_cb(subIdx) = r;
    end

    % Regression line
    x = Xh1_all(:);
    y = Xh2_all(:);
    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 200);
    yy = polyval(p, xx);
    plot(ax2, xx, yy, '-', 'Color', [0 0 0 0.6], 'LineWidth', 1.5);

    axis(ax2,'square');
    set(ax2,'tickdir','out','box','off', ...
        'ylim',[m M],'ytick',[m (m+M)/2 M], ...
        'xlim',[m M],'xtick',[m (m+M)/2 M], ...
        'fontsize',16);
    plot(ax2,[m M],[m M],'--k','LineWidth',1.5);

    [xlbl, ylbl] = axis_labels(pairIdx);
    xlabel(ax2, xlbl, 'FontSize',16);
    ylabel(ax2, ylbl, 'FontSize',16);

    meanICC = mean(R_cb,'omitnan');
    text(ax2, 1.05, 0.02, sprintf('Mean ICC = %.3f', meanICC), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontSize',15,'Clipping','off');
end

%% output figure as eps
exportgraphics(fig, 'scatterplot_T1wB0.eps');

end

%% Title definition
function ttl = panel_title(pairIdx)
switch pairIdx
    case 1
        ttl = {'Same Scanner,', 'Same Protocol'};
    case 2
        ttl = {'Different Scanners,', 'Same Protocol'};
    case 3
        ttl = {'Different Scanners,', 'Different Protocols'};
end
end

%% Label definition
function [xlbl, ylbl] = axis_labels(pairIdx)
switch pairIdx
    case 1
        xlbl = 'Protocol1 - Verio A [Run1]';
        ylbl = 'Protocol1 - Verio A [Run2]';
    case 2
        xlbl = 'Protocol1 - Verio A [Run1]';
        ylbl = 'Protocol1 - PrismaFit A [Run1]';
    case 3
        xlbl = 'Protocol1 - Verio A [Run1]';
        ylbl = 'Protocol2 - PrismaFit A [Run1]';
end
end