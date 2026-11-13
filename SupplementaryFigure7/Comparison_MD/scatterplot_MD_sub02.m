function scatterplot_MD_sub02

% This script aims to replicate a panel in Supplementary Figure 7 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Define the base directory and subject information
baseDir = '../../TractometryData/Protocol1/';
subjects = {'sub-02'};

%% Define session information and property file
sessionPairs = {
    {'ses-PrismaFitB1','ses-PrismaFitB2'}
    {'ses-PrismaFitB1','ses-PrismaFitA1'}
    {'ses-PrismaFitB1','ses-VerioA1'}
};

prop = 'md.csv';
propCols = 1:20;

%% Add path to external tools
addpath('../../ExternalTools');

%% Count the number
nSub   = numel(subjects);
nPairs = numel(sessionPairs);
nProp  = numel(propCols);

%% Define marker/color for sub-02
subColor  = [0.8500 0.3250 0.0980];
subMarker = 's';

fig = figure('Position',[100,100,1000,550]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

%% Main loop
for pairIdx = 1:nPairs
    sesa = sessionPairs{pairIdx}{1};
    sesb = sessionPairs{pairIdx}{2};

    ax = nexttile(pairIdx);
    hold(ax, 'on');

    R = nan(nSub,1);
    Xa_all = nan(nProp, nSub);
    Xb_all = nan(nProp, nSub);

    for subIdx = 1:nSub
        subjectID = subjects{subIdx};
        fileA = fullfile(baseDir, subjectID, sesa, prop);
        fileB = fullfile(baseDir, subjectID, sesb, prop);

        A = readmatrix(fileA);
        B = readmatrix(fileB);
        Xa = mean(A(1:100, propCols), 1).';
        Xb = mean(B(1:100, propCols), 1).';

        Xa_all(:,subIdx) = Xa;
        Xb_all(:,subIdx) = Xb;

        scatter(ax, Xa, Xb, 18, ...
    'Marker', subMarker, ...
    'MarkerEdgeColor', subColor, ...
    'MarkerFaceColor', 'none');

        T = [Xa, Xb];
        [r, ~, ~, ~, ~, ~, ~] = ICC(T, 'A-1', 0.05, 0);
        R(subIdx) = r;
    end

    %% Axis settings
    axis(ax, 'square');
    M = 1.0;
    m = 0.6;
    set(ax, 'TickDir', 'out', 'Box', 'off', ...
        'YLim', [m M], 'YTick', [m (m+M)/2 M], ...
        'XLim', [m M], 'XTick', [m (m+M)/2 M], ...
        'FontSize', 15);

    plot(ax, [m M], [m M], '--k', 'LineWidth', 1.5);

    %% Panel title
    ttlStr = panel_title(pairIdx);
    title(ax, ttlStr, ...
        'FontSize', 12, ...
        'FontSizeMode', 'manual', ...
        'Units', 'normalized', ...
        'Position', [0.52 0.91 0], ...
        'HorizontalAlignment', 'center');

    %% Labels
    [xlblStr, ylblStr] = axis_labels(pairIdx);
    xlabel(ax, xlblStr, 'FontSize', 15);
    ylabel(ax, ylblStr, 'FontSize', 15);

    %% ICC text
    meanICC = mean(R, 'omitnan');
    text(ax, 1.05, 0.02, sprintf('ICC = %.3f', meanICC), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 14, ...
        'Clipping', 'off');

    %% Regression line
    x = Xa_all(:);
    y = Xb_all(:);
    valid = ~isnan(x) & ~isnan(y);

    if sum(valid) >= 2
        p = polyfit(x(valid), y(valid), 1);
        xx = linspace(min(x(valid)), max(x(valid)), 200);
        yy = polyval(p, xx);
        plot(ax, xx, yy, '-', 'Color', [0 0 0 0.6], 'LineWidth', 1.5);
    end
end

%% output figure as eps
exportgraphics(fig, 'scatterplot_MD_sub02.eps');

end

%% Title definition
function ttl = panel_title(pairIdx)
switch pairIdx
    case 1
        ttl = {'Same Scanner', ''};
    case 2
        ttl = {'Different Scanners,', 'Same Scanner Model'};
    case 3
        ttl = {'Different Scanners,', 'Different Scanner Models'};
end
end

%% Label definition
function [xlbl, ylbl] = axis_labels(pairIdx)
switch pairIdx
    case 1
        xlbl = 'Protocol1 - PrismaFit B [Run1]';
        ylbl = 'Protocol1 - PrismaFit B [Run2]';
    case 2
        xlbl = 'Protocol1 - PrismaFit B [Run1]';
        ylbl = 'Protocol1 - PrismaFit A [Run1]';
    case 3
        xlbl = 'Protocol1 - PrismaFit B [Run1]';
        ylbl = 'Protocol1 - Verio A [Run1]';
end
end