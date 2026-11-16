function ScatterPlot_ScannerComparison_MD

% This script aims to replicate a panel in Figure 4 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Define the base directory and subject information
baseDir = '../../TractometryData/Protocol1/';
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};

%% Define session information and property file
sessionPairs = {
    {'ses-VerioA1','ses-VerioA2'}
    {'ses-VerioA1','ses-VerioB1'}
    {'ses-VerioA1','ses-PrismaFitA1'}
};
prop = 'md.csv';
propCols = 1:20;

%% Add path to external tools
addpath(genpath('../../ExternalTools'))

%% Count the number
nSub   = numel(subjects);
nPairs = numel(sessionPairs);
nProp  = numel(propCols);

fig = figure('Position',[100,100,1000,550]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

%% Main loop
for pairIdx = 1:nPairs
    sesa = sessionPairs{pairIdx}{1};
    sesb = sessionPairs{pairIdx}{2};

    ax = nexttile(pairIdx);
    % Panel title
    ttlStr = panel_title(pairIdx);
    title(ax, ttlStr, ...
        'FontSize', 15, ...
        'Units', 'normalized', ...
        'Position', [0.1 0.89 0], ...
        'HorizontalAlignment', 'left');

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
        Xa = mean(A(1:100,propCols),1).';
        Xb = mean(B(1:100,propCols),1).';

        Xa_all(:,subIdx) = Xa;
        Xb_all(:,subIdx) = Xb;

        subMarkers = {'o','s','^','d','p'};   % circle, square, triangle up, diamond, triangle down
        mk = subMarkers{subIdx};
        scatter(Xa, Xb, 18, 'Marker', mk);
        T = [Xa, Xb];
        [r, ~, ~, ~, ~, ~, ~] = ICC(T, 'A-1', 0.05, 0);
        R(subIdx) = r;
    end

    % Axis settings
    axis square; M=1.0; m=0.6;
    set(gca,'tickdir','out','box','off',...
        'ylim',[m M],'ytick',[m (m+M)/2 M],...
        'xlim',[m M],'xtick',[m (m+M)/2 M],'fontsize',15);
    plot([m M],[m M],'--k','LineWidth',1.5);

    % Panel title
    ttlStr = panel_title(pairIdx);
    h = title(ax, ttlStr, ...
        'FontSize', 12, ...
        'FontSizeMode', 'manual', ...
        'Units', 'normalized', ...
        'Position', [0.52 0.91 0], ...
        'HorizontalAlignment', 'center');

    % Labels
    [xlblStr, ylblStr] = axis_labels(pairIdx);
    hx = xlabel(ax, xlblStr, 'FontSize',15);
    hy = ylabel(ax, ylblStr, 'FontSize',15);

    % average ICC across subjects
    meanICC = mean(R,'omitnan');
    text(1.05, 0.02, sprintf('Mean ICC = %.3f', meanICC), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontSize',14,'Clipping','off');

    % Regression line (across all points in this panel)
    x = Xa_all(:);
    y = Xb_all(:);
    p = polyfit(x, y, 1);
    xmin = min(x);
    xmax = max(x);
    xx = linspace(xmin, xmax, 200);
    yy = polyval(p, xx);
    plot(ax, xx, yy, '-', ...
        'Color', [0 0 0 0.6], ...
        'LineWidth', 1.5);
end

%% output figure as eps
exportgraphics(fig, 'ScatterPlot_ScannerComparison_MD.eps');

end

%% Title deffinition
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

%% Label deffinition
function [xlbl, ylbl] = axis_labels(pairIdx)
switch pairIdx
    case 1
        xlbl = 'Protocol1 - Verio A [Run1]';
        ylbl = 'Protocol1 - Verio A [Run2]';
    case 2
        xlbl = 'Protocol1 - Verio A [Run1]';
        ylbl = 'Protocol1 - Verio B [Run1]';
    case 3
        xlbl = 'Protocol1 - Verio A [Run1]';
        ylbl = 'Protocol1 - PrismaFit A [Run1]';
end
end

