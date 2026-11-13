function ICCheatmap_ODI_sub02

% This script aims to replicate a panel in Supplementary Figure 10 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI


%% Define base directories and subject information
baseDirPairs = {
    '../../TractometryData/Protocol1/'
    '../../TractometryData/Protocol2/'
};

subjects = {'sub-02'};

%% Define sessions
sessions_p = { ...
    'ses-VerioA1','ses-VerioA2', ...
    'ses-VerioB1','ses-VerioB2', ...
    'ses-PrismaFitA1','ses-PrismaFitA2', ...
    'ses-PrismaFitB1','ses-PrismaFitB2'};

% Heatmap labels
sessions = [ ...
    strcat('Protocol1-', sessions_p), ...
    strcat('Protocol2-', sessions_p) ...
];

prop = 'odi.csv';
propCols = 1:20;

nFeat = numel(propCols);
nSub  = numel(subjects);
nSes  = numel(sessions);

%% Add external tools to the path
addpath('../../ExternalTools');

%% 1) Read all files once: data(sub, ses, feat)
data = nan(nSub, nSes, nFeat);

for sesIdx = 1:nSes
    if startsWith(sessions{sesIdx}, 'Protocol1-')
        baseDir = baseDirPairs{1};
        ses     = extractAfter(sessions{sesIdx}, 'Protocol1-');
    else
        baseDir = baseDirPairs{2};
        ses     = extractAfter(sessions{sesIdx}, 'Protocol2-');
    end

    for subIdx = 1:nSub
        subjectID = subjects{subIdx};
        filePath  = fullfile(baseDir, subjectID, ses, prop);

        if isfile(filePath)
            A = readmatrix(filePath);
            data(subIdx, sesIdx, :) = mean(A(1:100, propCols), 1);
        end
    end
end

%% 2) Compute session-to-session ICC (Original)
R_orig = nan(nSes, nSes);

for i = 1:nSes
    R_orig(i,i) = 1;
    for j = i+1:nSes

        all_r = nan(nSub,1);
        type  = 'A-1';
        alpha = 0.05;
        r0    = 0;

        for k = 1:nSub
            xi = squeeze(data(k, i, :))';   % 1 × nFeat
            xj = squeeze(data(k, j, :))';   % 1 × nFeat

            % Use only tracts with valid values in both sessions
            valid = ~isnan(xi) & ~isnan(xj);
            if sum(valid) < 2
                continue
            end

            T = [xi(valid)', xj(valid)'];   % nValidFeat × 2
            [r, ~, ~, ~, ~, ~, ~] = ICC(T, type, alpha, r0);
            all_r(k) = r;
        end

        R_orig(i,j) = mean(all_r, 'omitnan');
    end
end

R_orig = triu(R_orig) + triu(R_orig,1)';

%% 3) Plot ICC heatmap
fig = figure('Position',[100,100,1000,900], 'Renderer','opengl');

h1 = heatmap(sessions, sessions, R_orig);
h1.Title = 'ICC Heatmap for ODI';
h1.XLabel = 'Session';
h1.YLabel = 'Session';
h1.Colormap = jet;
h1.ColorLimits = [0 1];
h1.CellLabelFormat = '%.2f';

%% Output figure as eps
exportgraphics(fig, 'ICCheatmap_ODI_sub02.eps');

end