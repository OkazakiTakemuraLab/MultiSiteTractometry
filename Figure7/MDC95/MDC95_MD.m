function MDC95_MD
% This script aims to replicate a panel in Figure 7 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T.,
% Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Add path to external tools
addpath(genpath('../../ExternalTools'))

%% ===== Base dirs / subjects / sessions =====
baseDirPairs = {
    '../../TractometryData/Protocol1/'
    '../../TractometryData/Protocol2/'
};

subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};

sessions_p = { ...
    'ses-VerioA1', ...
    'ses-VerioA2', ...
    'ses-VerioB1', ...
    'ses-VerioB2', ...
    'ses-PrismaFitA1', ...
    'ses-PrismaFitA2' ...
};

sessions = [ ...
    strcat('Protocol1-', sessions_p), ...
    strcat('Protocol2-', sessions_p) ...
];

prop     = 'md.csv';
propCols = 1:20;

nTracts = numel(propCols);
nSub    = numel(subjects);
nSes    = numel(sessions);

%% ===== Bootstrap settings =====
nBoot = 5000;
rng(1);

%% ===== 1) Read data =====
data_orig = nan(nSub, nSes, nTracts);

for sesIdx = 1:nSes

    sessName = sessions{sesIdx};

    if startsWith(sessName, 'Protocol1-')
        baseDir = baseDirPairs{1};
        ses     = extractAfter(sessName, 'Protocol1-');
    else
        baseDir = baseDirPairs{2};
        ses     = extractAfter(sessName, 'Protocol2-');
    end

    for subIdx = 1:nSub

        fp = fullfile(baseDir, subjects{subIdx}, ses, prop);

        A = readmatrix(fp);

        data_orig(subIdx, sesIdx, :) = mean(A(1:100, propCols), 1, 'omitnan');
    end
end

%% ===== 2) ComBat for full dataset =====
% Batch labels:
% Each scanner/protocol condition is treated as one batch.
% Run 1 and run 2 from the same scanner/protocol condition share the same batch.

batch_ses = [1 1 2 2 3 3 4 4 5 5 6 6];
batch = repmat(batch_ses, 1, nSub);

T = reshape(permute(data_orig, [3 2 1]), nTracts, []);

T_cb = combat(T, batch, [], 1);

data_cb = permute(reshape(T_cb, nTracts, nSes, nSub), [3 2 1]);

%% ===== 3) MDC95 for full dataset =====
typeICC = 'A-1';
alpha = 0.05;
r0 = 0;

[MDC_orig, MDC_cb] = compute_mdc95(data_orig, data_cb, typeICC, alpha, r0);

%% ===== 4) Bootstrap MDC95 CI =====
MDC_orig_boot = nan(nBoot, nTracts);
MDC_cb_boot   = nan(nBoot, nTracts);

for b = 1:nBoot

    % Resample subjects with replacement
    boot_idx = randsample(nSub, nSub, true);

    data_b_orig = data_orig(boot_idx, :, :);

    % Re-run ComBat within each bootstrap sample
    T_b = reshape(permute(data_b_orig, [3 2 1]), nTracts, []);

    batch_b = repmat(batch_ses, 1, nSub);

    T_b_cb = combat(T_b, batch_b, [], 1);

    data_b_cb = permute(reshape(T_b_cb, nTracts, nSes, nSub), [3 2 1]);

    [mdc_o, mdc_c] = compute_mdc95(data_b_orig, data_b_cb, typeICC, alpha, r0);

    MDC_orig_boot(b,:) = mdc_o(:)';
    MDC_cb_boot(b,:)   = mdc_c(:)';
end

%% Percentile bootstrap CI
CI_orig_low  = prctile(MDC_orig_boot, 2.5, 1)';
CI_orig_high = prctile(MDC_orig_boot, 97.5, 1)';

CI_cb_low  = prctile(MDC_cb_boot, 2.5, 1)';
CI_cb_high = prctile(MDC_cb_boot, 97.5, 1)';

%% ===== 5) Save table =====
tract_labels = {
    'L TR','R TR', ...
    'L CST','R CST', ...
    'L CGC','R CGC', ...
    'L CGH','R CGH', ...
    'FMA','FMI', ...
    'L ILF','R ILF', ...
    'L SLF','R SLF', ...
    'L UF','R UF', ...
    'L AF','R AF', ...
    'L OR','R OR'
};

%% ===== 6) Plot: point estimates with 95% bootstrap CI =====
fig = figure('Position',[100 100 650 650], 'Renderer','opengl');
hold on;

y = 1:nTracts;

% Slight vertical offset for original and ComBat
offset = 0.15;
y_orig = y - offset;
y_cb   = y + offset;

colo = [0 0 0];

%% Error bars for original data
for tr = 1:nTracts

    x = MDC_orig(tr);
    xlow  = x - CI_orig_low(tr);
    xhigh = CI_orig_high(tr) - x;

    errorbar(x, y_orig(tr), xlow, xhigh, ...
        'horizontal', ...
        'o', ...
        'Color', colo, ...
        'MarkerFaceColor', colo, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerSize', 6, ...
        'LineWidth', 1.2, ...
        'CapSize', 8);
end

%% Error bars for ComBat data
for tr = 1:nTracts

    x = MDC_cb(tr);
    xlow  = x - CI_cb_low(tr);
    xhigh = CI_cb_high(tr) - x;

    errorbar(x, y_cb(tr), xlow, xhigh, ...
        'horizontal', ...
        'o', ...
        'Color', colo, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', colo, ...
        'MarkerSize', 6, ...
        'LineWidth', 1.2, ...
        'CapSize', 8);
end

%% Dummy handles for legend
h1 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', colo, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7);

h2 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', colo, ...
    'MarkerSize', 7);

%% Mean lines
m_orig = mean(MDC_orig, 'omitnan');
m_cb   = mean(MDC_cb,  'omitnan');

yl = [0.5 nTracts+0.5];

h3 = plot([m_orig m_orig], yl, '--', ...
    'Color', colo, ...
    'LineWidth', 2);

h4 = plot([m_cb m_cb], yl, ':', ...
    'Color', colo, ...
    'LineWidth', 2);

%% Axis settings
yticks(y);
yticklabels(tract_labels);

set(gca, ...
    'YDir', 'reverse', ...
    'FontSize', 14);

xlabel('MDC_{95} (MD)', 'FontSize', 14);
ylabel('Tract', 'FontSize', 14);

xMax = max([CI_orig_high; CI_cb_high], [], 'omitnan');
xlim([0 xMax*1.15]);
ylim(yl);

grid on; box on;

legend([h1 h2 h3 h4], ...
    {'The original data','The data after ComBat', ...
     'Mean of the original data','Mean of the data after ComBat'}, ...
    'Location','eastoutside');
legend boxoff

exportgraphics(fig, ...
    'MDC95_MD_bootstrap_CI.eps', ...
    'ContentType','vector');

end

%% ============================================================
% Local function: compute MDC95
% ============================================================
function [MDC_orig, MDC_cb] = compute_mdc95(data_orig, data_cb, typeICC, alpha, r0)

nTracts = size(data_orig, 3);

MDC_orig = nan(nTracts,1);
MDC_cb   = nan(nTracts,1);

for tr = 1:nTracts

    %% Original data
    v_orig = data_orig(:,:,tr);
    v_orig = v_orig(:);

    SD_orig = std(v_orig, 0, 'omitnan');

    M_orig = squeeze(data_orig(:,:,tr));

    ICC_orig = ICC(M_orig, typeICC, alpha, r0);

    MDC_orig(tr) = 1.96 * sqrt(2) * SD_orig * sqrt(1 - ICC_orig);

    %% ComBat data
    v_cb = data_cb(:,:,tr);
    v_cb = v_cb(:);

    SD_cb = std(v_cb, 0, 'omitnan');

    M_cb = squeeze(data_cb(:,:,tr));

    ICC_cb = ICC(M_cb, typeICC, alpha, r0);

    MDC_cb(tr) = 1.96 * sqrt(2) * SD_cb * sqrt(1 - ICC_cb);
end

end