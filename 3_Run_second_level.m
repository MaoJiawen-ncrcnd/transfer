clear; clc; close all;

%% 参数设置
model_tag      = 'full';
group_filter   = 'stroke'; % 只分析卒中人群，我认为不能合并分析卒中和健康人，可能一个正一个负导致连接不显著，但根据假设改

analysis_tag   = sprintf('%s_%s', model_tag, group_filter);
prob_threshold = 0.95;
force_recompute = false; % true 时强制重跑 PEB/BMR；false 时优先复用已有结果。

% 关闭 SPM 在 PEB/BMR 过程中自动弹出的诊断图窗。
% 最终的显著连边热图仍由本脚本单独生成和保存。
try
    spm_get_defaults('cmdline', true);
catch
end

% 使用脚本所在目录作为项目根目录，避免依赖当前工作目录或示例工程路径。
project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end

% 2_Run_first_level.m 估计后应将 GCM_full.mat 保存到 analyses 目录。
out_dir = fullfile(project_root, 'analyses');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

start_dir = fullfile(project_root, 'group_glm_task');
metadata_file = fullfile(project_root, 'TMS-MRI.csv');
gcm_file = fullfile(out_dir, sprintf('GCM_%s.mat', model_tag));

% 区域顺序必须与 2_Run_first_level.m 中 DCM 的 ROI 顺序一致。
% 后续矩阵采用 B(target, source) 方向：行是目标区，列是来源区。
% ROI 顺序必须和 1_ROI_extract.m和2_Run_first_level.m 的mask_name一致，需要改：
region_labels = {'', '', '', '', ...
                 '', '', '', ''};  

nregions = numel(region_labels);

%% 读取一级估计后的 GCM

if ~exist(gcm_file, 'file')
    error(['Missing first-level estimated GCM: %s\n' ...
           'Run 2_Run_first_level.m until it saves analyses/GCM_%s.mat. ' ...
           'Do not use subject-level DCM_%s.mat files here; they may be unestimated.'], ...
           gcm_file, model_tag, model_tag);
end

S = load(gcm_file);
if ~isfield(S, 'GCM')
    error('Expected variable ''GCM'' in %s.', gcm_file);
end
if isfield(S, 'roi_names') && numel(S.roi_names) == numel(region_labels)
    region_labels = cellstr(S.roi_names(:))';
end
if isfield(S, 'GCM_paths')
    gcm_paths = cellstr(S.GCM_paths(:));
else
    gcm_paths = {};
end

% 兼容两类 GCM 保存方式：
% 1) cell array of DCM structs：直接用于 PEB；
% 2) cell array of DCM file paths：先用 spm_dcm_load 载入 DCM 结构体。
GCM = S.GCM;
if iscell(GCM)
    is_path_cell = cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), GCM);
    if all(is_path_cell(:))
        gcm_paths = cellfun(@char, GCM(:), 'UniformOutput', false);
        GCM = spm_dcm_load(gcm_paths);
    end
elseif isstruct(GCM)
    GCM = num2cell(GCM(:));
else
    error('Variable ''GCM'' must be a cell array of DCM structures or DCM file paths.');
end

GCM = GCM(:);
nsubjects = numel(GCM);
if nsubjects == 0
    error('Loaded GCM is empty: %s.', gcm_file);
end

% PEB 必须使用已经估计过的一阶 DCM，因此这里检查 Ep/Cp/F。
% 被试目录下的 DCM_full.mat 可能只是指定好的模型，不能直接用于二阶。
validate_estimated_gcm(GCM, gcm_file);

%% 根据 TMS-MRI.csv 第二列筛选 stroke 被试
% CSV 第 1 列是 ID，如 TMS-002；第 2 列是分组，如 stroke/health。
% DCM 文件夹名通常是 sub-TMS002ses01，因此这里会把两边统一成数字 ID 后匹配。
[GCM, subject_ids, subject_groups, keep_subject] = filter_gcm_by_group( ...
    GCM, group_filter, metadata_file, start_dir, model_tag, gcm_paths);

nsubjects = numel(GCM);
if nsubjects == 0
    error('No subjects remained after filtering group "%s".', group_filter);
end

fprintf('\nSecond-level group filter: %s\n', group_filter);
fprintf('Subjects kept: %d / %d\n', nsubjects, numel(keep_subject));

%% 构建截距项二阶设计矩阵
% 当前只估计组平均连接效应，不加入临床或行为协变量。

X        = ones(nsubjects, 1);
X_labels = {'Group mean'};

M = struct();
M.Q      = 'all';
M.X      = X;
M.Xnames = X_labels;
M.maxit  = 256;

peb_file = fullfile(out_dir, sprintf('PEB_B_%s.mat', analysis_tag));
bma_file = fullfile(out_dir, sprintf('BMA_search_B_%s.mat', analysis_tag));

%% 对任务调制 B 参数建立 PEB

if ~force_recompute && exist(peb_file, 'file') && exist(bma_file, 'file')
    fprintf('\nLoading existing second-level PEB/BMA results for %s.\n', analysis_tag);
    fprintf('Set force_recompute = true to rerun PEB/BMR from GCM.\n');

    peb_data = load(peb_file, 'PEB_B', 'RCM_B');
    bma_data = load(bma_file, 'BMA_B');

    if ~isfield(peb_data, 'PEB_B') || ~isfield(bma_data, 'BMA_B')
        error('Existing PEB/BMA files are incomplete. Delete them or set force_recompute = true.');
    end

    PEB_B = peb_data.PEB_B;
    RCM_B = peb_data.RCM_B;
    BMA_B = bma_data.BMA_B;
else
    fprintf('\nRunning second-level PEB on %s as the full model.\n', analysis_tag);
    fprintf('Subjects: %d\n', nsubjects);

    [PEB_B, RCM_B] = spm_dcm_peb(GCM, M, {'B'});

    save(peb_file, 'PEB_B', 'RCM_B', 'M', 'X', 'X_labels', ...
         'model_tag', 'analysis_tag', 'group_filter', 'subject_ids', ...
         'subject_groups', 'keep_subject', 'gcm_file', 'gcm_paths', ...
         'metadata_file', '-v7.3');

    %% 自动 Bayesian model reduction / Bayesian model averaging
    % BMR/BMA 会在 full model 的参数空间中自动筛选有证据的连接。

    BMA_B = spm_dcm_peb_bmc(PEB_B);

    save(bma_file, 'BMA_B', 'prob_threshold', 'model_tag', 'analysis_tag', ...
         'group_filter', 'subject_ids', 'subject_groups', 'keep_subject', ...
         'gcm_file', 'gcm_paths', 'metadata_file', '-v7.3');
end

%% 生成显著 B 连边矩阵和热图
% 显著性标准：posterior probability Pp > prob_threshold。
% B_Ep_matrix 保存未阈值化效应量；B_Pp_matrix 保存每条边的后验概率。
% B_sig_matrix 只保留显著连边的效应量，不显著连边置 0。

[B_Ep_matrix, B_Pp_matrix, B_sig_matrix, B_sig_mask, ...
    B_param_names, B_param_indices] = build_b_matrices(BMA_B, nregions, prob_threshold, GCM{1}.M.pE);

matrix_file = fullfile(out_dir, sprintf('B_sig_matrix_%s.mat', analysis_tag));
save(matrix_file, 'B_sig_matrix', 'B_sig_mask', 'B_Ep_matrix', 'B_Pp_matrix', ...
     'B_param_names', 'B_param_indices', 'region_labels', 'prob_threshold', ...
     'model_tag', 'analysis_tag', 'group_filter', 'subject_ids', 'subject_groups', ...
     'keep_subject', 'gcm_file', 'gcm_paths', 'metadata_file', 'X', 'X_labels');

png_file = fullfile(out_dir, sprintf('B_sig_matrix_%s.png', analysis_tag));
fig_file = fullfile(out_dir, sprintf('B_sig_matrix_%s.fig', analysis_tag));
plot_sig_b_matrix(B_sig_matrix, region_labels, prob_threshold, analysis_tag, png_file, fig_file);

nsig = nnz(B_sig_mask);
if nsig == 0
    fprintf('No B connections passed Pp > %.2f. Saved an all-zero heatmap.\n', prob_threshold);
else
    fprintf('Significant B connections with Pp > %.2f: %d\n', prob_threshold, nsig);
end
print_top_b_connections(B_Ep_matrix, B_Pp_matrix, region_labels, min(10, numel(B_param_names)));

fprintf('\nDone.\n');
fprintf('PEB saved: %s\n', peb_file);
fprintf('BMA saved: %s\n', bma_file);
fprintf('Matrix saved: %s\n', matrix_file);
fprintf('Heatmap saved: %s\n', png_file);

%% 局部函数

function validate_estimated_gcm(GCM, gcm_file)
    % 检查 GCM 中每个 DCM 是否已经完成一阶估计。
    % Ep/Cp/F 分别对应参数后验均值、后验协方差和模型自由能。
    required_fields = {'Ep', 'Cp', 'F'};

    for s = 1:numel(GCM)
        if ~isstruct(GCM{s})
            error('GCM{%d} in %s is not a DCM structure.', s, gcm_file);
        end

        missing = required_fields(~cellfun(@(f) isfield(GCM{s}, f), required_fields));
        if ~isempty(missing)
            error(['GCM{%d} in %s is missing estimated DCM field(s): %s.\n' ...
                   'Run 2_Run_first_level.m until spm_dcm_fit has completed and ' ...
                   'analyses/GCM_full.mat contains estimated DCMs.'], ...
                   s, gcm_file, strjoin(missing, ', '));
        end
    end
end

function [GCM_filtered, subject_ids, subject_groups, keep_subject] = filter_gcm_by_group( ...
        GCM, group_filter, metadata_file, start_dir, model_tag, gcm_paths)
    % 根据 TMS-MRI.csv 的第 2 列分组筛选 GCM。
    % 优先使用一级脚本保存的 GCM_paths 重建 ID；如果没有，则从 DCM_full.mat 路径重建。

    if ~exist(metadata_file, 'file')
        error('Missing metadata CSV: %s.', metadata_file);
    end

    [metadata_ids, metadata_groups] = read_subject_group_csv(metadata_file);
    metadata_norm_ids = normalize_tms_ids(metadata_ids);

    if numel(unique(metadata_norm_ids)) ~= numel(metadata_norm_ids)
        error('Duplicate subject IDs were found in %s.', metadata_file);
    end

    metadata_groups = cellfun(@lower, strtrim(metadata_groups), 'UniformOutput', false);
    group_map = containers.Map(metadata_norm_ids, metadata_groups);

    all_subject_ids = infer_gcm_subject_ids(numel(GCM), gcm_paths, start_dir, model_tag);
    all_norm_ids = normalize_tms_ids(all_subject_ids);

    all_groups = cell(size(all_norm_ids));
    missing = false(size(all_norm_ids));
    for i = 1:numel(all_norm_ids)
        if isKey(group_map, all_norm_ids{i})
            all_groups{i} = group_map(all_norm_ids{i});
        else
            all_groups{i} = '';
            missing(i) = true;
        end
    end

    if any(missing)
        missing_ids = strjoin(all_subject_ids(missing), ', ');
        error('These GCM subjects were not found in %s: %s.', metadata_file, missing_ids);
    end

    keep_subject = strcmpi(all_groups, group_filter);
    GCM_filtered = GCM(keep_subject);
    subject_ids = all_subject_ids(keep_subject);
    subject_groups = all_groups(keep_subject);
end

function [subject_ids, subject_groups] = read_subject_group_csv(metadata_file)
    % 只读取 CSV 的前两列：第 1 列被试 ID，第 2 列分组。
    % 文件中中文列名可能是 GBK/ISO-8859-1 显示乱码，但 stroke/health 值是 ASCII。
    fid = fopen(metadata_file, 'r', 'n', 'ISO-8859-1');
    if fid < 0
        error('Cannot open metadata CSV: %s.', metadata_file);
    end

    cleaner = onCleanup(@() fclose(fid));
    C = textscan(fid, '%s%s%*[^\n]', 'Delimiter', ',', 'HeaderLines', 1, ...
                 'ReturnOnError', false);

    subject_ids = strtrim(C{1});
    subject_groups = strtrim(C{2});

    valid = ~cellfun(@isempty, subject_ids) & ~cellfun(@isempty, subject_groups);
    subject_ids = subject_ids(valid);
    subject_groups = subject_groups(valid);
end

function subject_ids = infer_gcm_subject_ids(nsubjects, gcm_paths, start_dir, model_tag)
    % 优先使用 2_Run_first_level.m 保存到 GCM_full.mat 中的 GCM_paths。
    % 如果旧结果中没有 GCM_paths，则用 spm_select/dir 重新生成 DCM 文件顺序。

    if ~isempty(gcm_paths)
        if numel(gcm_paths) ~= nsubjects
            error(['GCM_paths contains %d paths, but GCM contains %d subjects. ' ...
                   'Rerun 2_Run_first_level.m or check GCM_full.mat.'], ...
                   numel(gcm_paths), nsubjects);
        end
        subject_ids = extract_subject_ids_from_text(gcm_paths);
        return
    end

    pattern = sprintf('DCM_%s.mat', model_tag);
    paths = {};

    if exist(start_dir, 'dir')
        try
            selected = spm_select('FPListRec', start_dir, pattern);
            paths = cellstr(selected);
        catch
            paths = {};
        end

        if isempty(paths)
            d = dir(fullfile(start_dir, 'sub-*', pattern));
            paths = fullfile({d.folder}', {d.name}');
            paths = sort(paths);
        end
    end

    paths = paths(~cellfun(@isempty, paths));
    if numel(paths) ~= nsubjects
        error(['Cannot infer subject IDs for GCM. Found %d DCM_%s files, ' ...
               'but GCM contains %d subjects.'], numel(paths), model_tag, nsubjects);
    end

    subject_ids = extract_subject_ids_from_text(paths);
end

function subject_ids = extract_subject_ids_from_text(values)
    % 从路径或字符串中提取 TMS 编号，并统一保存为 TMS-XXX。
    if ischar(values)
        values = cellstr(values);
    end

    subject_ids = cell(size(values));
    for i = 1:numel(values)
        txt = char(values{i});
        tok = regexpi(txt, 'TMS[-_]?(\d+)', 'tokens', 'once');
        if isempty(tok)
            error('Cannot extract TMS subject ID from: %s.', txt);
        end
        subject_ids{i} = sprintf('TMS-%03d', str2double(tok{1}));
    end
    subject_ids = subject_ids(:);
end

function norm_ids = normalize_tms_ids(subject_ids)
    % 把 TMS-002、TMS002、sub-TMS002ses01 等格式统一成 002。
    if ischar(subject_ids)
        subject_ids = cellstr(subject_ids);
    end

    norm_ids = cell(size(subject_ids));
    for i = 1:numel(subject_ids)
        txt = char(subject_ids{i});
        tok = regexpi(txt, 'TMS[-_]?(\d+)', 'tokens', 'once');
        if isempty(tok)
            tok = regexp(txt, '(\d+)', 'tokens', 'once');
        end
        if isempty(tok)
            error('Cannot normalize TMS subject ID from: %s.', txt);
        end
        norm_ids{i} = sprintf('%03d', str2double(tok{1}));
    end
    norm_ids = norm_ids(:);
end

function [B_Ep_matrix, B_Pp_matrix, B_sig_matrix, B_sig_mask, ...
          B_param_names, B_param_indices] = build_b_matrices(BMA_B, nregions, prob_threshold, template_pE)
    % 将 BMA_B 中的一维参数向量映射回 8x8 的 B 连接矩阵。
    % 只解析第 1 个 B 条件，其中 i 是 target，j 是 source。
    % SPM 不同版本/模型会把 B 参数命名为 B(i,j,1)、B{1}(i,j) 或 B(i,j)。

    if ~isfield(BMA_B, 'Pnames') || ~isfield(BMA_B, 'Ep') || ~isfield(BMA_B, 'Pp')
        error('BMA_B must contain Pnames, Ep, and Pp fields.');
    end

    pnames = BMA_B.Pnames;
    if ischar(pnames)
        pnames = cellstr(pnames);
    end
    pnames = pnames(:);

    % 如果 Pnames 不是连接名而是 P1/P2/...，尝试用 Pind0 回到 DCM.M.pE 生成真实字段名。
    if ~any(cellfun(@(x) ~isempty(regexp(x, 'B(\(|\{)', 'once')), pnames)) && ...
            nargin >= 4 && isfield(BMA_B, 'Pind0')
        try
            pnames = cell(numel(BMA_B.Pind0), 1);
            for i = 1:numel(BMA_B.Pind0)
                pnames{i} = spm_fieldindices(template_pE, BMA_B.Pind0(i));
            end
        catch
            % 如果 fallback 失败，保留原始 Pnames，后面会给出明确错误。
        end
    end

    Ep = BMA_B.Ep(:);
    Pp = BMA_B.Pp(:);

    nparams = numel(pnames);
    if numel(Ep) < nparams || numel(Pp) < nparams
        error('BMA_B Ep/Pp vectors are shorter than BMA_B.Pnames.');
    end

    % 截距-only 设计下，第一个 Pnames block 对应组平均效应。
    Ep = Ep(1:nparams);
    Pp = Pp(1:nparams);

    B_Ep_matrix  = zeros(nregions, nregions);
    B_Pp_matrix  = nan(nregions, nregions);
    B_sig_matrix = zeros(nregions, nregions);
    B_sig_mask   = false(nregions, nregions);

    B_param_names   = {};
    B_param_indices = zeros(0, 4); % [target, source, condition, BMA 向量索引]

    for p = 1:nparams
        name = strtrim(pnames{p});
        [target, source, condition, ok] = parse_b_parameter_name(name);
        if ~ok
            continue;
        end

        if condition ~= 1
            continue;
        end

        if target < 1 || target > nregions || source < 1 || source > nregions
            error('Parsed out-of-range B parameter name: %s.', name);
        end

        % 保存未阈值化结果，便于之后检查所有候选连接。
        B_Ep_matrix(target, source) = Ep(p);
        B_Pp_matrix(target, source) = Pp(p);

        % 只把 Pp 超过阈值的连边写入显著矩阵。
        if Pp(p) > prob_threshold
            B_sig_matrix(target, source) = Ep(p);
            B_sig_mask(target, source)   = true;
        end

        B_param_names{end + 1, 1} = name; %#ok<AGROW>
        B_param_indices(end + 1, :) = [target, source, condition, p]; %#ok<AGROW>
    end

    if isempty(B_param_names)
        error(['No B condition-1 parameters were found in BMA_B.Pnames.\n' ...
               'Expected names like B(i,j,1), B{1}(i,j), or B(i,j).']);
    end
end

function [target, source, condition, ok] = parse_b_parameter_name(name)
    % 兼容 SPM 常见的 B 参数命名格式：
    %   B(i,j,1)   : 3D numeric array，第三维是条件
    %   B{1}(i,j)  : cell array，第 1 个 cell 是条件
    %   B(i,j)     : 单条件 2D B 矩阵，默认 condition = 1

    target = NaN;
    source = NaN;
    condition = NaN;
    ok = false;

    tok = regexp(name, 'B\((\d+)\s*,\s*(\d+)\s*,\s*(\d+)\)', 'tokens', 'once');
    if ~isempty(tok)
        target    = str2double(tok{1});
        source    = str2double(tok{2});
        condition = str2double(tok{3});
        ok = true;
        return
    end

    tok = regexp(name, 'B\{(\d+)\}\((\d+)\s*,\s*(\d+)\)', 'tokens', 'once');
    if ~isempty(tok)
        condition = str2double(tok{1});
        target    = str2double(tok{2});
        source    = str2double(tok{3});
        ok = true;
        return
    end

    tok = regexp(name, 'B\((\d+)\s*,\s*(\d+)\)', 'tokens', 'once');
    if ~isempty(tok)
        target    = str2double(tok{1});
        source    = str2double(tok{2});
        condition = 1;
        ok = true;
    end
end

function plot_sig_b_matrix(B_sig_matrix, region_labels, prob_threshold, model_tag, png_file, fig_file)
    % 绘制以 0 为中心的发散色图：红色为正效应，蓝色为负效应。
    % 使用 painters + print 保存 PNG，避免 exportgraphics 在部分远程/WSL 图形环境中导出黑屏。
    fig = figure('Name', sprintf('Significant B effects - %s', model_tag), ...
                 'NumberTitle', 'off', 'Color', 'w', 'Visible', 'off', ...
                 'Renderer', 'painters', 'InvertHardcopy', 'off');
    ax = axes('Parent', fig, 'Color', 'w');

    imagesc(ax, B_sig_matrix);
    axis(ax, 'square');
    colormap(ax, diverging_colormap(256));

    cmax = max(abs(B_sig_matrix(:)));
    if isempty(cmax) || isnan(cmax) || cmax == 0
        cmax = 1;
    end
    caxis(ax, [-cmax cmax]);

    cb = colorbar(ax);
    ylabel(cb, 'BMA Ep');

    nregions = numel(region_labels);
    set(ax, 'XTick', 1:nregions, 'XTickLabel', region_labels, ...
             'YTick', 1:nregions, 'YTickLabel', region_labels, ...
             'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none', ...
             'FontSize', 10, 'Color', 'w', 'Box', 'on', ...
             'XColor', 'k', 'YColor', 'k');

    xlabel(ax, 'Source region');
    ylabel(ax, 'Target region');
    if nnz(B_sig_matrix) == 0
        title(ax, sprintf('%s full model: no B effects with Pp > %.2f', model_tag, prob_threshold), ...
              'Interpreter', 'none');
    else
        title(ax, sprintf('%s full model: B effects with Pp > %.2f', model_tag, prob_threshold), ...
              'Interpreter', 'none');
    end

    hold(ax, 'on');
    for k = 0.5:1:(nregions + 0.5)
        line(ax, [0.5, nregions + 0.5], [k, k], 'Color', [0.85 0.85 0.85]);
        line(ax, [k, k], [0.5, nregions + 0.5], 'Color', [0.85 0.85 0.85]);
    end
    hold(ax, 'off');

    set(fig, 'Position', [100, 100, 900, 780]);
    drawnow;
    print(fig, png_file, '-dpng', '-r300');
    savefig(fig, fig_file);
    close(fig);
end

function print_top_b_connections(B_Ep_matrix, B_Pp_matrix, region_labels, ntop)
    % 打印后验概率最高的 B 连边，帮助判断“无显著”是完全没有证据，
    % 还是最高 Pp 接近但未达到 0.95。
    valid = ~isnan(B_Pp_matrix);
    if ~any(valid(:))
        fprintf('No valid B Pp values were mapped to the matrix.\n');
        return
    end

    target_idx = repmat((1:size(B_Pp_matrix, 1))', 1, size(B_Pp_matrix, 2));
    source_idx = repmat(1:size(B_Pp_matrix, 2), size(B_Pp_matrix, 1), 1);

    target_idx = target_idx(valid);
    source_idx = source_idx(valid);
    pp_values  = B_Pp_matrix(valid);
    ep_values  = B_Ep_matrix(valid);

    [pp_values, order] = sort(pp_values(:), 'descend');
    ep_values  = ep_values(order);
    target_idx = target_idx(order);
    source_idx = source_idx(order);

    ntop = min(ntop, numel(pp_values));
    fprintf('\nTop %d B connections by posterior probability:\n', ntop);
    for i = 1:ntop
        fprintf('  %2d. %s <- %s: Pp = %.4f, Ep = %.4f\n', ...
                i, region_labels{target_idx(i)}, region_labels{source_idx(i)}, ...
                pp_values(i), ep_values(i));
    end
end

function cmap = diverging_colormap(n)
    % 简单自定义蓝-白-红色图，避免依赖额外工具箱。
    if nargin < 1
        n = 256;
    end

    blue  = [0.230, 0.299, 0.754];
    white = [1.000, 1.000, 1.000];
    red   = [0.706, 0.016, 0.150];

    x = linspace(0, 1, n)';
    cmap = zeros(n, 3);
    for c = 1:3
        cmap(:, c) = interp1([0, 0.5, 1], [blue(c), white(c), red(c)], x);
    end
end
