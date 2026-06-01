clear;clc;close all;
spm('Defaults','fMRI');
spm_jobman('initcfg');
%%
workdir = '/data02/linzhenzong/TMS-TASK/taskfmri'; % 不用改
outdir = '/data02/linzhenzong/TMS-TASK_2'; % 输出路径，需要改
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% 定义本轮分析的mask，务必改
mask1 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask1.nii,1';
mask2 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask2.nii,1';
mask3 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask3.nii,1';
mask4 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask4.nii,1';
mask5 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask5.nii,1';
mask6 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask6.nii,1';
mask7 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask7.nii,1';
mask8 = '/data02/linzhenzong/TMS-TASK_2/ROI/mask8.nii,1';

% 定义本轮分析的mask后续的名称，务必改
mask1_name = '';
mask2_name = '';
mask3_name = '';
mask4_name = '';
mask5_name = '';
mask6_name = '';
mask7_name = '';
mask8_name = '';

task_cond_name = 'tap_lesion'; % 定义本轮分析的任务条件名称，根据情况改，可选tap_health

% 下面的参数不改
TR = 0.4;
dropN = 0;

group_outdir = fullfile(outdir, 'group_glm_task');
if ~exist(group_outdir, 'dir')
    mkdir(group_outdir);
end

event_template = fullfile(workdir, 'task-finger_events.tsv');
meta_csv = fullfile(workdir, 'TMS-MRI.csv');


%%
info = readtable(meta_csv);

info.ID_norm = upper(strrep(strtrim(string(info.ID)), '-', ''));
info.lesion_norm = upper(strtrim(string(info.('lesion'))));
lesion_map = containers.Map();
for i = 1:height(info)
    if strlength(info.ID_norm(i)) > 0
        lesion_map(char(info.ID_norm(i))) = char(info.lesion_norm(i));
    end
end
E_template = readtable(event_template, 'FileType', 'text', 'Delimiter', '\t');

name_subj = dir(fullfile(workdir, 'sourcedata', 'sub-TMS002*'));
name_subj = name_subj([name_subj.isdir]);

for s0 = 1:length(name_subj)
    subj_dir = fullfile(workdir, 'sourcedata', name_subj(s0).name);

    %% -------------------------
    % 1) 解析被试ID，用于匹配 TMS-MRI.csv
    % 例：sub-TMS001ses01 -> TMS001
    % --------------------------
    subj_name = name_subj(s0).name;
    tok = regexp(subj_name, '^sub-?([A-Za-z]+[0-9]+)', 'tokens', 'once');

    subj_id_norm = upper(strrep(strtrim(string(tok{1})), '-', ''));
    lesion_side = lesion_map(char(subj_id_norm));

    E = E_template;

    tt = string(E.trial_type);

    if strcmpi(lesion_side, 'L')
        % 右手是患手
        tt(tt == "tap_right") = "tap_lesion";
        tt(tt == "tap_left")  = "tap_health";
    else 
        % 左手是患手
        tt(tt == "tap_left")  = "tap_lesion";
        tt(tt == "tap_right") = "tap_health";
    end

    E.trial_type = tt;

    % 可选：把重编码后的事件表另存一份，方便核查
    out_event = fullfile(subj_dir, 'task-finger_events_relabeled.tsv');
    %writetable(E, out_event, 'FileType', 'text', 'Delimiter', '\t');
end

%%
name_subj = dir(fullfile(workdir, 'sourcedata', 'sub*'));
name_subj = name_subj([name_subj.isdir]);

for s0 = 1 : length(name_subj)

    subj_dir = fullfile(workdir, 'sourcedata', name_subj(s0).name);
    subj_outdir = fullfile(outdir, 'sourcedata', name_subj(s0).name);
    if ~exist(subj_outdir, 'dir')
        mkdir(subj_outdir);
    end

    f = spm_select('ExtFPList', subj_dir, '^.*MNI152NLin6Asym_res-02_desc-preproc_bold\.nii$', Inf);
    
    if dropN > 0
        if size(f,1) <= dropN
            warning('scan 数不足，跳过: %s', name_subj(s0).name);
            continue;
        end
        f = f(dropN+1:end, :);
    end
    
    nscan = size(f,1);

    conf = dir(fullfile(subj_dir, '*finger*_desc-confounds_timeseries.tsv'));
    T = readtable(fullfile(subj_dir, conf(1).name), 'FileType','text', 'Delimiter','\t');

    X = [];
    motion_cols = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    for k = 1:numel(motion_cols)
        if ismember(motion_cols{k}, T.Properties.VariableNames)
            x = T.(motion_cols{k});
            x(isnan(x)) = 0;
            X = [X, x];
        end
    end

    % 所有 motion_outlier* 列
    all_cols = T.Properties.VariableNames;
    outlier_idx = startsWith(all_cols, 'motion_outlier');
    outlier_cols = all_cols(outlier_idx);
    
    for k = 1:numel(outlier_cols)
        x = T.(outlier_cols{k});
        x(isnan(x)) = 0;
        X = [X, x];
    end
 
    % 与扫描数对齐
    if dropN > 0
        if size(X,1) > dropN
            X = X(dropN+1:end, :);
        else
            warning('confounds 行数不足，跳过: %s', name_subj(s0).name);
            continue;
        end
    end

    if size(X,1) > nscan
        X = X(1:nscan, :);
    elseif size(X,1) < nscan
        X = [X; zeros(nscan-size(X,1), size(X,2))];
    end

    rp_path = fullfile(subj_dir, 'multi_reg_task.txt');
    %writematrix(X, rp_path, 'Delimiter', '\t');

    evfile = fullfile(subj_dir, 'task-finger_events_relabeled.tsv');
    E = readtable(evfile, 'FileType', 'text', 'Delimiter', '\t');
    
    trial_types = string(E.trial_type);
    trial_types = strtrim(trial_types);
    
    % 去掉空值
    valid_idx = ~ismissing(trial_types) & trial_types ~= "";
    E = E(valid_idx, :);
    trial_types = trial_types(valid_idx);
    
    % 所有任务都显性建模
    cond_names = unique(trial_types, 'stable');
    nCond = numel(cond_names);

    glmdir = fullfile(subj_outdir, 'glm_task');
    if ~exist(glmdir, 'dir')
        mkdir(glmdir);
    end

    clear matlabbatch
	close all

    
    % First GLM specification
    %-----------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.dir = {glmdir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(f);

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct([]);
    
    for c = 1:nCond
        cond_idx = strcmp(trial_types, cond_names(c));
    
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name = char(cond_names(c));
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset = E.onset(cond_idx);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = E.duration(cond_idx);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth = 1;
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rp_path};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % First GLM estimation
    %-----------------------------------------------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(glmdir,'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % Contrast manager
    % contrast 1: tapping > baseline
    % contrast 2: Effects of interest
    % --------------------------
    matlabbatch{3}.spm.stats.con.spmmat = {fullfile(glmdir, 'SPM.mat')};
    
    lesion_idx = find(cond_names == "tap_lesion", 1);   
    % tap_lesion > implicit baseline
    w = zeros(1, nCond);
    w(lesion_idx) = 1;
    
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'tap_lesion > baseline';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = w;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    % F contrast 
    matlabbatch{3}.spm.stats.con.consess{2}.fcon.name = 'Effects of interest';
    matlabbatch{3}.spm.stats.con.consess{2}.fcon.weights = eye(nCond);
    matlabbatch{3}.spm.stats.con.consess{2}.fcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.delete = 0;
    %% -------------------------
    % 提取 ROI VOI
    % --------------------------
    % ROI 1
    matlabbatch{4}.spm.util.voi.spmmat = {fullfile(glmdir, 'SPM.mat')};
    matlabbatch{4}.spm.util.voi.adjust = 2;
    matlabbatch{4}.spm.util.voi.session = 1;
    matlabbatch{4}.spm.util.voi.name = mask1_name;

    matlabbatch{4}.spm.util.voi.roi{1}.mask.image = {mask1};
    matlabbatch{4}.spm.util.voi.roi{1}.mask.threshold = 0.5;

    matlabbatch{4}.spm.util.voi.roi{2}.spm.spmmat = {fullfile(glmdir, 'SPM.mat')};
    matlabbatch{4}.spm.util.voi.roi{2}.spm.contrast = 1;      % tapping > baseline
    matlabbatch{4}.spm.util.voi.roi{2}.spm.conjunction = 1;
    matlabbatch{4}.spm.util.voi.roi{2}.spm.threshdesc = 'none';
    matlabbatch{4}.spm.util.voi.roi{2}.spm.thresh = 0.001;
    matlabbatch{4}.spm.util.voi.roi{2}.spm.extent = 0;

    matlabbatch{4}.spm.util.voi.roi{3}.mask.image = {fullfile(glmdir,'mask.nii')};
    matlabbatch{4}.spm.util.voi.expression = 'i1';

    % ROI 2
    matlabbatch{5} = matlabbatch{4};
    matlabbatch{5}.spm.util.voi.name = mask2_name;
    matlabbatch{5}.spm.util.voi.roi{1}.mask.image = {mask2};

    % ROI 3
    matlabbatch{6} = matlabbatch{4};
    matlabbatch{6}.spm.util.voi.name = mask3_name;
    matlabbatch{6}.spm.util.voi.roi{1}.mask.image = {mask3};

    % ROI 4
    matlabbatch{7} = matlabbatch{4};
    matlabbatch{7}.spm.util.voi.name = mask4_name;
    matlabbatch{7}.spm.util.voi.roi{1}.mask.image = {mask4};

    % ROI 5
    matlabbatch{8} = matlabbatch{4};
    matlabbatch{8}.spm.util.voi.name = mask5_name;
    matlabbatch{8}.spm.util.voi.roi{1}.mask.image = {mask5};

    % ROI 6
    matlabbatch{9} = matlabbatch{4};
    matlabbatch{9}.spm.util.voi.name = mask6_name;
    matlabbatch{9}.spm.util.voi.roi{1}.mask.image = {mask6};

    % ROI 7
    matlabbatch{10} = matlabbatch{4};
    matlabbatch{10}.spm.util.voi.name = mask7_name;
    matlabbatch{10}.spm.util.voi.roi{1}.mask.image = {mask7};

    % ROI 8
    matlabbatch{11} = matlabbatch{4};
    matlabbatch{11}.spm.util.voi.name = mask8_name;
    matlabbatch{11}.spm.util.voi.roi{1}.mask.image = {mask8};


    %% -------------------------
    % 9) 运行
    % --------------------------
    try
        spm_jobman('run', matlabbatch);
    catch ME
        warning('SPM batch 运行失败: %s | %s', name_subj(s0).name, ME.message);
        continue;
    end

    %% -------------------------
    % 10) 复制结果到 group 目录
    % --------------------------
    group_subj_dir = fullfile(group_outdir, name_subj(s0).name);
    if ~exist(group_subj_dir, 'dir')
        mkdir(group_subj_dir);
    end

    copyfile(fullfile(glmdir, 'SPM.mat'), group_subj_dir);
    copyfile(fullfile(glmdir, 'VOI_*.mat'), group_subj_dir);

    fprintf('完成: %s\n', name_subj(s0).name);
end

fprintf('\n全部完成。\n');

