clear;clc;close all;
%% Settings

% Experiment settings
start_dir = "/data02/linzhenzong/TMS-TASK_2/group_glm_task"; % 这是1的输出路径，务必改
cd(start_dir);

subjdir = dir(start_dir);
subjdir(1:2,:) = [];

% ROI 顺序必须和 1_ROI_extract.m 的mask_name一致，需要改：
roi_names = {'', '', '', '', ...
             '', '', '', ''};

nsubjects   = length(subjdir);
nregions    = numel(roi_names);
nconditions = 1;   % DCM中纳入的任务输入个数；当前只纳入 tap_lesion，看情况改 
model_tag   = 'full';

% Index of each condition in the DCM ，看情况改 
TASK = 1;          % tap_lesion 在DCM输入中的编号。
                   % tap_lesion 在SPM.Sess(1).U中的位置由下面的 tap_idx 动态查找。
                   % 如果要选进多个任务，需要同步修改 nconditions、include、b、c。

% ROI 在 DCM 矩阵中的索引，就是1中的顺序，需要改：
ilM1  = 1; clM1  = 2;
ilSMA = 3; clSMA = 4;
ilPMd = 5; clPMd = 6;
ilPMv = 7; clPMv = 8;

prem6 = [ilSMA ilPMd ilPMv clSMA clPMd clPMv];  % 任务输入 tap_lesion 直接作用于 SMA 和 PMd/v，需要改

% MRI scanner settings
TR = 0.4;    % Repetition time (secs)
TE = 0.028;  % Echo time (secs)

out_dir = '/data02/linzhenzong/TMS-TASK_2/analyses';
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

%% Specify DCMs (one per subject)

% A-matrix (on / off)，这里设置只分析ilM1与其他所有区域之间的连接，其他连接全开。根据假设改。
a = zeros(nregions, nregions);
a(1:nregions+1:end) = 1;
for r = 1:nregions
    if r ~= ilM1
        a(r, ilM1) = 1;   % ilM1 -> r
        a(ilM1, r) = 1;   % r -> ilM1
    end
end

% B-matrix: full model. tap_lesion 调制当前 A 结构中的所有非自连接，一般不用改。
b = zeros(nregions, nregions, nconditions);
full_B = a;
full_B(1:nregions+1:end) = 0;
b(:, :, TASK) = full_B;

figure('Name', 'B full model', 'Color', 'w');
imagesc(full_B);
axis square;
set(gca, 'XTick', 1:nregions, 'XTickLabel', roi_names, ...
         'YTick', 1:nregions, 'YTickLabel', roi_names, ...
         'XTickLabelRotation', 45);
title('Full B model', 'FontWeight', 'bold');
colorbar;

% C-matrix，这是任务输入矩阵，当前设置 tap_lesion 直接作用于 SMA 和 PMd/v，根据假设改，改的是prem6和TASK的对应关系。
c = zeros(nregions,nconditions);
c(prem6,TASK) = 1;

% D-matrix (disabled)
d = zeros(nregions,nregions,0);

for subject = 1:nsubjects
    
    name = subjdir(subject).name;
    
    % Load SPM
    glm_dir = fullfile(start_dir,name);
    SPM     = load(fullfile(glm_dir,'SPM.mat'));
    SPM     = SPM.SPM;
    
    % Load ROIs
    f = cellfun(@(roi) fullfile(glm_dir, sprintf('VOI_%s_1.mat', roi)), ...
                roi_names, 'UniformOutput', false);
    xY_cell = cell(1, nregions);
    for r = 1:nregions
        if ~exist(f{r}, 'file')
            error('Subject %s: missing ROI file %s', name, f{r});
        end
        XY = load(f{r});
        xY_cell{r} = XY.xY;
        xY_cell{r}.name = roi_names{r};
    end
    xY = [xY_cell{:}];
    
    % Move to output directory
    cd(glm_dir);
    
    % 读取一阶 GLM 中所有条件名
    all_cond_names = cell(numel(SPM.Sess(1).U), 1);
    for u = 1:numel(SPM.Sess(1).U)
        all_cond_names{u} = char(SPM.Sess(1).U(u).name{1});
    end

    % 找到 tap_lesion 在一阶条件列表 SPM.Sess(1).U 中的位置，根据任务改tap_lesion
    tap_idx = find(strcmp(all_cond_names, 'tap_lesion'), 1);
    if isempty(tap_idx)
        error('Subject %s: tap_lesion not found in SPM.Sess(1).U', name);
    end

    % 只把 tap_lesion 选进 DCM
    include = zeros(numel(SPM.Sess(1).U), 1);
    include(tap_idx) = 1;

    % Specify. Corresponds to the series of questions in the GUI.
    s = struct();
    s.name       = model_tag;
    s.u          = include;                 % Conditions
    s.delays     = repmat(TR,1,nregions);   % Slice timing for each region
    s.TE         = TE;
    s.nonlinear  = false;
    s.two_state  = false;
    s.stochastic = false;
    s.centre     = true;
    s.induced    = 0;
    s.a          = a;
    s.b          = b;
    s.c          = c;
    s.d          = d;
    DCM = spm_dcm_specify(SPM,xY,s);
    save(fullfile(glm_dir, sprintf('DCM_%s.mat', model_tag)), 'DCM');

    % Return to script directory
    cd(start_dir);
end

%% Collate into a GCM file and estimate

use_parfor = true;

fprintf('\nEstimating GCM for %s model ...\n', model_tag);

pattern = sprintf('DCM_%s.mat', model_tag);
dcms = spm_select('FPListRec', start_dir, pattern);

if isempty(dcms)
    error('No DCM files found for %s', model_tag);
end

GCM_paths = cellstr(dcms);
GCM = spm_dcm_load(GCM_paths);
GCM = spm_dcm_fit(GCM, use_parfor);

save(fullfile(out_dir, sprintf('GCM_%s.mat', model_tag)), ...
     'GCM', 'GCM_paths', 'model_tag', 'roi_names', '-v7.3');

%% Run diagnostics
try
    spm_dcm_fmri_check(GCM);
catch ME
    warning('DCM:DiagnosticFailed', 'spm_dcm_fmri_check failed: %s', ME.message);
end

fprintf('\nDone.\n');
fprintf('Estimated GCMs saved in: %s\n', out_dir);
