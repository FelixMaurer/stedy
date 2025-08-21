%% demo_run.m
% Minimal end-to-end demo for the STED membrane profiling pipeline.
% - Adds repo code + Bio-Formats to MATLAB path
% - Processes example/data/sample_sted_pcExc.msr
% - Writes outputs into example/output
% - Verifies that segmentation files were produced

clc; clear; close all;

%% ------------------------------------------------------------------------
% Resolve paths
% -------------------------------------------------------------------------
thisFile = mfilename('fullpath');
exampleDir = fileparts(thisFile);
repoRoot   = fileparts(exampleDir);
exampleDir = [exampleDir,'\example'];

codeFiles = { ...
    fullfile(repoRoot,'process_sted_file.m'), ...
    fullfile(repoRoot,'traceSkel.m'), ...
    fullfile(repoRoot,'run_sted_batch.m') ...
};

% Check required code files exist
missing = codeFiles(~cellfun(@exist, codeFiles, repmat({"file"},size(codeFiles))));
if ~isempty(missing)
    warning('Some core files are missing from the repo:\n%s', strjoin(missing, newline));
end

% Add repo root to path
addpath(repoRoot);

% Add Bio-Formats (bfmatlab) to path if user set BF_PATH env var or common locations
bfCandidates = {
    getenv('BF_PATH') ...
    fullfile(repoRoot,'third_party','bfmatlab') ...
    fullfile(userpath, 'bfmatlab') ...
};
bfAdded = false;
for k = 1:numel(bfCandidates)
    if ~isempty(bfCandidates{k}) && exist(bfCandidates{k},'dir')
        addpath(genpath(bfCandidates{k}));
        bfAdded = true;
        break;
    end
end

% Final check for bfopen
if ~exist('bfopen','file')
    fprintf(2, ['[ERROR] Bio-Formats (bfmatlab) not found on the MATLAB path.\n' ...
        'Please download and add it, e.g.:\n' ...
        '  addpath(genpath(''C:\\path\\to\\bfmatlab''));\n\n' ...
        'Download: https://www.openmicroscopy.org/bio-formats/\n\n']);
    % Continue anyway—process_sted_file will fail gracefully if missing.
end

%% ------------------------------------------------------------------------
% Input example data and output directory
% -------------------------------------------------------------------------
msrPath = fullfile(exampleDir,'data','sample_sted_pcExc.msr');
obfPath = fullfile(exampleDir,'data','sample_sted_pcExc_deconv.obf');
outDir  = fullfile(exampleDir,'output');

if ~exist(msrPath,'file')
    fprintf(2, ['[ERROR] Example .msr not found:\n  %s\n\n' ...
        'Please place a small test file there or adjust msrPath in demo_run.m.\n'], msrPath);
    return;
end

% Create output folder
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------------------------------------------------------
% Run per-file processing
% -------------------------------------------------------------------------
fprintf('[INFO] Processing example file:\n  %s\n', msrPath);
try
    process_sted_file(msrPath, outDir);
catch ME
    fprintf(2, '[ERROR] process_sted_file failed: %s\n', ME.message);
    rethrow(ME);
end

%% ------------------------------------------------------------------------
% Verify outputs
% -------------------------------------------------------------------------
segFiles = dir(fullfile(outDir,'**','*segmentation.mat'));
pngFiles = dir(fullfile(outDir,'**','*.png'));

fprintf('[INFO] Found %d segmentation .mat file(s) and %d PNG figure(s) under:\n  %s\n', ...
    numel(segFiles), numel(pngFiles), outDir);

if isempty(segFiles)
    fprintf(2, ['[WARN] No segmentation files found. Check that your example .msr\n' ...
        '       contains the expected channels and that a matching *_deconv*.obf\n' ...
        '       exists near the .msr (discovered recursively).\n']);
else
    % Peek into the first segmentation file
    segPath = fullfile(segFiles(1).folder, segFiles(1).name);
    S = load(segPath);
    if isfield(S,'segmentation')
        seg = S.segmentation;
        nProfiles = numel(seg.profiles);
        szConf = size(seg.confImg);
        szSted = size(seg.stedImg);
        fprintf('[INFO] Example segmentation summary:\n');
        fprintf('  - Profiles traced: %d\n', nProfiles);
        fprintf('  - confImg size:    %dx%d\n', szConf(1), szConf(2));
        fprintf('  - stedImg size:    %dx%d\n', szSted(1), szSted(2));
    else
        fprintf(2, '[WARN] segmentation struct not found in %s\n', segPath);
    end
end

%% ------------------------------------------------------------------------
% Optional: show one of the generated figures if present
% -------------------------------------------------------------------------
if ~isempty(pngFiles)
    try
        % Prefer a *_mem_unwrap.png for quick visual inspection
        unwrapIdx = find(contains({pngFiles.name}, 'mem_unwrap'), 1);
        if isempty(unwrapIdx), unwrapIdx = 1; end
        imgPath = fullfile(pngFiles(unwrapIdx).folder, pngFiles(unwrapIdx).name);
        I = imread(imgPath);
        figure('Name','Demo Output Preview'); imshow(I,[]);
        title(sprintf('Preview: %s', pngFiles(unwrapIdx).name), 'Interpreter','none');
    catch ME
        fprintf(2, '[WARN] Could not preview PNG: %s\n', ME.message);
    end
end

fprintf('[DONE] Demo finished.\n');
