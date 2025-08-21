%% ------------------------------------------------------------------------
% Script Name : run_sted_batch.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   Batch driver script for processing multiple .msr files with
%   process_sted_file(). This script:
%     1) Recursively discovers all *.msr files under a given root folder.
%     2) Filters files to those containing "sted" (case-insensitive)
%        AND "pcExc" in the filename.
%     3) Uses Bio-Formats (bfopen) to probe each file and determine
%        whether it is a single image (≤2 planes) or a stack (≥3 planes).
%     4) Prints a summary of discovered datasets (counts of singles/stacks).
%     5) Loops over all accepted files and calls process_sted_file()
%        to carry out segmentation and profile extraction.
%
%   The output of each call is written to the directory specified by
%   `outDir`, preserving per-file subdirectories. Figures and segmentation
%   .mat files are produced by the per-file function.
%
% Usage :
%   % Set the root folder where *.msr files are located
%   stedRoot = 'C:\path\to\data';
%
%   % Run the batch script
%   run_sted_batch
%
%   % Notes:
%   %   - Adjust 'outDir' at the top of this script if desired.
%   %   - Calls process_sted_file(filePath, outDir) for each match.
%
% Dependencies :
%   - MATLAB R2021b+ (tested with imrotate3/bwskel/fit/denoiseImage)
%   - Image Processing Toolbox
%   - Curve Fitting Toolbox
%   - Deep Learning Toolbox
%   - Bio-Formats for MATLAB (bfopen)
%   - process_sted_file.m (this project)
%
% Reference :
%   - process_sted_file.m header for algorithmic details
%   - Bio-Formats documentation: openmicroscopy.org/bio-formats/
%
% License :
%   MIT
%   Copyright (c) 2025 Felix Maurer
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%     The above copyright notice and this permission notice shall be included in
%     all copies or substantial portions of the Software.
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%% ------------------------------------------------------------------------
% Output file system (matches original default)
% -------------------------------------------------------------------------
outDir = 'SegmentationResults3D_deconv_separateAdj';
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------------------------------------------------------
% Files: discover candidate .msr files (recursive)
% -------------------------------------------------------------------------
clc
stedRoot = '';
files = dir(fullfile([stedRoot, filesep, '**', filesep, '*.msr']));

% Preallocate/holders to mirror original behavior
clear sizes; k_stedFiles = 1;
clear stedData;

%% ------------------------------------------------------------------------
% Pre-scan: load headers to classify into single images vs stacks
% (stores data/file descriptors; filters by name condition)
% -------------------------------------------------------------------------
for fileIdx = 1:length(files)
    fileName   = files(fileIdx).name;
    fileFolder = files(fileIdx).folder;
    filePath   = [fileFolder, filesep, fileName];

    % Selection rule: must contain "sted" (case-insensitive) AND "pcExc"
    condition = contains(lower(fileName),'sted') & contains(fileName,'pcExc');

    if condition
        fprintf('%s\n',filePath);

        % Load to measure number of planes (single vs stack) like original
        data = bfopen(filePath);

        sizes(k_stedFiles)            = length(data{1});
        stedData(k_stedFiles).data    = data;
        stedData(k_stedFiles).filePath= filePath;
        stedData(k_stedFiles).fileName= fileName;

        k_stedFiles = k_stedFiles + 1;
    end
end

%% ------------------------------------------------------------------------
% Summary: report how many files and how many singles vs stacks
% -------------------------------------------------------------------------
clc
fprintf('%d files with sted data found.\n', length(stedData));
fprintf('%d single, %d stacks found.\n', sum(sizes < 3), sum(sizes >= 3));
fprintf('--------------- analyzing stacks -----------------\n');

%% ------------------------------------------------------------------------
% Process: loop over selected files and call the per-file function
% -------------------------------------------------------------------------
for idx = 1:length(sizes)
    close all;  % keep parity with original workflow
    fprintf('%s\n', stedData(idx).filePath);
    process_sted_file(stedData(idx).filePath, outDir);
end
