function process_sted_file(msrFilePath, outDir)
%% ------------------------------------------------------------------------
% Script Name : process_sted_file.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   Process a single Leica/Abberior .msr file and write analysis outputs
%   (overlay plots + segmentation .mat) into outDir. The pipeline:
%     1) Load confocal & STED channels via Bio‑Formats and the matching
%        deconvolved *.obf (auto‑located).
%     2) Denoise with DnCNN and binarize.
%     3) Morphological cleanup (open → area‑filter → close → thicken).
%     4) Fuse masks: membraneBW = BWConf | BWSted.
%     5) Extract and smooth membranes with cubic smoothing splines.
%     6) Skeletonize the membrane, trace centerlines, and re‑smooth.
%     7) Arc‑length parametrization of centerlines:
%        then resample at fixed Δs = 4 px to get p(s) and n(s).
%     8) Sample intensity profiles along normals
%        from raw and deconvolved (confocal + STED) images.
%     9) Save figures (membrane/normal overlays, unwrapped profiles) and
%        a 'segmentation' struct with profiles, images, and geometry.
%
% Usage :
%   process_sted_file(msrFilePath, outDir)
%   % Example:
%   %   process_sted_file('C:\data\sample.msr', 'SegmentationResults3D_deconv_separateAdj')
%   %
%   % Batch usage (see companion script):
%   %   run_sted_batch.m discovers *.msr and calls this function per file.
%
% Dependencies :
%   - MATLAB R2021b+ (tested with imrotate3/bwskel/fit/denoiseImage)
%   - Image Processing Toolbox (imbinarize, bwmorph, bwboundaries, bwskel, etc.)
%   - Curve Fitting Toolbox (fit with 'smoothingspline')
%   - Deep Learning Toolbox (denoisingNetwork('DnCNN'), denoiseImage)
%   - Bio‑Formats for MATLAB (bfopen) — https://www.openmicroscopy.org/bio-formats/
%   - traceSkel.m — helper that traces skeleton polylines
%   - Matching deconvolved *.obf located recursively next to the .msr
%
% Reference :
%   - Zhang & Zuo et al., "Beyond a Gaussian Denoiser: Residual Learning of
%     Deep CNN for Image Denoising (DnCNN)", IEEE TIP, 2017.
%   - Bio‑Formats documentation: openmicroscopy.org/bio-formats/
%   - Curvature and arc‑length basics in planar curves (any standard text).
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
% Output file system
% -------------------------------------------------------------------------
suffix = 'out';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Save the running script for provenance (best-effort)
try
    copyfile([mfilename('fullpath')],[outDir,filesep,'script',suffix,'.m']);
catch
    % Non-fatal if this fails (e.g., permissions)
end

%% ------------------------------------------------------------------------
% Load input (.msr) and prepare descriptors
% -------------------------------------------------------------------------
clc
data = bfopen(msrFilePath);
sizes = length(data{1});    % for parity with original code

fileFolder = fileparts(msrFilePath);
[~,fileName,ext] = fileparts(msrFilePath);
fileName = [fileName,ext];

%% ------------------------------------------------------------------------
% Select data: split into confocal and STED stacks
% -------------------------------------------------------------------------
confStack = data{1,1};
stedStack = data{2,1};

%% ------------------------------------------------------------------------
% Take deconv: locate matching deconvolved .obf (recursive from msr folder)
% -------------------------------------------------------------------------
decFiles = dir(fullfile(fileFolder,['**',filesep, fileName(1:end-5), '*deconv*.obf']));
decPath = [decFiles(1).folder,filesep,decFiles(1).name];
deconvData = bfopen(decPath);

decConfStack = deconvData{1,1};
decStedStack = deconvData{2,1};

%% ------------------------------------------------------------------------
% Build 3D arrays (confocal/STED raw and deconvolved)
% -------------------------------------------------------------------------
conf3D    = zeros(size(confStack{1,1},1),size(confStack{1,1},2),size(confStack,1),'uint8');
sted3D    = zeros(size(confStack{1,1},1),size(confStack{1,1},2),size(confStack,1),'uint8');
decConf3D = zeros(size(confStack{1,1},1),size(confStack{1,1},2),size(confStack,1),'uint8');
decSted3D = zeros(size(confStack{1,1},1),size(confStack{1,1},2),size(confStack,1),'uint8');

for stackIdx = 1:size(confStack,1)
    % 3D arrays
    conf3D(:,:,stackIdx)    = confStack{stackIdx,1};
    sted3D(:,:,stackIdx)    = stedStack{stackIdx,1};
    % deconv
    decConf3D(:,:,stackIdx) = decConfStack{stackIdx,1};
    decSted3D(:,:,stackIdx) = decStedStack{stackIdx,1};
end

%% ------------------------------------------------------------------------
% Go through 3 dimensions (original code executed dim = 1)
% -------------------------------------------------------------------------
for dim = 1  % original script processed only the native orientation
    if dim == 1
        % leave stack
        workingConf = conf3D;
        workingSted = sted3D;
    elseif dim == 2
        % rotate stack
        workingConf = imrotate3(conf3D,90,[0 1 0],"nearest","loose");
        workingSted = imrotate3(sted3D,90,[0 1 0],"nearest","loose");
    elseif dim == 3
        % rotate stack
        workingConf = imrotate3(conf3D,90,[1 0 0],"nearest","loose");
        workingSted = imrotate3(sted3D,90,[1 0 0],"nearest","loose");
    end

    stackN = size(workingSted,3);

    for stackIdx = 1:stackN

        %% ----------------------------------------------------------------
        % Select image & prepare output paths
        % -----------------------------------------------------------------
        % check if output exists
        subDirName = strrep(strrep(strrep(fileName(1:end-1),'\',''),':',''),'.','');
        outFolder  = [outDir,filesep,subDirName];
        if ~exist(outFolder,'dir'), mkdir(outFolder); end

        outPathPre = [outFolder,filesep,suffix, sprintf('slice_%d',stackIdx)];
        segmentationOutputFiles = dir([outPathPre,'*segmentation.mat']);

        if isempty(segmentationOutputFiles)  % only analyse if new
            tic;

            confImg    = uint8(workingConf(:,:,stackIdx));   % confocal
            stedImg    = uint8(workingSted(:,:,stackIdx));   % sted
            decConfImg = uint8(decConf3D(:,:,stackIdx));     % deconv confocal
            decStedImg = uint8(decSted3D(:,:,stackIdx));     % deconv sted

            %% ------------------------------------------------------------
            % Denoising (DnCNN) + binarization
            % -------------------------------------------------------------
            net = denoisingNetwork('DnCNN');

            if dim == 1
                denConfImg = denoiseImage(uint8(confImg),net);
                denStedImg = denoiseImage(uint8(stedImg),net);
                % overwrite stacks with denoised
                conf3D(:,:,stackIdx) = denConfImg;
                sted3D(:,:,stackIdx) = denStedImg;

                binConfImg = imbinarize(denConfImg,'adaptive','Sensitivity',.5);
                binStedImg = imbinarize(denStedImg,'adaptive','Sensitivity',.5);

            elseif dim == 2
                % rotate-stack cases: resize then smooth/adjust and binarize
                confImg = imresize(confImg,[size(conf3D,1),size(conf3D,3)*10]);
                stedImg = imresize(stedImg,[size(conf3D,1),size(conf3D,3)*10]);

                denConfImg = imgaussfilt(imadjust((confImg)),10);
                denStedImg = imgaussfilt(imadjust((stedImg)),10);

                binConfImg = imbinarize(denConfImg);
                binStedImg = imbinarize(denStedImg);

            elseif dim == 3
                confImg = imresize(confImg,[size(conf3D,2),size(conf3D,3)*10]);
                stedImg = imresize(stedImg,[size(conf3D,2),size(conf3D,3)*10]);

                denConfImg = imgaussfilt(imadjust((confImg)),10);
                denStedImg = imgaussfilt(imadjust((stedImg)),10);

                binConfImg = imbinarize(denConfImg);
                binStedImg = imbinarize(denStedImg);
            end

            %% ------------------------------------------------------------
            % Remove noise & morphological cleanup
            % ------------------------------------------------------------
            BWConf = bwmorph(binConfImg,'open',1);
            BWSted = bwmorph(binStedImg,'open',1);

            % remove small components
            BWConf = bwareaopen(BWConf,1000);
            BWSted = bwareaopen(BWSted,1000);

            % close
            se = strel('disk',7,8);
            BWConf = imclose(BWConf,se);
            BWSted = imclose(BWSted,se);

            % thicken
            BWConf = bwmorph(BWConf,'thicken',4);
            BWSted = bwmorph(BWSted,'thicken',4);

            % merge both for better result and comparability
            membraneBW = BWConf | BWSted;

            %% ------------------------------------------------------------
            % Compute boundaries (membrane contours) and smooth
            % ------------------------------------------------------------
            memBounds = bwboundaries(membraneBW);
            bdyBW = zeros(size(stedImg,1),size(stedImg,2),'logical'); 
            clear smBdy;

            % filter boundaries by minimum length
            bdyIdc = 1:length(memBounds);
            sel = ones(1,length(memBounds));
            for boundIdx = 1:length(memBounds)
                bdy = memBounds{boundIdx};
                if size(bdy,1) < 15
                    sel(boundIdx) = 0;
                end
            end
            memBounds = {memBounds{bdyIdc(logical(sel))}}';

            % smooth boundaries and compute normals
            for boundIdx = 1:length(memBounds)
                bdy = memBounds{boundIdx};
                x = bdy(:,1);
                y = bdy(:,2);

                % smooth boundary (periodic)
                extra = 3;
                pntIdc = (-extra+1:length(x)+extra)';
                xfit = fit(pntIdc,[x(end+1-extra:end);x;x(1:extra)],'smoothingspline','smoothingparam',0.01);
                yfit = fit(pntIdc,[y(end+1-extra:end);y;y(1:extra)],'smoothingspline','smoothingparam',0.01);
                xsm = xfit(1:length(x));
                ysm = yfit(1:length(x));

                % at each point of the boundary find the normal vector
                dX = diff(xsm);
                dY = diff(ysm);
                % change orientation based on curvature
                ddX = diff(dX);
                ddY = diff(dY);
                ddX = [ddX; ddX(1)];
                ddY = [ddY; ddY(1)];
                curv = mean((dX.*ddY-ddX.*dY)./sqrt(dX.*dX+dY.*dY).^3);
                if curv >=0
                    normals = [-dY,dX];
                else
                    normals = [dY,-dX];
                end
                normals = normals./sqrt(normals(:,1).^2+normals(:,2).^2);

                smBdy(boundIdx).pnts    = [xsm,ysm];
                smBdy(boundIdx).normals = normals;

                % sample in pixels (render smoothed boundary mask)
                for pntIdx = 1:length(xsm)-1
                    % sample on pixel grid
                    bdyBW(round(xsm(pntIdx)),round(ysm(pntIdx))) = 1;
                end
            end

            %% ------------------------------------------------------------
            % Add skeleton and trace central lines
            % ------------------------------------------------------------
            skelBW = bwskel(membraneBW);
            % ctrBdy = traceSkel(skelBW) returns list of centerline polylines
            ctrBdy = traceSkel(skelBW);

            ctrBW = zeros(size(stedImg,1),size(stedImg,2),'logical'); 
            clear smCtrBdy;

            % filter short centerlines
            nPnts = zeros(1,length(ctrBdy));
            for ctrIdx = 1:length(ctrBdy)
                nPnts(ctrIdx) = size(ctrBdy(ctrIdx).pnts,1);
            end
            selIdc = 1:length(ctrBdy);
            selIdc = selIdc(nPnts>20);
            ctrBdy = ctrBdy(selIdc);

            % smooth centerlines and compute normals
            for ctrIdx = 1:length(ctrBdy)
                ctrPnts = ctrBdy(ctrIdx).pnts;
                x = ctrPnts(:,1);
                y = ctrPnts(:,2);

                extra = 0;
                pntIdc = (-extra+1:length(x)+extra)';
                xfit = fit(pntIdc,[x(end-extra+1:end);x;x(1:extra)],'smoothingspline','smoothingparam',0.01);
                yfit = fit(pntIdc,[y(end-extra+1:end);y;y(1:extra)],'smoothingspline','smoothingparam',0.01);

                xsm = xfit(1:length(x));
                ysm = yfit(1:length(x));

                % at each point of the boundary find the normal vector
                dX = diff(xsm);
                dY = diff(ysm);
                % change orientation based on curvature
                ddX = diff(dX);
                ddY = diff(dY);
                ddX = [ddX; ddX(1)];
                ddY = [ddY; ddY(1)];
                curv = mean((dX.*ddY-ddX.*dY)./sqrt(dX.*dX+dY.*dY).^3);
                if curv >=0
                    normals = [-dY,dX];
                else
                    normals = [dY,-dX];
                end
                normals = normals./sqrt(normals(:,1).^2+normals(:,2).^2);

                smCtrBdy(ctrIdx).pnts    = [xsm,ysm];
                smCtrBdy(ctrIdx).normals = normals;

                % render centerline mask (for visualization parity)
                for pntIdx = 1:length(xsm)-1
                    ctrBW(round(xsm(pntIdx)),round(ysm(pntIdx))) = 1;
                end
            end

            %% ------------------------------------------------------------
            % Trace central lines with length parametrization
            % ------------------------------------------------------------
            clear smParamBdy;
            for ctrIdx = 1:length(smCtrBdy)
                ctrPnts = smCtrBdy(ctrIdx).pnts;

                % arc-length parameterization
                x = ctrPnts(:,1);
                y = ctrPnts(:,2);
                dx = diff(x); dy = diff(y);
                ds = sqrt(dx.^2+dy.^2);
                s  = cumsum(ds);
                s  = [0; s];

                extra = 0;
                sfit = [s(end+1-extra:end);s;s(1:extra)];
                xfit = fit(sfit,[x(end+1-extra:end);x;x(1:extra)],'smoothingspline','smoothingparam',0.01);
                yfit = fit(sfit,[y(end+1-extra:end);y;y(1:extra)],'smoothingspline','smoothingparam',0.01);

                dSample  = 4;                      % desired fixed length increment
                NSamples = max(s)/dSample;
                fprintf('%d length samples.\n',NSamples);

                sampleS = linspace(min(s),max(s),NSamples);
                sampleS = [sampleS,sampleS(1)];
                xsm = xfit(sampleS);
                ysm = yfit(sampleS);

                % normals from length-parameterized curve
                dX = diff(xsm);
                dY = diff(ysm);
                ddX = diff(dX);
                ddY = diff(dY);
                ddX = [ddX; ddX(1)];
                ddY = [ddY; ddY(1)];
                curv = mean((dX.*ddY-ddX.*dY)./sqrt(dX.*dX+dY.*dY).^3);
                if curv >=0
                    normals = [-dY,dX];
                else
                    normals = [dY,-dX];
                end
                normals = normals./sqrt(normals(:,1).^2+normals(:,2).^2);

                smParamBdy(ctrIdx).pnts    = [xsm(1:end-1),ysm(1:end-1)];
                smParamBdy(ctrIdx).normals = normals;
            end

            %% ------------------------------------------------------------
            % Visualization: membrane & normals overlay (plot lines)
            % ------------------------------------------------------------
            fac = 30;                 % normal sampling length (pixels)
            close all;
            plot_lines = true;
            if plot_lines
                fig = figure;
                surf(imadjust(stedImg'),'EdgeColor','none')
                view(2); colormap('gray')

                % figure geometry similar to original script
                fig.Units = 'centimeters';
                fig.Position(2) = 3;
                fig.Position(3) = 10;
                fig.Position(4) = fig.Position(3)*size(stedImg,2)/size(stedImg,1);
                ax = gca;
                ax.Position(1:2) = 0; ax.Position(3:4) = 1;

                hold on
                % membrane boundaries in light color
                for memIdx = 1:length(smBdy)
                    pnts = smBdy(memIdx).pnts;
                    plot3(pnts(:,1),pnts(:,2),ones(1,length(pnts(:,2)))*1000,'color',[0.9 0.9 0.8],'LineWidth',1.4)
                end

                % centerlines + normals
                colors = lines(length(smParamBdy));
                for ctrIdx = 1:length(smParamBdy)
                    ctrPnts  = smParamBdy(ctrIdx).pnts;
                    ctrNorms = smParamBdy(ctrIdx).normals;

                    % exclude boundary points for plotting
                    ctrPnts  = ctrPnts(2:end-1,:);
                    ctrNorms = ctrNorms(2:end-1,:);

                    for pntIdx = 1:size(ctrPnts,1)
                        plot3(ctrPnts(pntIdx,1),ctrPnts(pntIdx,2),1000,'color',colors(ctrIdx,:), ...
                              'LineStyle','none','Marker','.');
                        plot3(ctrPnts(pntIdx,1)+[0 ctrNorms(pntIdx,1)]*fac, ...
                              ctrPnts(pntIdx,2)+[0 ctrNorms(pntIdx,2)]*fac, ...
                              1000*ones(1,2),'-r')
                        plot3(ctrPnts(pntIdx,1)-[0 ctrNorms(pntIdx,1)]*fac, ...
                              ctrPnts(pntIdx,2)-[0 ctrNorms(pntIdx,2)]*fac, ...
                              1000*ones(1,2),'-g')
                    end
                end
                hold off; axis equal
                toc
                xlim([1,size(stedImg,1)]); ylim([1,size(stedImg,2)]);
                print(fig,[outPathPre '.png'],'-dpng');
            end

            %% ------------------------------------------------------------
            % Sample image profiles across normals (confocal/STED & deconv)
            % ------------------------------------------------------------
            clear profiles
            for ctrIdx = 1:length(smParamBdy)
                ctrPnts  = smParamBdy(ctrIdx).pnts;
                ctrNorms = smParamBdy(ctrIdx).normals;

                % exclude boundary
                ctrPnts  = ctrPnts(2:end-1,:);
                ctrNorms = ctrNorms(2:end-1,:);

                % sample image profiles
                stedSamples    = [];
                confSamples    = [];
                decStedSamples = [];
                decConfSamples = [];

                for pntIdx = 1:size(ctrPnts,1)
                    % build line along normal; unify direction via start->end
                    linePnts = [ ...
                        ctrPnts(pntIdx,2)+linspace(-ctrNorms(pntIdx,2)*fac, ctrNorms(pntIdx,2)*fac, 1+2*ceil(fac)); ...
                        ctrPnts(pntIdx,1)+linspace(-ctrNorms(pntIdx,1)*fac, ctrNorms(pntIdx,1)*fac, 1+2*ceil(fac)) ...
                    ]';
                    linePnts = round(linePnts);  % pixel grid

                    % sample from images (try/catch per-point as in original)
                    samples    = zeros(2,size(linePnts,1));
                    decSamples = zeros(2,size(linePnts,1));
                    for linePntIdx = 1:size(linePnts,1)
                        try
                            samples(1,linePntIdx)    = stedImg(   linePnts(linePntIdx,2), linePnts(linePntIdx,1));
                            samples(2,linePntIdx)    = confImg(   linePnts(linePntIdx,2), linePnts(linePntIdx,1));
                            decSamples(1,linePntIdx) = decStedImg(linePnts(linePntIdx,2), linePnts(linePntIdx,1));
                            decSamples(2,linePntIdx) = decConfImg(linePnts(linePntIdx,2), linePnts(linePntIdx,1));
                        end
                    end

                    stedSamples    = [stedSamples,    samples(1,:)']; 
                    confSamples    = [confSamples,    samples(2,:)']; 
                    decStedSamples = [decStedSamples, decSamples(1,:)']; 
                    decConfSamples = [decConfSamples, decSamples(2,:)']; 
                end

                profiles(ctrIdx,1).sted    = stedSamples;
                profiles(ctrIdx,1).conf    = confSamples;
                profiles(ctrIdx,1).decSted = decStedSamples;
                profiles(ctrIdx,1).decConf = decConfSamples;
            end

            %% ------------------------------------------------------------
            % Plot unwrapped images (stacked confocal/STED along profiles)
            % ------------------------------------------------------------
            plot_unwrapped_images = true;
            if plot_unwrapped_images
                stedUnwrapped = [];
                confUnwrapped = [];
                colors = lines(length(profiles));

                % overlay centerlines on STED image with labels
                fig = figure;
                surf(imadjust(stedImg'),'EdgeColor','none')
                view(2); colormap('gray')
                xlim([1,size(stedImg,1)]); ylim([1,size(stedImg,2)])

                fig.Units = 'centimeters';
                fig.Position(2) = 3;
                fig.Position(3) = 10;
                fig.Position(4) = fig.Position(3)*size(stedImg,2)/size(stedImg,1);
                ax = gca; ax.Position(1:2) = 0; ax.Position(3:4) = 1;

                hold on
                for ctrIdx = 1:length(profiles)
                    ctrPnts = smParamBdy(ctrIdx).pnts;
                    plot3(ctrPnts(:,1),ctrPnts(:,2),ones(1,size(ctrPnts,1))*1000, ...
                          'color',colors(ctrIdx,:),'LineWidth',1.6)

                    % use deconvolved profiles for unwrapping (as in original)
                    use_deconv = true;
                    if use_deconv
                        stedUnwrapped = [stedUnwrapped,(profiles(ctrIdx).decSted)]; 
                        confUnwrapped = [confUnwrapped,(profiles(ctrIdx).decConf)]; 
                    else
                        stedUnwrapped = [stedUnwrapped,(profiles(ctrIdx).sted)];    
                        confUnwrapped = [confUnwrapped,(profiles(ctrIdx).conf)];    
                    end

                    % annotate centerline id next to its mid-point
                    pos = ax.Position;
                    xPosition = ctrPnts(round(length(ctrPnts)/2),1);
                    yPosition = ctrPnts(round(length(ctrPnts)/2),2);
                    posx = (xPosition-min(xlim))/(max(xlim)-min(xlim))*pos(3)+pos(1);
                    posy = (yPosition-min(ylim))/(max(ylim)-min(ylim))*pos(4)+pos(2);
                    annotation('textbox','Position',[posx,posy,0,0], ...
                               'String',sprintf('%02d',ctrIdx), ...
                               'EdgeColor','y','color',colors(ctrIdx,:), ...
                               'FontSize',16);
                end
                hold off
                print(gcf,[outPathPre, 'mem.png'],'-dpng');

                % build "membrane unwrap" as [conf; sted] and mark segment IDs
                memUnWrap = [imadjust(uint8(confUnwrapped)); imadjust(uint8(stedUnwrapped))];

                fig = figure;
                surf(flip(memUnWrap,1),'EdgeColor','none')
                view(2); colormap('gray')

                % draw colored separators for each traced segment
                clear lengths;
                for ctrIdx = 1:length(profiles)
                    lengths(ctrIdx) = size(profiles(ctrIdx).sted,2); 
                end
                linePos = [0,cumsum(lengths)];

                hold on
                for ctrIdx = 1:length(profiles)
                    plot3([linePos(ctrIdx:ctrIdx+1)],[1,1],[1000,1000], ...
                          'color',colors(ctrIdx,:),'LineWidth',10);
                end
                hold off

                xlim([1,size(memUnWrap,2)]); ylim([1,size(memUnWrap,1)]);
                fig.Units = 'centimeters';
                fig.Position(2) = 3;
                fig.Position(3) = 10;
                fig.Position(4) = fig.Position(3)*size(memUnWrap,1)/size(memUnWrap,2);
                ax = gca; ax.Position(1:2) = 0; ax.Position(3:4) = 1;

                % label segment indices along top
                pos = ax.Position; 
                dlinePos = diff(linePos); 
                for ctrIdx = 1:length(profiles)
                    xPosition = linePos(ctrIdx);
                    yPosition = 30;
                    posx = (xPosition-min(xlim))/(max(xlim)-min(xlim))*ax.Position(3)+ax.Position(1);
                    posy = (yPosition-min(ylim))/(max(ylim)-min(ylim))*ax.Position(4)+ax.Position(2);
                    annotation('textbox','Position',[posx,posy,0,0], ...
                               'String',sprintf('%02d',ctrIdx), ...
                               'EdgeColor','y','color',colors(ctrIdx,:), ...
                               'FontSize',16);
                end

                print(gcf,[outPathPre,'mem_unwrap.png'],'-dpng');
            end

            %% ------------------------------------------------------------
            % Save everything
            % ------------------------------------------------------------
            clear segmentation;
            segmentation.profiles    = profiles;
            segmentation.confImg     = confImg;
            segmentation.stedImg     = stedImg;
            segmentation.membrane    = smParamBdy;
            segmentation.membraneBdy = smBdy;
            save([outPathPre 'segmentation.mat'],'segmentation');

        end % if ~exist(segmentation)
    end     % for stackIdx
end         % for dim
end         % function
