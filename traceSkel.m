function lines = traceSkel(skelBW, doPlot)
% TRACE SKEL
% Description :
%   Trace skeletonized membrane centerlines from a binary skeleton image.
%   The procedure:
%     1) Detect bifurcation centers (≥3 neighbors) and remove them
%        to split skeleton into simple branches.
%     2) Detect endpoints (1 neighbor).
%     3) For each connected component (branch), pick an endpoint (or
%        arbitrary pixel if none exist) and iteratively walk to the
%        nearest remaining skeleton pixel until distance exceeds ~2 px.
%     4) Return each traced polyline as lines(k).pnts with coordinates
%        stored as [x, y].
%     5) Optionally plot traced lines over the skeleton if doPlot=true.
%
% Usage :
%   % Given a binary skeleton image skelBW:
%   lines = traceSkel(skelBW);           % no plotting
%   lines = traceSkel(skelBW, true);     % with plotting
%
% Dependencies :
%   - MATLAB Image Processing Toolbox (makelut, bwlookup, regionprops)
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
% Default arguments
% -------------------------------------------------------------------------
if nargin < 2
    doPlot = false;
end

%% ------------------------------------------------------------------------
% Detect bifurcation centers (neighbor count >= 3) and remove them
% -------------------------------------------------------------------------
numberNeighboringPixels = 3;
lut   = makelut(@(x) sum(x(:)) >= (numberNeighboringPixels+1), 3);
lutBW = bwlookup(skelBW, lut) & skelBW;
stats = regionprops(lutBW, 'Centroid');
bifPos = [stats.Centroid];

% take out bifurcations
skelBWTrace = skelBW;
bifX = bifPos(1:2:end);
bifY = bifPos(2:2:end);
for bifIdx = 1:length(bifX)
    for idx1 = -1:1
        for idx2 = -1:1
            skelBWTrace(round(bifY(bifIdx))+idx1, round(bifX(bifIdx))+idx2) = 0;
        end
    end
end

%% ------------------------------------------------------------------------
% Detect endpoints (neighbor count == 1)
% -------------------------------------------------------------------------
numberNeighboringPixels = 1;
lut   = makelut(@(x) sum(x(:)) == (numberNeighboringPixels+1), 3);
lutBW = bwlookup(skelBWTrace, lut) & skelBWTrace;
stats = regionprops(lutBW, 'Centroid');
endPos = [stats.Centroid];
endX   = endPos(1:2:end);
endY   = endPos(2:2:end);
endPnt = [endX', endY'];

%% ------------------------------------------------------------------------
% Fallback when no endpoints exist: choose arbitrary skeleton pixel
% -------------------------------------------------------------------------
[y, x] = find(skelBWTrace);
skelPnts = [x, y];
if isempty(endPnt)
    endPnt = skelPnts(1,:);
end

%% ------------------------------------------------------------------------
% Trace each connected component as a polyline via greedy nearest neighbor
% -------------------------------------------------------------------------
stats = regionprops(skelBWTrace, 'PixelList');
for statsIdx = 1:length(stats)
    pnts = stats(statsIdx).PixelList;  % [x y] integer pixel coordinates

    % find end point that is in this component
    endPntFound = false;
    for endPntIdx = 1:size(endPnt,1)
        thisEndPnt = endPnt(endPntIdx,:);
        relVec = round(thisEndPnt) - round(pnts);
        dists  = sqrt(relVec(:,1).^2 + relVec(:,2).^2);
        minIdx = find(dists==min(dists), 1);
        if dists(minIdx) < 1.6
            endPntFound = true;
            break
        end
    end
    % if there was no endpoint in this component choose arbitrary
    if ~endPntFound
        thisEndPnt = pnts(1,:);
    end

    % iterative walk along nearest neighbors
    linePoints = [];
    cPnt = thisEndPnt;
    skelDistPnts = pnts;

    % remove first point
    distVec = cPnt - skelDistPnts;
    dists   = sqrt(distVec(:,1).^2 + distVec(:,2).^2);
    pntIdx  = find(dists==min(dists), 1, 'first');
    cond = true(1, size(skelDistPnts,1));
    cond(pntIdx) = false;
    skelDistPnts = skelDistPnts(cond,:);

    % save first point
    linePoints = [linePoints; cPnt];
    minDist = 1;

    while minDist < 2
        distVec = cPnt - skelDistPnts;
        dists   = sqrt(distVec(:,1).^2 + distVec(:,2).^2);
        minDist = min(dists);
        pntIdx  = find(dists==min(dists), 1, 'first');

        % determine new point
        cPnt = skelDistPnts(pntIdx,:);

        % remove visited point
        cond = true(1, size(skelDistPnts,1));
        cond(pntIdx) = false;
        skelDistPnts = skelDistPnts(cond,:);

        % append
        linePoints = [linePoints; cPnt];
    end

    % store as [x y] (flip from PixelList order)
    lines(statsIdx).pnts = flip(linePoints, 2);
end

%% ------------------------------------------------------------------------
% Optional visualization of traced lines
% -------------------------------------------------------------------------
if doPlot
    hold on;
    for lineIdx = 1:length(lines)
        pnts = lines(lineIdx).pnts; % [x y]
        plot(pnts(:,1), pnts(:,2), ...
            'color', [lineIdx/length(lines) 0 1-lineIdx/length(lines)], ...
            'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 20);
    end
    hold off;
end
end
