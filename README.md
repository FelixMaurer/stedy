# stedy
Image processing of 3D confocal and STED stacks. Developed for the tracing and geometrical analysis of biological membranes in edge-on view.
STED Membrane Profiling (MATLAB)
================================

Segment membranes from 3D confocal/STED stacks, trace centerlines, and
extract intensity profiles along geometric normals — using Bio-Formats,
DnCNN denoising, morphological cleanup, skeletonization, spline smoothing,
and arc-length reparameterization.

Overview
--------
This project refactors an analysis workflow into:
- process_sted_file.m  : Function to process a single .msr file (and
  matching deconvolved .obf) and write figures + segmentation.mat
- run_sted_batch.m     : Driver script that recursively discovers .msr
  files, filters by name, summarizes, and calls the function
- traceSkel.m          : Utility to trace skeletonized centerlines
  (with optional plotting flag)

Pipeline (high level)
---------------------
1. Load confocal & STED channels from .msr with Bio-Formats, plus the
   matching deconvolved .obf.
2. Denoise (DnCNN) and binarize each slice.
3. Morphological cleanup (open → area filter → close → thicken);
   fuse masks:
      membraneBW = BWConf | BWSted
4. Contours: bwboundaries → smoothing splines.
5. Normals from tangents, sign chosen by mean curvature
6. Skeletonize membranes (bwskel), trace centerlines (traceSkel),
   re-smooth, then arc-length reparameterize,
   resample at Δs=4 px.
7. Sample profiles along normals
   for both raw and deconvolved channels.
8. Save overlay figures and segmentation.mat.

Repository Structure
--------------------
process_sted_file.m
run_sted_batch.m
traceSkel.m
example/
  data/
    sample_sted_pcExc.msr
    sample_sted_pcExc_deconv.obf
    sample_sted_pcExc_conf.gif
    sample_sted_pcExc_sted.gif
  demo_run.m
README.txt

Example Data
--------------------
- scan of fixed acanthocyte of patient stained with Atto647N
- deconvoluted sted image
- animated .gif visualization of confocal and sted stack

Requirements
------------
- MATLAB R2021b or newer
- Toolboxes: Image Processing, Curve Fitting, Deep Learning
- Bio-Formats for MATLAB (bfopen)
- Tested on Windows/macOS/Linux supported by MATLAB

Installation
------------
1. Clone or download the repository.
2. Install Bio-Formats (bfmatlab) and add it to MATLAB path.
3. Open MATLAB in the repo root folder.

Quick Start
-----------
A. Process a single file:
   msr = 'example\data\sample_sted_pcExc.msr';
   out = 'example\processed';
   process_sted_file(msr, out);

B. Batch process a folder tree:
   Edit run_sted_batch.m and set:
      outDir   = 'processed';
      stedRoot = 'data';
   Run run_sted_batch.

Example Data & Testing
----------------------
The example/ folder contains:
- example\data\sample_sted_pcExc.msr
- example\data\sample_sted_pcExc_deconv.obf

Run demo_run.m inside example/:
   cd example
   run('demo_run.m')

Expected outputs per slice:
- out.png          : STED overlay with contours + centerlines
- mem.png          : Centerline overlay and labels
- mem_unwrap.png   : Unwrapped profile image ([conf; STED])
- segmentation.mat : Struct with profiles, images, and geometry

Usage Details
-------------
process_sted_file(msrFilePath, outDir)
   Inputs: .msr file path, output directory
   Behavior: loads data, finds deconv .obf, processes native orientation,
             skips slices already processed.

run_sted_batch.m
   Discovers .msr files under stedRoot, filters by name, calls
   process_sted_file() for each.

traceSkel(skelBW, doPlot)
   skelBW : logical skeleton image (from bwskel)
   doPlot : optional boolean (default false)
   Returns: array of structs with fields lines(k).pnts = [x y].

Troubleshooting
---------------
- Undefined function 'bfopen':
  Bio-Formats not on MATLAB path.
- Cannot find deconvolved .obf:
  Ensure naming matches pattern and file is nearby.
- fit errors:
  Requires Curve Fitting Toolbox.
- Out of memory:
  Try fewer slices or disable figures.

Citation
--------
- Zhang & Zuo et al., "Beyond a Gaussian Denoiser: Residual Learning of
  Deep CNN for Image Denoising (DnCNN)", IEEE TIP, 2017.
- Bio-Formats: https://www.openmicroscopy.org/bio-formats/

License
-------
MIT, © 2025 Felix Maurer
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
