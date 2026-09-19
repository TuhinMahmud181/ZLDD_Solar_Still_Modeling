%==========================================================================
%  ZLDD CASCADE -- DRIVER
%
%  Author : Tuhin Mahmud
%
%  Files:
%    Baseline_Input_Variable.m   design case (single declaration)
%    ZLDD_complete_modeling.m    model
%    plot_baseline.m             Results-section figures
%    Sensitivity.m               local (OAT) sensitivity, factor declaration
%    Morris_Global.m             global (Morris) screening
%    pvgis_irradiance_data.m     plane-of-array irradiance series
%    OAT.mat                     stored OAT sweep results
%    Global.mat                  stored Morris results, r = 20 trajectories
%
%  This file declares run controls only. All derived quantities, energy
%  accounts and validity screens are computed and reported inside the model.
%==========================================================================

clc; clear; close all

%% ---- Baseline run -------------------------------------------------------
%  Nx = 40, the grid verified by the GCI study. The model draws the
%  Results-section figures itself when FC.quiet is false, so plot_baseline
%  is not called again here.

FC      = Baseline_Input_Variable();
%   results = ZLDD_complete_modeling(FC);
%  h = plot_baseline(results);                    % display only
%  Sensitivity('export', 'svg', h);               % every open figure -> vector
%% ---- Parametric sweep ---------------------------------------------------
%  Sensitivity is a FUNCTION; call it from the Command Window. The sweep
%  overrides the grid to Nx = 20 for run economy, so its results are read
%  RELATIVE to run 1 and absolute values are not comparable with the
%  baseline figures above.
%
%      Sensitivity('dryrun', true);                % list runs, solve nothing
%      T = Sensitivity('only', {'Tfeed'});         % single factor
%      T = Sensitivity();                          % full sweep -> OAT.mat
%      load('OAT.mat');                            % or reload a finished sweep
%
%      Sensitivity('figs', T);                     % figure suite
%      Sensitivity('tornado', T, 'R_still_pct');   % tornado on any response
%      Sensitivity('pinch', T);                    % SEC against cascade recovery
%      Sensitivity('validity', T);                 % validity envelope, on request
%      E = Sensitivity('elast', T);                % elasticity ranking
%      Sensitivity('export', 'svg');               % every open figure -> vector

%% ---- Morris screening ---------------------------------------------------
%  What OAT cannot give: sigma, the dependence of a factor's effect on where
%  the other factors sit. Cost is (k+1)*r solves.
%
%      Morris_Global('dryrun', true);              % factor set and cost, seconds
%      M = Morris_Global('r', 20, 'timeout', 900); % r = 20, roughly 12 hours
%      load('Global.mat');                         % or reload a finished screening
%
%      h = Morris_Global('figs', M);               % display plots
%      Sensitivity('export', 'svg', h);            % save plots
%
%  Global.mat carries 240 runs, of which 16 returned no solution and a
%  further 7 failed a validity screen. Those entries stay in place as NaN so
%  that the affected elementary effects are reported as unavailable. Do not
%  compact the arrays: rows, runs, done, Y, V, Xapp, OK and WALL are aligned
%  by row index, and deleting entries from some but not all of them silently
%  mispairs every trajectory downstream.
%==========================================================================
