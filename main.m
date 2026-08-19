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
%    OAT.mat                     contain OAT run data 
%    Global.mat                  Contain morris runing data for 20 trajectory 
%  This file declares run controls only. All derived quantities, energy
%  accounts and validity screens are computed and reported inside the model.
%==========================================================================

clc; clear; close all

% FC = Baseline_Input_Variable();  % Baseline variable 

%% ---- Baseline run -------------------------------------------------------
%  Nx = 40, the grid verified by the GCI study.

%  The model draws the Results-section figures itself when FC.quiet is
%  false, so plot_baseline is not called again here.

% results = ZLDD_complete_modeling(FC);

%% ---- Parametric sweep ---------------------------------------------------
%  Sensitivity is a FUNCTION; call it from the Command Window. The sweep
%  overrides the grid to Nx = 20 for run economy, so its results are read
%  RELATIVE to run 1 and absolute values are not comparable with the
%  baseline figures above.
%
%  Sensitivity('dryrun', true);        % list runs, solve nothing
%      T = Sensitivity('only', {'Tfeed'});  % single factor
%      T = Sensitivity();                  % full sweep

%        load('OAT.mat');
% 
%      Sensitivity('figs',  T);                    % figure suite
%       Sensitivity('tornado', T, 'R_still_pct');   % tornado on any response
%       Sensitivity('pinch',  T);                   % SEC against cascade recovery
%       Sensitivity('validity', T);                 % validity envelope, on request only
%      E = Sensitivity('elast', T);                % elasticity ranking
% 
% 
%     Sensitivity('export', 'svg');                 % every open figure -> vector






%% ---- Morris screening ---------------------------------------------------
%  What OAT cannot give: sigma, the dependence of a factor's effect on where
%  the other factors sit. Cost is (k+1)*r solves.
%
%      Morris_Global('dryrun', true);                % factor set and cost, seconds



% M = Morris_Global('r', 20, 'timeout', 900);  % r=20 , take almost 12 hours

%    load('Global.mat');                      % load morris data 
% h = Morris_Global('figs', M);                 % display plot 
% Sensitivity('export','svg',h);                % Save plot 

%==========================================================================
