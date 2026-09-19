function FC = Baseline_Input_Variable()
%BASELINE_INPUT_VARIABLE  Canonical baseline configuration for the ZLDD cascade.
%
%   FC = Baseline_Input_Variable() returns the configuration struct consumed
%   by ZLDD_complete_modeling(). This file is the SINGLE declaration of the
%   design case: main.m, the parametric sweep and the TEA all read it, and no
%   other copy of these numbers should exist anywhere in the project.
%
%   FLOWSHEET: CLOSED air loop, zero make-up. The dry-air flow is conserved,
%   no ambient air enters, and seawater is the only route by which water
%   enters the plant.
%
%   USAGE
%       FC  = Baseline_Input_Variable();
%       res = ZLDD_complete_modeling(FC);
%
%   To vary a parameter, COPY AND EDIT THE RETURNED STRUCT. Never edit this
%   file for a one-off run -- an edit here silently changes every result the
%   project has produced, including figures already in the manuscript:
%       FC = Baseline_Input_Variable();
%       FC.Qfan = 0.18;                     % a study, not a new baseline
%
%   GRID. Nx = 40 here, the grid the GCI study verified. zldd_OAT overrides it
%   to Nx = 20 for run economy and records the value used on every row, so the
%   published baseline and the sweep run on different grids by design, from
%   this one declaration. Sweep results are read RELATIVE to run 1.
%
%   WHY THIS RETURNS A VALUE. The original declaration had no output argument,
%   so FC was built in the function's own workspace and destroyed on return.
%   Do not reach for GLOBAL or ASSIGNIN instead: the baseline has to be a
%   returnable value so it can be copied, perturbed and RECORDED. The sweep
%   stores the whole struct as provenance (prov.baseline_FC).
%
%   TWO ITEMS BELOW ARE MARKED "DECISION": the grid and the fan flow. Each is
%   a place where a value here disagreed with another part of the project.
%   They are flagged rather than silently resolved, because picking one
%   arbitrarily is how two baselines drift apart without an error being raised.
%
%   Units are SI throughout except where a field is explicitly labelled
%   otherwise (Vfeed in L/day, angles in degrees, t_operating in hours).
%
%   See also ZLDD_COMPLETE_MODELING, ZLDD_OAT, PVGIS_IRRADIANCE_DATA.

FC = struct();

% ---- Geometry ------------------------------------------------------------
FC.L        = 1.0;        % [m]   plate length (streamwise)
FC.W        = FC.L;        % [m]   plate width
FC.theta    = 2;          % [deg] film-flow tilt. Sets the gap through
                          %       h_bar = hcomp + tan(theta)*L_ch, so theta --
                          %       not hcomp -- is the larger gap-height lever:
                          %       at this tilt the taper contributes ~0.038 m
                          %       against the 0.020 m offset.
FC.beta     = 21;         % [deg] glazing tilt. MUST equal the PVGIS slope in
                          %       the irradiance query at the foot of this
                          %       file. Changing one alone decouples the
                          %       glazing from the incident series.
FC.hcomp    = 0.02;       % [m]   vapour-gap compartment height offset

% DECISION 1 -- GRID. CONFIRMED: Nx = 40 for the published baseline, Nx = 20 for
% the sweep. Was annotated "(grid-converged)" at 20, which the verification does
% not support -- the GCI study was run at Nx = 40 and reported GCI_80 = 0.24 %,
% certifying 40, not 20. Advection uses first-order upwinding, so truncation error
% scales as dx ~ 1/Nx and halving the grid roughly doubles it. zldd_OAT overrides
% this to 20 for run economy and records the value used on every row.
FC.Nx       = 40;         % [-]   nodes per plate. VERIFIED GRID: GCI_80 = 0.24 %.
FC.Np       = 10;         % [-]   number of absorber-plate stages
FC.fandia   = 0.25;       % [m]   humid-air outlet pipeline diameter

% ---- Optical / material --------------------------------------------------
FC.kappa_w  = 300;        % [1/m] solar-weighted grey absorption coefficient of
                          %       the film. KNOWN BIAS, conceded not hidden: a
                          %       lumped grey fit to a strongly spectral medium,
                          %       over-predicting absorption in the deep stages
                          %       relative to a spectrally resolved model. Swept
                          %       /2 to x2 in the uncertainty block for this reason.

% ---- Ambient (Saudi Arabia, Rayed dataset, March) ------------------------
% Same site as the PVGIS query at the foot of this file. Site, ambient block
% and irradiance query are one set and must be changed together.
FC.Ta       = 308;        % [K]   mean over the operating window, NOT the daily max
FC.Vwind    = 2.0;        % [m/s] mean wind speed; sets the external convective loss
FC.Tground  = 305;        % [K]   site soil temperature

% ---- Feed / operating conditions -----------------------------------------
FC.Tfeed    = 300;        % [K] seawater AS SUPPLIED, at the HX cold side. With
                          %     tfeed_dynamic = true this is a BOUNDARY
                          %     CONDITION; the cascade inlet temperature is
                          %     solved, not specified. Seeding the initial state
                          %     cold is deliberate -- at t = 0 no brine has yet
                          %     been raised.
FC.Vfeed    = 100;        % [L/day] feed volumetric flow rate
FC.TDSfeed  = 35;         % [kg/m3] seawater salinity. PER-VOLUME (Sm) basis, NOT
                          %     g/kg solution (Sg) -- the two are not
                          %     interchangeable and the property correlations are
                          %     fitted on the mass basis.

% DECISION 2 -- FAN FLOW. This file had Qfan = 0.30 while zldd_OAT declares the
% Qfan base as 0.25, and the sweep's factor note is written around 0.25 (gap-1
% humidity w = 0.0238 against w_sat ~ 0.0295, gap Re ~ 2.8e4). Both values appear
% in the level list, so neither is obviously a typo -- but they cannot both be the
% baseline, or every normalised curve and elasticity is anchored to a run 1 that
% does not sit at the declared base. Set to 0.25 to agree with the sweep. If 0.30
% is correct, change it here AND in the 'Qfan' entry of local_factors(); the
% base-consistency assertion in zldd_OAT refuses to run until the two agree.
FC.Qfan     = 0.25;       % [m3/s] fan volumetric flow rate

% ---- Air path: CLOSED LOOP (still -> AWG -> recuperator -> heaters -> still)
% Airtight, zero make-up: dry-air flow is conserved and no ambient air enters.
% Confirmed by the baseline run -- ambient moisture captured S9 = 0, moisture
% rejected to atmosphere S11 = 0, and the thermal-input datum is T_coil rather
% than Ta, because no ambient stream crosses the boundary.
FC.RH_amb   = 0.10;       % [-] ambient relative humidity. RETAINED BUT NOT
                          %     CONSUMED BY THE AIR PATH: the still inlet
                          %     humidity ratio is w_in = w_sat(T_coil), set by
                          %     the coil alone. RH_amb survives only in the
                          %     radiative/sky terms and in diagnostics.

% ---- Coil temperature: FIXED SETPOINT, and it must be --------------------
% The coil sets the still inlet humidity ratio, w_in = w_sat(T_coil), so it feeds
% the cascade directly. Deriving T_coil from the SOLVED exhaust dew point would be
% circular and would need an outer fixed-point iteration; the model rejects
% coil_mode = 'auto' rather than quietly reinterpreting it.
%
% T_coil is a first-order design variable and acts twice: it sets w_in, and it
% sets the heat-pump lift. Lowering it dries the loop air and raises the driving
% potential (Psat(Tw) - Pv) throughout the cascade, at the cost of AWG duty.
% Expect it to move SEC and recovery in OPPOSITE directions.
FC.coil_mode    = 'fixed';
FC.T_coil       = 295;    % [K] AWG coil surface temperature
FC.T_coil_floor = 275;    % [K] frost guard. The model has NO frost logic, so this
                          %     is a hard floor, not a correlation limit.

% ---- Electric air heater: SPECIFIED BY OUTLET TEMPERATURE ---------------
% An independent electrically driven unit. The evaporator vapour heats the
% SEAWATER FEED; the air is heated separately. Because the unit is specified by
% outlet temperature the DUTY is an output, not an input -- which is why this
% setpoint drives the single largest purchased load in the plant.
FC.T_air_in       = 325;  % [K] still inlet air temperature setpoint
FC.eta_air_heater = 0.98; % [-] resistance element efficiency

% ---- Feed preheat HX outlet: SOLVED INSIDE THE MODEL --------------------
% Recomputed at every RHS call from the instantaneous cascade brine:
%     brine -> evaporator vapour -> condensing duty -> feed temperature.
% Closed WITHIN the integration, so outlet temperature and condensed fraction are
% both OUTPUTS. The loop carries negative feedback and is therefore stable:
%   hotter feed -> more evaporation -> less brine -> less vapour -> cooler feed
FC.tfeed_dynamic   = true;
FC.dT_preheat_cap  = 60;  % [K] STARTUP-ONLY cap on the computed feed rise. A
                          %     SOLVER AID, NOT A PHYSICAL LIMIT, released after
                          %     t_preheat_guard so the steady answer carries no
                          %     tuning constant. If the report shows it still
                          %     clamping at t_end, the reported feed temperature
                          %     is (Tfeed + cap) -- a constant you chose, not a
                          %     value the model solved. The ok_preheatcap screen
                          %     in zldd_OAT tests exactly this.
FC.t_preheat_guard = 1800;% [s] window over which the startup cap applies
FC.dT_pinch_HX     = 5;   % [K] minimum HX approach. Sets a THERMODYNAMIC ceiling
                          %     on the feed at T_sat_vap - dT_pinch_HX = 368.15 K;
                          %     inactive while Tfeed_target sits below it.
FC.Tfeed_target    = 320; % [K] DESIGN target at the cascade inlet, set well below
                          %     the pinch ceiling so the film does not enter near
                          %     flash and psat_saline() is not worked at the edge
                          %     of its fit. Vapour the feed does not absorb is
                          %     SURPLUS and is routed to the air preheater.
                          %     MUST STAY IN STEP with the Tfeed_target base in
                          %     zldd_OAT/local_factors().
FC.dT_pinch_cond   = 5;   % [K] minimum approach in the vapour-fired air preheater

% ---- Air-to-air recuperator (still exhaust <-> AWG outlet) --------------
% PASSIVE, and it removes an irreversibility rather than working around one:
% without it the plant cools the loop air across a temperature span and reheats it
% across the same span, paying a compressor for one and a resistance element for
% the other. Both duties shrink together -- the coil sees pre-cooled air, the
% heater pre-warmed air.
%
% eps is defined on the COLD side, which carries Cmin because the hot pass
% condenses. The recoverable ceiling is (Tv_exit - T_coil), so a narrow spread
% caps recovery whatever the effectiveness. Raising eps TIGHTENS the AWG dew-point
% margin, so watch the ok_dewpoint screen rather than assuming more is better.
% Area is large because both sides are gas; trading that area against PV capacity
% for the avoided load is a TEA question (Paper 2), not a process one.
%
% eps = 0 disables the unit entirely and quantifies its value directly.
FC.eps_recup = 0.75;      % [-] counterflow effectiveness, COLD side

% ---- Brine evaporator / crystalliser: CONTINUOUS, at NaCl saturation ----
% THREE DISTINCT TEMPERATURES, easily conflated:
%   T_sat_vap    = CONDENSING temperature in the feed HX
%   T_boil       = T_sat_vap + BPE, the boiling liquor
%   T_vap_heater = T_boil + dT_drive, the element/jacket
%
% The vapour leaving is PURE STEAM: it desuperheats within the first centimetres
% of the exchanger, then condenses ISOTHERMALLY at T_sat_vap regardless of brine
% concentration. BPE therefore never reaches the feed HX; it is a pure heater-side
% penalty. Confirmed numerically -- a +/-1.5 K perturbation moves plant SEC by
% 0.16 % and leaves every cascade metric bit-identical.
%
% BPE is held CONSTANT despite a variable inlet concentration because continuous
% operation pins the liquor at NaCl saturation: past saturation, removing water
% precipitates solid NaCl rather than raising the dissolved concentration, so the
% liquid phase and its water activity are fixed. w_salt_target is the SOLIDS
% FRACTION OF THE DISCHARGED SLURRY, not the concentration of the boiling liquor
% -- these are routinely confused.
FC.bl.w_salt_target = 0.99;   % [-]  salt mass fraction of the discharged slurry
FC.bl.P_evap        = 101325; % [Pa] evaporator operating pressure
FC.bl.T_sat_vap     = 373.15; % [K]  Tsat(P_evap); condensing T in the feed HX
FC.bl.BPE           = 8.5;    % [K]  boiling-point elevation at NaCl saturation
                              %      (Tb ~ 108.5 C at 1 atm). CITE A SOURCE in
                              %      Section 2 -- currently asserted without one.
FC.bl.dT_drive      = 13.0;   % [K]  element driving DT above the boiling liquor
FC.bl.h_cryst       = 65e3;   % [J/kg salt] NaCl crystallisation enthalpy, a duty
                              %      CREDIT. Thermochemical constant, not an
                              %      estimate, hence excluded from the
                              %      uncertainty block.
FC.bl.eta_heater    = 0.98;   % [-]  resistance heater efficiency
FC.bl.f_heatloss    = 0.15;   % [-]  vessel loss as a fraction of useful duty. A
                              %      PURE ENGINEERING ESTIMATE WITH NO STATED
                              %      BASIS, multiplying the second-largest
                              %      purchased load. Swept 0.05-0.25; a parameter
                              %      you invented needs a stated band.

% ---- Operating window ----------------------------------------------------
FC.t_day_start = 8*3600;      % [s]  start of the operating window
FC.t_day_end   = 18*3600;     % [s]  end of the operating window
FC.t           = FC.t_day_end - FC.t_day_start;   % [s]  integration horizon
FC.t_operating = FC.t/3600;   % [hr] DERIVED, so it cannot disagree with the window

% ---- Solver --------------------------------------------------------------
FC.stepSize  = 5000;      % [-] number of OUTPUT time points. Sampling density
                          %     only -- the integrator chooses its own internal
                          %     steps, so this does not set accuracy.
FC.plotstage = 1;         % [-] stage index for the per-stage diagnostic plots

% ---- Reporting depth -----------------------------------------------------
% false -> headline performance, stream table, unit specifications.
% true  -> also the solver audit: per-plate tables, mass and energy closure,
%          interface reciprocity, Reynolds regimes, quasi-steady checks.
% TRUE here because a baseline run should always be auditable. Batch callers such
% as zldd_OAT override it to false; 80 verbose reports are unreadable, and the
% audit is not what a sweep is looking at.
FC.verbose = false;

% ---- PVGIS plane-of-array irradiance. SLOPE MUST MATCH FC.beta ----------
% Same site as the ambient block above. Site, glazing tilt and this query form ONE
% consistent set and must be changed together, never one alone. The same POA
% series drives the PV array supplying the plant load, so generation and demand
% share a single resource -- which is the point of reporting a solar-driven SEC
% rather than a grid one.
%
% Month argument: 1-12 selects the calendar month; 13 returns zero irradiance
% (dark case, used for the storage and quasi-steady limit tests).
irr = pvgis_irradiance_data();
[FC.t_irr_data, FC.Gi_irr_data] = irr.get('Saudi', 3);    % 3 --> March

% Alternate site, retained for the sensitivity discussion. Switching it also
% requires FC.beta and the ambient block to change:
% [FC.t_irr_data, FC.Gi_irr_data] = irr.get('Dhaka', 3);

% ---- Return-value validation ---------------------------------------------
% Cheap, and it converts two otherwise SILENT failure modes into an error at the
% point of definition rather than a NaN two hours into a sweep: a mistyped field
% name creating a new field instead of setting an existing one, and an irradiance
% query that returned empty.
required = {'L','W','theta','beta','hcomp','Nx','Np','kappa_w','Ta','Vwind', ...
            'Tground','Tfeed','Vfeed','TDSfeed','Qfan','RH_amb','coil_mode', ...
            'T_coil','T_air_in','eta_air_heater','tfeed_dynamic', ...
            'dT_preheat_cap','dT_pinch_HX','Tfeed_target','dT_pinch_cond', ...
            'eps_recup','bl','t','stepSize','verbose', ...
            't_irr_data','Gi_irr_data'};
missing = required(~isfield(FC, required));
assert(isempty(missing), 'Baseline_Input_Variable:missingField', ...
       'Baseline is missing required field(s): %s', strjoin(missing,', '));

assert(~isempty(FC.t_irr_data) && ~isempty(FC.Gi_irr_data), ...
       'Baseline_Input_Variable:noIrradiance', ...
       'PVGIS query returned an empty series -- check pvgis_irradiance_data().');

% Feed preheat target must sit below the pinch ceiling, or the HX regime switches
% to pinch-limited and every result above the ceiling is inert for a reason
% unrelated to the design intent.
T_ceiling = FC.bl.T_sat_vap - FC.dT_pinch_HX;
assert(FC.Tfeed_target <= T_ceiling, 'Baseline_Input_Variable:pinch', ...
       'Tfeed_target = %.1f K exceeds the pinch ceiling %.1f K.', ...
       FC.Tfeed_target, T_ceiling);

assert(FC.T_coil >= FC.T_coil_floor, 'Baseline_Input_Variable:frost', ...
       'T_coil = %.1f K is below the frost guard %.1f K and the model has no frost logic.', ...
       FC.T_coil, FC.T_coil_floor);

end
