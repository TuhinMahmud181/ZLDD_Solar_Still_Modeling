function results = ZLDD_complete_modeling(FC)
% FULL_REPORT  Transient thermal-fluid-mass simulation of a multi-stage,
% tilted "ZLDD" (zero-liquid/zero-discharge-diverted) cascade solar still.
%
% results = Full_Report(FC) builds the model parameters from the input
% struct FC, assembles the initial condition, integrates the coupled
% PDE/ODE system forward in time with ode15s, extracts a fully
% post-processed results struct, and prints/plots the standard set of
% summary, mass-balance, energy-balance, and diagnostic reports.
%
% ---- PHYSICAL MODEL SUMMARY ----------------------------------------
% The still is a cascade of Np tilted absorber plates stacked one above
% the other, each Nx nodes long in the streamwise (downslope) direction,
% followed by ONE MORE film-bearing stage: the floor (stage/gap index
% Ns = Np+1). The floor carries its own water film exactly like a real
% plate -- fed by Plate Np's outlet near the bottom -- but with its own
% optical/thermal properties and footprint (floor_L x W, not L x W).
% Solar energy reaching the floor is absorbed by that water layer first
% (Beer-Lambert), then by the floor material itself, same as every other
% stage. Whatever survives the film as concentrated brine leaves the
% plant at the floor's own outlet. The compartment above the floor is an
% ACTIVE gap in the counter-current fan-air path -- the fan draws air in
% there, at the bottom of the stack.
%
% Saline feedwater flows down each plate as a thin film, is heated by
% solar radiation transmitted through the glass cover, and partially
% evaporates into a counter-current, fan-driven air stream that flows
% back up through the gaps between plates (including the gap above the
% floor) before exiting near the glass.
%
% Every derived quantity, energy accounting, and validity screen is
% computed and printed inside this file; the calling script supplies
% only the run controls in FC.
%
% ---- FLOWSHEET (single closed loop -- zero make-up, zero discharge) ---
%
%  WATER SIDE
%      seawater -> FEED PREHEAT HX -> CASCADE STILL -> brine
%                       ^                                   |
%                       | condenses, warms feed              v
%                       +-------------------------  BRINE EVAPORATOR
%                                                    (continuous, at
%                                                     NaCl saturation)
%                                                            |
%                                                     surplus vapour
%                                                     (feeds preheater,
%                                                      see AIR SIDE)
%
%  AIR SIDE  (counter-current, recirculating; every value below is a
%             model OUTPUT that moves with Qfan, T_coil, T_air_in,
%             eps_recup -- read from the report, never assume a number)
%
%      still exhaust                                   still inlet
%      Tv_exit, wv_exit                                T_air_in (setpt)
%           |                                               ^
%           v                                               |
%    +------------------ RECUPERATOR (passive) -------------------+
%    |  hot pass  ----------------------------------------->      |
%    |     Tv_exit -> T_recup_hot_out        (may drain condensate)
%    |                 heat crosses the wall                      |
%    |     T_coil  -> T_air_recup_out                              |
%    |  <-----------------------------------------  cold pass      |
%    +---------------------------------------------------------------+
%           |                                               ^
%           v  T_recup_hot_out, wv_to_coil                  | T_air_recup_out
%       [ AWG COIL ]                             [ VAPOUR-FIRED PREHEATER ] <- surplus
%        dries to T_coil                           lifts to T_int            vapour (evap.)
%        w = w_sat(T_coil)                                |
%           |                                               v
%           +------------------------------------> [ ELECTRIC AIR HEATER ]
%              T_coil, saturated                      trims T_int -> T_air_in
%              condensate |                                   |
%                         v                                   v
%                    PLASMA CHAMBER  <---------------   still bottom
%                    (terminal sink; also collects
%                     condensate from the recuperator
%                     hot pass and preheater steam side)
%
%  ORDERING CONSTRAINT (a consequence of the exchanger, not a coincidence):
%      T_coil < T_recup_hot_out   and   T_air_recup_out < Tv_exit
%  The cold outlet can never exceed the hot INLET; the hot outlet can
%  never fall below the cold INLET. A violation of either inequality
%  indicates an inconsistent recuperator closure, not an infeasible
%  operating point.
%
%  Cmin IS THE COLD SIDE. The hot pass condenses and sheds latent heat
%  at near-constant temperature (very large effective capacity rate), so
%  eps is defined on the cold (dry-air) side; the hot outlet is SOLVED
%  from the energy balance rather than assumed.
%
%  The two air-side heating units are specified DIFFERENTLY, by design:
%  the preheater takes whatever surplus steam exists (duty in, temp
%  out), while the electric heater holds the setpoint (temp in, duty
%  out). The heater therefore absorbs a steam supply that swings through
%  the day and falls to ZERO whenever the feed HX is supply-limited --
%  RATE the element for the FULL lift.
%
%  The loop is AIRTIGHT: no ambient air enters or leaves, so dry-air
%  flow is conserved around the circuit and the ATMOSPHERIC HARVEST
%  term is identically zero -- all product is seawater-derived. The
%  harvest term is still formed and reported, so that the split is
%  visible as an explicit zero rather than an omission.
%
%  ENERGY SUPPLY: all four electrical loads -- blower, brine
%  evaporator/crystalliser, electric air heater, and AWG compressor --
%  are supplied by on-site solar PV, so the plant is fully solar-driven
%  and off-grid. "Electric" denotes the form energy is delivered to the
%  unit, not a grid import. SEC is reported as purchased-work-equivalent
%  so it stays comparable with grid-connected RO/MED benchmarks.
%
%  Consequences worth remembering when reading results:
%    * results.mfw is the CASCADE EVAPORATION rate, not plant product.
%      Use results.product.mdot_total for output/SEC calculations.
%    * Inlet humidity is set by the AWG COIL (w_in = w_sat(T_coil)),
%      not by ambient; P.RH_amb is retained but unused by the air path.
%      Results are NOT comparable with once-through runs at the same
%      RH_amb.
%    * results.product.mdot_harvest is a structural zero (see above).
%    * The FEED PREHEAT HX creates a genuine RECYCLE (cascade -> brine
%      -> evaporator -> vapour -> HX -> feed) closed INSIDE the
%      integration (feed_preheat_T, called from rhs at k == 1) -- it's
%      negative feedback and stable, but it's inside the Jacobian.
%    * The air side carries NO recycle: T_coil and T_air_in are fixed
%      setpoints, so both inlet properties are known at build time.
%      Downstream units are evaluated once, after the PDE solve
%      (preheater_heater_awg).
%
%  SITE: all ambient/irradiance data are the Saudi Arabia (Rayed) March
%  dataset. The ambient block, glazing tilt, and PVGIS query are all on
%  this one site and are mutually consistent only for that site.
%
%==========================================================================
%
% States solved simultaneously (see build_parameters -> P.idx for the
% exact packing order inside the state vector Y):
%   Tg          - glass cover temperature                       (scalar)
%   Tv, wv      - vapor-gap air temperature & humidity ratio    (1 x Ns)
%   Twall       - side-wall temperature, one lumped node/gap    (1 x Ns)
%   M,  Tw      - film areal mass holdup & temperature          (Ns x Nx)
%   Ms          - film areal salt holdup                        (Ns x Nx)
%   Tp          - absorber-plate/floor temperature              (Ns x Nx)
%   u           - film velocity, a genuine momentum-PDE state   (Ns x Nx)
% (Ns = Np + 1: stage/gap Ns is the floor. The floor carries no separate
% lumped scalar node pair; its temperature is simply Tp(Ns,:) and its
% wall segment is Twall(Ns).)
%
% Film thickness (delta) and salt concentration (C) are NOT states -- they
% are algebraically recovered from the conserved holdups M, Ms via
% invert_holdup(), which iterates because water density depends on both.
%
% A few FC-provided or hard-coded values in build_parameters() are
% carried through the parameter structure without being consumed by any
% active equation (P.Tfan_in_seed, P.h_bottom_gap, P.RH_amb on the air
% path). Each is marked as such at its point of definition.


if nargin < 1, FC = struct(); end
t_start =tic;
default_stage = FC.plotstage;

% ---------------- 1. PARAMETERS ----------------------------------
P = build_parameters();

% ---------------- 2. INITIAL STATE --------------------------------
Y0 = build_initial_state(P);

% ---------------- 3. TIME INTEGRATION ------------------------------
tspan = linspace(0, P.t_sim,FC.stepSize);

if ~P.quiet
    fprintf('Estimating Jacobian sparsity pattern (one-time cost, speeds up the whole solve)...\n');
end
Jpat = estimate_jacobian_pattern(Y0, P);
if ~P.quiet
    fprintf('  Jacobian is %.2f%% dense (%d of %d possible entries) -- ode15s will only\n', ...
        100*nnz(Jpat)/numel(Jpat), nnz(Jpat), numel(Jpat));
    fprintf('  evaluate/factor the nonzero entries instead of treating it as fully dense.\n');
end

opts = odeset('RelTol',1e-6, ...
              'AbsTol', [1e-4, ...                     % Tg
                         1e-4*ones(1,P.Ns), ...        % Tv (Ns gaps, including the gap above the floor)
                         1e-8*ones(1,P.Ns), ...        % wv
                         1e-4*ones(1,P.Ns), ...        % Twall
                         1e-8*ones(1,P.n_field), ...   % M   (Ns stages x Nx, incl. floor film)
                         1e-4*ones(1,P.n_field), ...   % Tw
                         1e-8*ones(1,P.n_field), ...   % Ms
                         1e-4*ones(1,P.n_field), ...   % Tp  (Ns stages x Nx; Tp(Ns,:) is the floor temperature field)
                         1e-6*ones(1,P.n_field)]', ...  % u
              'Stats', ternary(P.quiet,'off','on'),...
              'NonNegative', find_nonneg_indices(P), ...
              'JPattern', Jpat, ...
              'Events', @(t,Y) zldd_events(t,Y,P));

if ~P.quiet
    opts = odeset(opts, 'OutputFcn', @(t,y,flag) progress_report(t,y,flag,P.t_sim));
    fprintf('Solving ZLDD solar still model with Np = %d plates, Nx = %d nodes/plate...\n', P.Np, P.Nx);
end


[tsol, Ysol, te, ye, ie] = ode15s(@(t,Y) rhs(t, Y, P), tspan, Y0, opts);

for jp = 1:numel(te)
    Sj = unpack_state(ye(jp,:).', P);
    [delta_j, C_j] = invert_holdup(Sj.M(:), Sj.Ms(:), Sj.Tw(:), P);
    [Cmax, idx] = max(C_j);
    [k_hit, x_hit] = ind2sub([P.Ns, P.Nx], idx);
    fprintf('Event %d at t=%.4f s: max C = %.2f kg/m3 at plate %d, node %d (delta=%.4e m)\n', ...
        ie(jp), te(jp), Cmax, k_hit, x_hit, delta_j(idx));
end


% ---------------- 4. EXTRACT ALL DATA -------------------------------
fprintf('Extracting full result set...\n');
results = extract_results(tsol, Ysol, P);

% ---------------- 5. DEFAULT RESULTS -----------------------------------

% Reporting is suppressed (P.quiet) when this file is driven from a
% parameter sweep, which would otherwise emit one full report and one
% figure set per run.
%
% Every printed figure carries its BASIS and its SOURCE beside the
% number: quantities that differ between blocks (still-only vs. plant,
% final-state vs. window-integral) differ because their basis differs.
if ~P.quiet
    print_headline(results);            % PART A
    print_stream_table(results);        % PART B
    print_procurement_spec(results);    % PART C
    check_baseline_validity(results);

    if P.verbose                        % PART D
        print_verification(results);
        % WINDOW SCALED TO THE SLOWEST PHYSICAL MODE, NOT HARD-CODED.
        % A fixed window cannot resolve anything, because the cascade
        % residence time moves with the geometry and the air flow, and
        % the salt field is the slowest mode by design. A convergence
        % check shorter than the transient it tests certifies nothing,
        % so the window is scaled to the solved residence time: two
        % residence times is the minimum that can see the salt field at
        % all, capped at half the run so the window never swallows the
        % startup transient.
        win_ss = min(max(2*results.retention_time_total, 100), 0.5*results.t(end));
        check_steady_state(results, default_stage, win_ss, 1e-3); % requested plot stage
    else
        fprintf(['\n  (Model-verification detail suppressed. Set FC.verbose = true\n' ...
                 '   for per-plate tables, closure residuals and solver audits.)\n']);
    end

    plot_performance(results);
    plot_baseline(results);  
end

results.runtime_total_s = toc(t_start);
fprintf('Total run time (solve + extract): %.2f s\n', results.runtime_total_s)
%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function P = build_parameters()

% ==================== FC Geometry ====================
P.L            = FC.L;                 % [m] plate length in streamwise (downslope) direction
P.W            = FC.W;                 % [m] plate width (cross-stream, assumed uniform for all Np plates)
P.theta        = deg2rad(FC.theta);    % [rad] plate tilt angle from horizontal (film flow driving angle)
P.beta         = deg2rad(FC.beta);     % [rad] glass cover tilt angle from horizontal 
P.Nx           = FC.Nx;                % [-] number of streamwise discretization nodes per plate 
P.Np           = FC.Np;                % [-] number of REAL absorber-plate stages
P.Ns           = P.Np + 1;             % [-] TOTAL film-bearing stages including the floor (stage Ns = the floor,
                                        %     modeled as a film-bearing plate with its own optical/thermal/
                                        %     geometric properties).
                                        %     Correspondingly there are Ns vapor gaps (gap k sits above stage k;
                                        %     gap 1 is under the glass, gap Ns is the compartment above the floor,
                                        %     an ACTIVE gap in the counter-current fan air path).




% ===================================================================
% MODEL CLOSURE SELECTION
% Each flag selects between alternative constitutive closures or
% numerical treatments. The defaults are the formulations used for all
% reported results; the alternatives are retained so that the
% sensitivity of any result to a given modelling choice can be
% quantified by re-running with a single flag changed.
% ===================================================================
P.opt.solar_film_model     = 'beer';    % 'beer'  = Beer-Lambert in-depth absorption
                                        % 'lumped'= additional surface absorptivity applied
P.opt.evap_clamp           = 'smooth';  % 'smooth' = C^inf softplus positive part
                                        % 'twoway' = signed driving force (condensation permitted)
                                        % 'hard'   = max(dP,0), non-differentiable at dP = 0
P.opt.evap_clamp_eps       = 1.0;       % [Pa] smoothing width for 'smooth'
P.opt.picard_iters         = 3;         % Picard sweeps for the local duct mass flow
P.opt.humidity_balance     = 'volume';  % 'volume' = V*rho_da*dw/dt = mdot_da*(win-w) + mevap
                                        % 'relax'  = first-order lag toward the same steady state
P.opt.tau_w_relax          = 1.0;       % [s] lag time constant used only by 'relax'
P.opt.nusselt_model        = 'smooth';  % 'smooth'   = blended laminar/turbulent, fixed Pr exponent
                                        % 'switched' = discontinuous Re = 2300 transition
P.opt.momentum_beta        = 6/5;       % Boussinesq momentum correction factor (1.0 = plug profile)
P.opt.interface_consistent = true;      % true  = conjugate pairs share a single nodal flux
                                        % false = each control volume evaluates its own lumped form
P.opt.route_solar_spill    = false;     % true  = aperture spill absorbed by wall node 1
                                        % false = tracked in the energy report only (absent from the dynamics)
P.opt.alpha_wall_solar     = 0.30;      % [-] wall solar absorptivity, used only if route_solar_spill
P.opt.alternating_plates   = true;      % true  = serpentine stack; consecutive plates run in
                                        %         opposite streamwise directions, so nodal
                                        %         vectors crossing a gap are index-reversed
                                        % false = all plates co-directional

% ----- Irradiance forcing mode (see solar_irradiance for the full note) -----
P.opt.constant_irradiance  = false;     % true  = Ig == P.G_const for all t. Makes the RHS
                                        %         AUTONOMOUS, so a genuine steady state exists
                                        %         and mfw_daily becomes an exact design-point
                                        %         evaluation rather than a rate extrapolation.
                                        % false = measured PVGIS record (non-autonomous).
                                        % This flag ALSO switches the steady-state criterion
                                        % (absolute drift vs. forcing-tracking ratio) and the
                                        % Esolar basis. Do not set it without also setting
                                        % P.G_const and P.t_sim below.
if isfield(FC,'constant_irradiance'), P.opt.constant_irradiance = FC.constant_irradiance; end

% Convergence criterion for the autonomous (constant-G) case. mfw_ss is
% averaged over the final frac_stationary_window of the run, and the
% least-squares drift of mfw across that window must fall below
% tol_stationary for the reported yield to be meaningful. Quote BOTH
% numbers in the manuscript: a yield without a stationarity residual is
% an unverified claim.
P.frac_stationary_window = 0.25;    % [-] trailing fraction of t_sim used as the stationary window
P.tol_stationary         = 1e-3;    % [-] max admissible relative drift of mfw over that window

if P.opt.constant_irradiance
    % Design-point plane-of-array irradiance. STATE WHICH G THIS IS in the
    % manuscript (annual-mean tilted-plane / design-day noon / seasonal
    % mean) and keep the choice consistent across all reported cases --
    % every headline yield scales with it.
    if isfield(FC,'G_const')
        P.G_const = FC.G_const;
    else
        % Energy-equivalent default: the mean of the measured record over
        % the nominal operating window ONLY (not the full 24 h record --
        % averaging over the night would halve G for no physical reason),
        % so that the constant-G case and the diurnal case receive the
        % SAME total incident energy. Any difference in daily distillate
        % between the two runs is then attributable purely to the
        % nonlinearity of mfw(G), which is the quantity worth reporting.
        t_win  = linspace(FC.t_day_start, ...
                          FC.t_day_start + FC.t_operating*3600, 2001);
        G_win  = max(interp1(FC.t_irr_data, FC.Gi_irr_data, ...
                             mod(t_win,86400), 'pchip', 0), 0);
        P.G_const = trapz(t_win, G_win) / (FC.t_operating*3600);
    end
else
    P.G_const = NaN;   % not consumed in diurnal mode
end

% -------------------------------------------------------------------
% DESIGN-PARAMETER DIAGNOSTIC NOTICES
% Each flag below GATES a MATLAB warning() issued at the END of
% build_parameters() -- i.e. BEFORE the Jacobian/solve, cluttering the
% console head with stack traces. All four are set FALSE so that the
% pre-solve console stays clean.
%
% Setting a flag true here emits the MATLAB warning() form for that
% condition; the conditions themselves are evaluated regardless, and the
% quantities behind them are printed in the reporting blocks (the
% floor-to-wall and floor-to-ground conductances in the gap diagnostics,
% the optical partition in the solar closure, the aperture spill in the
% solar accounting).
%
% Suppressing the warning() form is a presentation choice, NOT a
% resolution of the underlying conditions: every one is a disclosed
% modelling assumption and each belongs in the supplementary material.
% The checks are diagnostic only: no flag here enters any equation.
%
% Status of each condition:
%   warn_fan_preheat    - condition holds; reheat is genuinely recirculated.
%   warn_thermal_bridge - the floor-to-wall path carries the contact
%                         conductance and the lateral spreading resistance
%                         of the floor slab in series (P.R_spread_floor;
%                         see the construction of P.UA_floor_wall).
%                         Spreading dominates: R_floor_wallstrip alone is
%                         proportional to t_wall and, taken by itself,
%                         overstates the conductance by 3-4 orders of
%                         magnitude. The flag is left ACTIVE so the ratio is
%                         checked and reported on every run.
%   warn_plate_optics   - condition holds by design; plates are transmissive.
%   warn_aperture_spill - condition holds; spill is accounted, not routed.
% -------------------------------------------------------------------
P.opt.warn_fan_preheat   = false;  % Inlet air is reheated across the AWG condenser.
                                   % This enthalpy is recovered from within the air
                                   % loop rather than supplied from outside it, so it
                                   % is reported as a recirculated stream and is not
                                   % counted as an external energy input. The external
                                   % inputs are the absorbed solar radiation, the
                                   % blower work and the AWG compressor work.
P.opt.warn_thermal_bridge = false; % Conductance formed with the spreading resistance in series;
                                   % Floor-to-wallstrip conductance is orders of magnitude
                                   % larger than floor-to-ground, so the wall strip can
                                   % deliver more heat into the floor than the floor's own
                                   % solar absorption. Both conductances and the resulting
                                   % flux are printed in the gap diagnostics; read them
                                   % there rather than assuming a magnitude. The asymmetry
                                   % arises because R_floor_wallstrip = t_wall/k_wall carries NO
                                   % insulation in series, whereas R_floor_ground carries
                                   % t_ins/k_ins. The wall is a thin metallic skin
                                   % (P.t_wall = 2 mm at P.k_wall = 205 W/m-K), and because
                                   % R_floor_wallstrip is PROPORTIONAL to t_wall, a thin
                                   % skin gives a LARGER contact conductance than a thick
                                   % one: the bridge strengthens as the wall is thinned.
                                   % The lateral spreading resistance of the floor slab is
                                   % therefore what limits this path, and it is carried in
                                   % series in P.UA_floor_wall. The flag is suppressed to
                                   % keep console output clean; the ratio is still checked
                                   % and reported on every run.
P.opt.warn_plate_optics  = false;  % Cascade plates are transmissive acrylic rather
                                   % than absorbing surfaces. RF_p = 0.08 corresponds
                                   % to Fresnel reflection at two air/acrylic
                                   % interfaces and alpha_p = 0.01 to the weak bulk
                                   % absorption of PMMA in the solar band, giving a
                                   % per-stage beam survival of (1-a_w)*tau_p = 0.851.
P.opt.warn_aperture_spill = false; % The glass aperture exceeds the plate footprint,
                                   % so a fraction of the admitted beam falls outside
                                   % stage 1. This is quantified each run in the solar
                                   % closure block of print_energy_balance. Setting
                                   % P.opt.route_solar_spill = true deposits it on the
                                   % gap-1 wall node instead of accounting only.

% ----------------------------- FC Material Properties -----------------------------
P.kappa_w      = FC.kappa_w;           % [1/m] water absorption coefficient for in-depth solar absorption (Beer-Lambert decay through film) 

% ----------------------------- FC Ambient Conditions -----------------------------
P.Ta           = FC.Ta;                % [K] ambient (outdoor) air temperature. KELVIN, absolute. Enforced by the
                                       %     input-validation block; every radiative term (sigma*T^4) and every
                                       %     property correlation in this file assumes absolute temperature, so a
                                       %     degC value here fails loudly rather than degrading silently.
P.Vwind        = FC.Vwind;             % [m/s] ambient wind speed, used in the McAdams-type forced-convection correlation for hwind

% ----------------------------- FC Feed Conditions -----------------------------
P.Tfeed        = FC.Tfeed;             % [K] saline feedwater inlet temperature
P.Vfeed        = FC.Vfeed;             % [L] total feed volume supplied over the operating window 
P.TDSfeed      = FC.TDSfeed;           % [kg/m3] feed salinity (TDS). kg/m3 and g/L are numerically identical, so the
                                       %     two labels denote the same quantity; kg/m3 is used throughout for SI
                                       %     consistency with the salt holdup Ms [kg/m2] and concentration C = Ms/delta.
                                       %     Must be consistent with P.C_saturation, which uses the same units.

% ----------------------------- FC Air-Loop Conditions -----------------------------
% CLOSED AIR LOOP, ZERO MAKE-UP. The air stream is recirculated:
%
%     still -> AWG (dehumidify at T_coil) -> ELECTRIC AIR HEATER -> still
%
% The loop is airtight. No ambient air enters at any point and no loop
% air is discharged, so the DRY-AIR flow is strictly conserved around
% the circuit. Three consequences follow, and they determine how the
% inlet state must be specified:
%
%   (i)  Humidity ratio is set by the AWG COIL, not by ambient. Air
%        leaves the AWG saturated at T_coil and the air heater adds only
%        sensible heat, so the humidity ratio entering the cascade is
%              w_in = w_sat(T_coil).
%        With T_coil a FIXED setpoint this is computable once, at build
%        time, with no state dependence and no iteration -- which is the
%        entire reason the coil is specified rather than derived.
%
%   (ii) There is NO ATMOSPHERIC HARVEST. A closed loop with zero
%        make-up admits no atmospheric water at any point, so that term
%        is identically zero and ALL product is seawater-derived. The
%        harvest term is still formed and printed, so the seawater /
%        harvest split remains visible as an explicit zero.
%
%   (iii) P.RH_amb does not influence the air path. Ambient humidity
%        survives only in radiative/sky terms. It is retained as an
%        input, and is NOT consumed by the evaporative driving force.
%
% The still inlet AIR TEMPERATURE is set by the ELECTRIC AIR HEATER,
% which is specified by its OUTLET TEMPERATURE (P.T_air_in). Its duty is
% therefore the output:
%
%     Q_air_heater = mdot_da*(Cp_da + w*Cp_v)*(T_air_in - T_coil)
%
% This unit is electrically driven and is independent of the vapour
% raised in the brine evaporator: that vapour heats the SEAWATER FEED
% (see the feed preheat HX below), not the air.
P.Qfan         = FC.Qfan;              % [m3/s] fan volumetric flow rate (moist air, at fan inlet conditions)

% ---- DYNAMIC FEED PREHEAT ------------------------------------------
% true  : the still FEED temperature is computed inside the RHS from the
%         instantaneous cascade brine (see feed_preheat_T). It becomes a
%         time-varying OUTPUT and the vapour recycle is closed within
%         the integration. FC.Tfeed is then only a seed for the initial
%         condition and for pre-solve property evaluation.
% false : feed held at the constant value FC.Tfeed throughout.
%
% NOTE ON DIRECTION OF CAUSALITY. The feed preheat HX creates a genuine
% RECYCLE that the once-through configuration did not have:
%     cascade -> brine -> evaporator -> vapour -> HX -> feed -> cascade
% It is negative feedback (hotter feed -> more evaporation -> less brine
% -> less vapour -> cooler feed) and therefore stable, but it is inside
% the Jacobian and estimate_jacobian_pattern() must carry the coupling.
P.tfeed_dynamic = true;
if isfield(FC,'tfeed_dynamic'), P.tfeed_dynamic = logical(FC.tfeed_dynamic); end

P.dT_preheat_cap = 60;    % [K] STARTUP-ONLY cap on the computed FEED rise.
                          %     A SOLVER aid, not a physical limit: it
                          %     softens the first minutes, when the cascade
                          %     has barely evaporated, the brine is nearly
                          %     the full feed, and the implied vapour is at
                          %     its largest.
                          %     It is released after P.t_preheat_guard so
                          %     that the steady-state answer contains no
                          %     tuning constant. The PINCH is the physical
                          %     ceiling and is active at all times.
if isfield(FC,'dT_preheat_cap'), P.dT_preheat_cap = FC.dT_preheat_cap; end

P.t_preheat_guard = 1800; % [s] window over which the startup cap applies.
                          %     30 min: long enough to cover the initial
                          %     transient, short enough that it cannot
                          %     influence any reported result. If a run
                          %     reports CAPPED at t_end, this window is too
                          %     long -- or the cap still influences the solution.
if isfield(FC,'t_preheat_guard'), P.t_preheat_guard = FC.t_preheat_guard; end

P.dT_pinch_HX = 5;        % [K] minimum approach in the feed preheat HX.
                          %     Caps T_feed_out at T_sat_vap - dT_pinch,
                          %     i.e. the feed cannot be driven closer to
                          %     the condensing steam than this.
if isfield(FC,'dT_pinch_HX'), P.dT_pinch_HX = FC.dT_pinch_HX; end

P.Tfeed_target = 343;     % [K] DESIGN target feed temperature at the
                          %     cascade inlet. Deliberately well below the
                          %     pinch ceiling (T_sat_vap - dT_pinch_HX =
                          %     368.15 K) so the film does not enter near
                          %     flash conditions and psat_saline() is not
                          %     asked to work at the edge of its range.
                          %     The pinch remains an INACTIVE BACKSTOP.
                          %     Vapour the feed does NOT absorb is surplus
                          %     and is routed to the air preheater.
if isfield(FC,'Tfeed_target'), P.Tfeed_target = FC.Tfeed_target; end

P.dT_pinch_cond = 5;      % [K] minimum approach in the VAPOUR-FIRED AIR
                          %     PREHEATER. Caps the air at T_sat - 5 K.
                          %     Never binds at the current Qfan, but the
                          %     unit must not be able to misbehave if the
                          %     air flow is reduced.
if isfield(FC,'dT_pinch_cond'), P.dT_pinch_cond = FC.dT_pinch_cond; end

% ---- AIR-TO-AIR RECUPERATOR (still exhaust <-> AWG outlet) ----------
% Effectiveness of the passive plate exchanger that lets the loop
% exchange heat WITH ITSELF at the one point where the circuit holds two
% streams at different temperatures: air leaving the still (warm, humid)
% and air leaving the AWG coil (cold, dry). The loop order is unchanged
% -- still -> AWG -> heater -> still -- the duct simply passes through
% one box twice, on the way out and on the way back, separated by a wall.
%
% It cuts BOTH duties at once: the coil sees pre-cooled air, and the
% heater sees pre-warmed air. It touches only the SENSIBLE load; the
% latent duty still belongs entirely to the coil, because only the coil
% can take the air below its dew point.
%
% 0 disables the unit: the coil and heater then see unrecuperated air.
P.eps_recup = 0.75;       % [-] counterflow plate-exchanger effectiveness
if isfield(FC,'eps_recup'), P.eps_recup = FC.eps_recup; end
if P.eps_recup < 0 || P.eps_recup >= 1
    error('ZLDD:RecupEps', ...
        ['eps_recup = %.3f must lie in [0,1). Unity would require ' ...
         'infinite exchanger area.'], P.eps_recup);
end

% ---- ELECTRIC AIR HEATER (AWG outlet -> still inlet) ----------------
% SPECIFIED BY OUTLET TEMPERATURE. The duty follows.
P.T_air_in     = 310;                  % [K] still inlet air temperature setpoint
if isfield(FC,'T_air_in'), P.T_air_in = FC.T_air_in; end
P.eta_air_heater = 0.98;               % [-] resistance element efficiency
if isfield(FC,'eta_air_heater'), P.eta_air_heater = FC.eta_air_heater; end

P.Tfan_in      = P.T_air_in;           % [K] still inlet air temperature (alias read by
                                       %     the pre-solve property evaluations and by
                                       %     the initial condition).
P.Tfan_in_seed = P.T_air_in;           % [K] UNUSED. No active equation or report reads this
                                       %     value; it is carried for interface completeness.

P.COP_awg      = 3.5;                  % [-] constant-COP comparator, carried so that the
                                       %     lift-dependent model below can be shown
                                       %     side by side. NOT used for any reported figure.
% ---- Lift-dependent AWG coefficient of performance -------------------
% A constant COP discards the principal benefit of a warm coil: the
% heat-pump lift (T_cond - T_coil) collapses, so COP rises. The cycle is
% represented by a second-law efficiency against the Carnot cooling COP,
%
%     COP = eta_II * T_coil / (T_cond - T_coil),   T_cond = Tfan_in + dT_app
%
% eta_II is a mid-range value for a small vapour-compression unit and is
% an ASSUMPTION, not a result. It should be carried as a sensitivity
% band on every SEC figure rather than quoted as a single number, since
% the compressor term scales inversely with it.
P.eta_II_awg        = 0.45;            % [-] second-law efficiency of the AWG cycle
P.dT_cond_approach  = 5;               % [K] condenser approach above the reheat target
% ---- COIL TEMPERATURE MODE -------------------------------------------
% 'fixed' ONLY. In the CLOSED loop the coil sets the still inlet
% humidity ratio, w_in = w_sat(T_coil), so it feeds the cascade
% directly. 'auto' -- deriving T_coil from the SOLVED exhaust dew point
% -- would therefore be circular:
%       exhaust dew point -> T_coil -> w_in -> cascade -> exhaust
% and would require an outer fixed-point iteration over the whole solve.
% A fixed setpoint breaks that loop at build time for no cost, which is
% why the configuration specifies it. 'auto' is rejected rather than
% silently reinterpreted.
P.coil_mode      = 'fixed';
if isfield(FC,'coil_mode'),      P.coil_mode      = lower(FC.coil_mode); end
if ~strcmp(P.coil_mode,'fixed')
    error('ZLDD:CoilMode', ...
        ['coil_mode = ''%s'' is not available in the closed-loop\n' ...
         'configuration: T_coil sets the still inlet humidity, so\n' ...
         'deriving it from the solved exhaust is circular. Use\n' ...
         'FC.coil_mode = ''fixed'' and set FC.T_coil explicitly.'], P.coil_mode);
end
P.dT_coil_margin = 4;                  % [K] below the exhaust dew point
if isfield(FC,'dT_coil_margin'), P.dT_coil_margin = FC.dT_coil_margin; end
P.T_coil_floor   = 275;                % [K] frost guard -- the model has
                                       %     no frost logic, so 'auto' is
                                       %     clamped here with a warning.
if isfield(FC,'T_coil_floor'),   P.T_coil_floor   = FC.T_coil_floor; end

P.T_coil       = FC.T_coil;            % [K] AWG coil surface temperature. In the CLOSED loop this is a
                                       %     FIRST-ORDER INPUT, not a downstream quantity: loop air leaves the
                                       %     AWG saturated at T_coil and re-enters the cascade at that humidity
                                       %     ratio, so T_coil sets the evaporative driving force directly.
                                       %     Lowering it dries the air and raises (Psat(Tw) - Pv) throughout the
                                       %     cascade, at the cost of AWG duty and a larger heat-pump lift.
                                       %     Below freezing the coil frosts (the model has no frost logic).

% ---- Ambient air humidity ----
% RETAINED INPUT, NOT CONSUMED BY THE AIR PATH. The loop is closed and
% admits no ambient air, so RH_amb sets no inlet state and does NOT
% influence the evaporative driving force -- w_in comes from T_coil
% alone. It is kept because the radiative/sky terms and the ambient
% dew-point diagnostics reference it. Reading this variable as a driver
% of evaporation would misinterpret every sensitivity result in this
% file.
P.RH_amb = 0.60;                       % [-] ambient relative humidity (0-1)
if isfield(FC,'RH_amb') && ~isempty(FC.RH_amb), P.RH_amb = FC.RH_amb; end
if P.RH_amb <= 0 || P.RH_amb > 1
    error('ZLDD:BadRH','FC.RH_amb must lie in (0,1]. Got %g.', P.RH_amb);
end

% ===================================================================
% DOWNSTREAM UNITS (brine evaporator -> feed preheat HX -> AWG)
% -------------------------------------------------------------------
% BRINE EVAPORATOR. CONTINUOUS operation: concentrated brine from the
% floor stage enters, a crystallising slurry leaves at w_salt_target,
% and the liquor in the vessel sits permanently at NaCl saturation. The
% duty is COMPUTED from the inlet brine condition; it is not an input.
%
% THREE DISTINCT TEMPERATURES, not one. Conflating them is the single
% easiest error to make here, and it propagates straight into the feed
% preheat duty:
%
%   T_sat_vap  = Tsat(P_evap)          CONDENSING temperature in the HX
%   T_boil     = T_sat_vap + BPE       liquor boiling temperature
%   T_vap_heater = T_boil + dT_drive   element/jacket driving temperature
%
% The vapour leaving the evaporator is PURE STEAM. It carries a
% superheat of exactly BPE degrees, desuperheats in the first
% centimetres of the exchanger, and then condenses ISOTHERMALLY at
% T_sat_vap -- independent of brine concentration. So:
%
%   * BPE does NOT reach the feed preheat HX. The condensing temperature
%     is Tsat(P_evap) whether the liquor is at 5 wt% or 26 wt%.
%   * The BPE superheat is cp_v*BPE ~ 2.0*7 = 14 kJ/kg against
%     hfg = 2257 kJ/kg, i.e. ~0.6% of the duty.
%   * BPE is therefore a PURE HEATER-SIDE PENALTY: it raises the
%     electrical duty and the required element temperature and returns
%     nothing on the recovery side.
%
% WHY BPE IS A CONSTANT HERE. The inlet concentration is whatever the
% cascade delivers and is not fixed, but in CONTINUOUS operation the
% liquor is held at NaCl saturation (~26 wt%) at all times. Past
% saturation, removing water precipitates solid NaCl rather than raising
% the dissolved concentration: the liquid phase stays at X_sat, its
% water activity stays at ~0.75, and BPE PINS at ~7 K. The 0.99 target
% describes the SOLIDS FRACTION of the discharged slurry, not the
% concentration of the boiling liquor. There is consequently no
% BPE-versus-concentration trajectory to integrate.
% ===================================================================
P.bl.w_salt_target  = 0.99;      % [-] salt MASS FRACTION of the discharged slurry (SOLIDS, not liquor)
P.bl.P_evap         = 101325;    % [Pa] evaporator operating pressure
P.bl.T_sat_vap      = 373.15;    % [K] Tsat(P_evap) -- CONDENSING temperature in the feed preheat HX
P.bl.BPE            = 7.0;       % [K] boiling-point elevation of NaCl-saturated liquor (~26 wt%, a_w ~ 0.75)
P.bl.dT_drive       = 13.0;      % [K] element/jacket driving DT above the boiling liquor
P.bl.h_cryst        = 65e3;      % [J/kg salt] NaCl crystallisation enthalpy (released -- a small duty CREDIT)
P.bl.eta_heater     = 0.98;      % [-] electric resistance heater efficiency
P.bl.f_heatloss     = 0.15;      % [-] heater vessel heat loss as a fraction of useful duty
P.bl.eta_II_awg     = P.eta_II_awg;   % [-] AWG second-law efficiency
% Derived, after any FC override below:
%   T_boil       = T_sat_vap + BPE
%   T_vap_heater = T_boil + dT_drive

% Allow every downstream parameter to be overridden from the driver.
if isfield(FC,'bl') && isstruct(FC.bl)
    fn = fieldnames(FC.bl);
    for ii = 1:numel(fn)
        P.bl.(fn{ii}) = FC.bl.(fn{ii});
    end
end
if P.bl.w_salt_target <= 0 || P.bl.w_salt_target >= 1
    error('ZLDD:BadSaltTarget', ...
        'FC.bl.w_salt_target must lie strictly in (0,1). Got %g.', P.bl.w_salt_target);
end

% ---- Derive the evaporator temperature triple --------------------
% Order matters: these are built AFTER the FC override loop so that a
% driver overriding BPE, P_evap or dT_drive gets a consistent set rather
% than a stale T_boil. A driver may still override T_boil/T_vap_heater
% directly, in which case the explicit value wins.
if ~isfield(P.bl,'T_boil') || isempty(P.bl.T_boil)
    P.bl.T_boil = P.bl.T_sat_vap + P.bl.BPE;               % [K] boiling liquor
end
if ~isfield(P.bl,'T_vap_heater') || isempty(P.bl.T_vap_heater)
    P.bl.T_vap_heater = P.bl.T_boil + P.bl.dT_drive;       % [K] element/jacket
end

% The feed cannot be heated above the condensing steam, less the pinch.
% This is a HARD ceiling on the feed preheat and is enforced in
% feed_preheat_T(); stated here so the limit is visible with the
% parameters that set it.
P.T_feed_max = P.bl.T_sat_vap - P.dT_pinch_HX;             % [K]

% Console/plot suppression, used when this file is driven from a
% parameter sweep so that intermediate runs do not each emit a full
% report and figure set.
P.quiet = false;
if isfield(FC,'quiet'), P.quiet = logical(FC.quiet); end

% Verification depth. The default report (PARTS A-C) is the deliverable:
% headline performance, stream table, unit specifications. PART D is the
% solver audit -- per-plate tables, mass/energy closure, interface
% reciprocity, Reynolds regimes, gap properties, quasi-steady checks. It
% verifies the model rather than reporting it, so it is OFF by default.
% Set FC.verbose = true to emit it (e.g. when responding to a reviewer
% who asks how the model was validated).
P.verbose = false;
if isfield(FC,'verbose'), P.verbose = logical(FC.verbose); end

% ----------------------------- FC Simulation Parameters -----------------------------
P.t_sim        = FC.t;                 % [s] total simulation time span passed to the ODE integrator
P.t_operating  = FC.t_operating;       % [h] operating duration 

% ----------------------------- FC Daytime Simulation Window -----------------------------
P.t_day_start  = FC.t_day_start;       % [s] clock-time offset for start of simulation, used to align solar position/irradiance lookup with wall-clock time

% ----------------------------- FC PVGIS Irradiance Data -----------------------------
P.t_irr_data   = FC.t_irr_data;        % [ s] time vector for tabulated PVGIS irradiance data, used for interpolation at each RHS call
P.Gi_irr_data  = FC.Gi_irr_data;       % [W/m2] global tilted-plane irradiance data corresponding to t_irr_data, interpolated onto the solver's internal time steps

% ----------------------------- Geometry (hard-coded) -----------------------------
P.h_top_gap       = 0.10;       % [m] height of the top vapor-gap compartment (above the topmost plate, below glass)
P.h_bottom_gap    = 0.10;       % [m] INERT: no consumer anywhere in this file. The compartment above
                                % the floor takes its THERMAL/FLUID geometry from P.hcomp, the same
                                % height as every other inter-plate gap (see the P.hgap_narrow/
                                % P.hgap_wide construction below), so h_bottom_gap enters no gap volume,
                                % no wall area, no view factor, and no heat/mass transfer coefficient.
                                % The ENVELOPE calculation ('bottom_middle' -> P.stack_mean ->
                                % P.stack_front/stack_back, the stack heights printed by print_summary)
                                % also uses P.hcomp, so the reported geometry and the simulated
                                % geometry are the same geometry. The value is retained only as a
                                % record of the original design intent.
P.wall_clearance  = 0.10;       % [m] lateral clearance between plate edge and side wall, used to inflate chamber_L beyond plate length L

% ----------------------------- Optical Properties -----------------------------
P.tau_g           = 0.90;      % [-] glass transmissivity to solar-band radiation
P.alpha_g         = 0.05;      % [-] glass solar absorptivity
P.RF_g            = 0.05;      % [-] glass solar reflectivity 
P.eps_g           = 0.94;      % [-] glass long-wave (thermal-IR) emissivity, distinct from solar-band optical properties above

P.alpha_p         = 0.01;      % [-] absorber plate SOLAR absorptivity
P.RF_p            = 0.08;      % [-] absorber plate solar reflectivity
P.eps_p           = 0.90;      % [-] absorber plate long-wave emissivity (radiative exchange with glass/walls)

P.alpha_w         = 0.05;      % [-] water film solar absorptivity (governs in-depth absorption alongside kappa_w)
P.eps_w           = 0.95;      % [-] water film long-wave emissivity (surface radiative exchange, e.g. film-to-glass)

% ----------------------------- Glass Properties -----------------------------
P.mCp_g           = 6300;      % [J/K] lumped glass cover thermal mass (areal mCp already integrated over Ag, or per-unit-area)

% ----------------------------- Wall Properties -----------------------------
P.t_wall          = 2e-3;      % [m] side-wall structural material thickness
P.k_wall          = 205;       % [W/m-K] side-wall thermal conductivity (this value, ~205 W/m-K, is aluminum-range)
P.rho_wall        = 2700;      % [kg/m3] side-wall material density
P.Cp_wall         = 900;       % [J/kg-K] side-wall material specific heat
P.eps_wall        = 0.9;       % [-] side-wall external emissivity, for radiative loss to ambient/sky

% ----------------------------- Insulation Properties -----------------------------
P.t_ins           = 50e-3;      % [m] insulation layer thickness (outboard of wall material, in series for R_wall_out)
P.k_ins           = 0.03;      % [W/m-K] insulation thermal conductivity (fiberglass/foam range)
P.rho_ins         = 40;        % [kg/m3] insulation density
P.Cp_ins          = 1200;      % [J/kg-K] insulation specific heat

% ----------------------------- Absorber Plate Properties -----------------------------
P.rho_p           = 1190;      % [kg/m3] plate material density (this is polymer-range, e.g. PMMA/acrylic)
P.Cp_p            = 1470;      % [J/kg-K] plate specific heat
P.k_p             = 0.19;      % [W/m-K] plate through-thickness conductivity 
P.t_p             = 0.002;     % [m] plate thickness (2 mm)

% ----------------------------- Water Properties -----------------------------
P.k_water         = 0.61;      % [W/m-K] water thermal conductivity )

% ----------------------------- Heat Transfer Coefficients -----------------------------
P.hwind           = 5.7 + 3.8*P.Vwind;  % [W/m2-K] external forced-convection coefficient from wind (McAdams/Watmuff-type correlation) -- valid range typically Vwind < ~5 m/s; 

% ----------------------------- Ambient Boundary Conditions -----------------------------
P.Tground         = FC.Tground ;       % [K] effective ground/sky-adjacent temperature sink for the floor's 
P.Tref            = 273.15;    % [K] reference temperature for enthalpy datum (0 degC) 

% ----------------------------- Physical Constants -----------------------------
P.g               = 9.81;      % [m/s2] gravitational acceleration, drives film momentum PDE body force term
P.sigma           = 5.67e-8;   % [W/m2-K4] Stefan-Boltzmann constant, for all radiative exchange terms
P.Le              = 0.845;     % [-] Lewis number for water vapor in air (~0.845 is the standard value used in Chilton-Colburn analogy for evaporative mass transfer -> heat transfer coefficient conversion)
P.Rgas            = 8.314;     % [J/mol-K] universal gas constant
P.Mwater          = 0.018;     % [kg/mol] molar mass of water (18 g/mol), used in Clausius-Clapeyron / psychrometric relations

% ----------------------------- Base (Floor) Material Properties -- DIFFERENT from wall -----------------------------
P.rho_base   = 2700;      % [kg/m3] floor material density
P.Cp_base    = 900;       % [J/kg-K] floor material specific heat
P.t_base     = 2e-3;      % [m] floor thickness
% The floor's areal thermal mass is carried by the per-stage vectors
% (rho_solid_stage(Ns)*Cp_solid_stage(Ns)*t_solid_stage(Ns)) consumed by
% the d_Tp balance. No separate lumped areal capacitance is defined for
% the floor: a second, independently-computed copy of the same quantity
% could drift apart from it.



% ----------------------------- Air Thermophysical Properties -----------------------------
P.Cp_air          = air_cp(P.Ta);           % [J/kg-K] dry-air specific heat, evaluated at ambient Ta 
P.k_air           = air_conductivity(P.Ta); % [W/m-K] dry-air thermal conductivity at Ta
P.rho_air         = air_density(P.Ta);      % [kg/m3] dry-air density at Ta (ambient reference)
P.rho_air_fan_in  = air_density(P.Tfan_in); % [kg/m3] air density evaluated at 
P.mu_air          = air_viscosity(P.Ta);    % [Pa-s] dry-air dynamic viscosity at Ta

% ----------------------------- Water Thermophysical Properties -----------------------------
P.rhow_in         = water_density(P.Tfeed, P.TDSfeed);  % [kg/m3] feedwater density at feed T and TDS

% ----------------------------- Derived Geometry -----------------------------
P.dx              = P.L/P.Nx;                     % [m] streamwise node spacing (uniform mesh)

P.hcomp           = FC.hcomp;             % [m] vertical height of one compartment (plate-to-plate gap spacing)
                                           % Assigned here, ahead of every P.hcomp use below

P.chamber_L       = P.L + P.wall_clearance;   % [m] effective chamber length including wall clearance, distinct from plate length L used in flux terms
                    % NOTE: the cos(theta) term above projects the GLASS
                    % cover's own footprint (tilted at angle beta) onto the
                    % plate stack's geometry (tilted at theta) -- it is a
                    % glass-specific correction. The floor sits without that
                    % angle mismatch, so its own footprint length is simply
                    % P.floor_L below, NOT P.chamber_L.
P.floor_L         = P.L + P.wall_clearance;   % [m] floor's own footprint length (stage/gap Ns): plate length + wall clearance, no cos(theta) projection

% ---- DERIVED clearance actually available at the taper stations ----
% P.wall_clearance is a DESIGN INPUT: it is the slack that would exist if
% the plate lay flat (theta = 0). Once the plate is tilted, its horizontal
% projection shrinks to L*cos(theta), so the slack the chamber must absorb
% is the closure constraint
%
%       L*cos(theta) + (total clearance) = chamber_L
%
% which gives P.clearance_tot below. Splitting it symmetrically between the
% two ends puts the plate's back edge at horizontal station s = clearance_end
% and its front edge at s = chamber_L - clearance_end, so the plate spans
% exactly L*cos(theta) as it must.
%
% The taper station geometry below MUST use P.clearance_end, not
% P.wall_clearance. Spending P.wall_clearance at BOTH ends consumes more
% than the total clearance budget P.clearance_tot, compressing the plate's
% horizontal span to chamber_L - 2*wall_clearance instead of L*cos(theta).
% That error does not cancel: it displaces the x = L station along the
% glazing, and because h_m ~ 1/h_local it understates the plate-1 outlet
% evaporation at exactly the node where the flux is highest.
% P.wall_clearance itself is a per-end quantity and defines P.chamber_L
% and P.floor_L above; only the streamwise station map uses the split
% below.
P.clearance_tot   = P.chamber_L - P.L*cos(P.theta);   % [m] total slack, from the closure constraint
P.clearance_end   = P.clearance_tot / 2;              % [m] per-end slack, symmetric split
if P.clearance_tot < 0
    error('build_parameters:geometryInfeasible', ...
          ['Plate projection L*cos(theta) = %.4f m exceeds chamber_L = %.4f m. ' ...
           'Increase P.wall_clearance or reduce P.L / P.theta.'], ...
           P.L*cos(P.theta), P.chamber_L);
end

P.hgap_narrow  = zeros(P.Ns,1);   % "x=0 slot" -- see note below on gap 1 being reversed
P.hgap_wide    = zeros(P.Ns,1);   % "x=L slot" -- see note below on gap 1 being reversed
P.hwall_narrow = zeros(P.Ns,1);   % "x=0 slot" -- see note below on gap 1 being reversed
P.hwall_wide   = zeros(P.Ns,1);   % "x=L slot" -- see note below on gap 1 being reversed
% NOTE ON NAMING: despite the "narrow"/"wide" names, every downstream
% consumer of these four fields (film_gap_coeffs' local interpolation,
% its two duplicated copies in extract_results/print_energy_balance, and
% Awall_vec/V_gap_vec/tapered_rect_view_factor, which use these purely as
% an unordered {value1,value2} pair via sum/mean) treats hgap_narrow(k)/
% hwall_narrow(k) POSITIONALLY as "the value at x=0" and hgap_wide(k)/
% hwall_wide(k) as "the value at x=L" -- NOT by numeric magnitude. That
% matters because gap 1 (under the glass) tapers the OPPOSITE direction
% from every other gap: physically WIDE at x=0 (the feed-inlet end,
% where the fan/duct connection sits) and NARROW at x=L (the outlet
% end). This orientation is taken from the physical unit's duct layout and
% is an INPUT ASSERTION, not a model result -- nothing in this file can
% verify it, and reversing it changes the gap-1 transfer coefficients and
% view factors. It must be cited to the build drawing wherever it is
% relied upon in a manuscript. Every other gap
% (k=2..Ns, including the floor's own gap Ns) is narrow at x=0, wide at
% x=L. So for k==1 only, the two computed endpoint heights
% are deliberately assigned into the OPPOSITE slots from every other k:
% P.hgap_narrow(1) ends up holding the numerically LARGER value (the
% true x=0/wide height) and P.hgap_wide(1) the numerically SMALLER value
% (the true x=L/narrow height). This ordering is deliberate: the
% variable names describe which FORMULA/SLOT a value came from, not its
% magnitude, for gap 1 specifically.
%
% The assignments for k = 1 below therefore place the larger height in
% the x=0 slot and the smaller in the x=L slot, opposite to every other
% gap.

for k = 1:P.Ns
    if k == 1
        % Gap 1 lies between the inclined glazing (beta) and the first
        % cascade plate (theta). Because the glazing is steeply inclined
        % relative to the plate, the enclosure is deepest at the elevated
        % end of the chamber -- the end carrying the feed distributor and
        % the fan duct -- and shallows toward the opposite end. Since the
        % film on plate 1 originates at the feed distributor, this
        % corresponds to a maximum gap height at x = 0 and a minimum at
        % x = L, opposite in sense to every inter-plate gap.
        %
        % Because hc = Nu*k_a/Dh, the convective and mass transfer
        % coefficients in gap 1 are consequently smallest at the film
        % inlet and increase monotonically downstream. This produces a
        % streamwise evaporation profile on plate 1 that rises with x,
        % and, in combination with the film energy balance, an interior
        % maximum in the film temperature where the growing evaporative
        % load overtakes the heat supplied by the plate and the absorbed
        % beam.
        %
        % Only the streamwise distribution is affected by the taper
        % SENSE; the endpoint mean is unchanged under a sense reversal, so
        % Awall_vec, V_gap_vec, P.Acomp and P.Dh are independent of the
        % sense. They are NOT independent of the endpoint VALUES: any
        % change to the clearance convention (see P.clearance_end above)
        % moves both endpoints and therefore the mean, which propagates
        % into Dh -> Nu -> hc/hm for the whole gap. Do not edit the two
        % assignments below without re-checking Acomp/Dh.
        %
        % Convention check on the plate term: the plate's horizontal run is
        % L*cos(theta), so its vertical drop is L*cos(theta)*tan(theta) =
        % L*sin(theta). The x=L expression is therefore already written in
        % the same horizontal-run convention as the clearance_end*tan(theta)
        % term at x=0 -- the two forms look different but are consistent.
        offset = P.h_top_gap;
        P.hgap_narrow(k)   = offset + tan(P.beta)*(P.chamber_L - P.clearance_end) + P.clearance_end*tan(P.theta);    % x = 0 : WIDE (feed + fan end)
        P.hgap_wide(k)     = offset + tan(P.beta)*P.clearance_end + P.L*sin(P.theta);                                % x = L : NARROW
        P.hwall_narrow(k)  = offset + P.chamber_L*tan(P.beta);   % x = 0 : tall wall at the back
        P.hwall_wide(k)    = offset;                              % x = L : short wall at the front
    else
       
        % Same closure constraint as gap 1: the plate's back edge sits at
        % horizontal station s = clearance_end and its front edge at
        % s = chamber_L - clearance_end, spanning L*cos(theta).
        offset = P.hcomp;
        slope  = 2*tan(P.theta);
        P.hgap_narrow(k)   = offset + slope*P.clearance_end;
        P.hgap_wide(k)     = offset + slope*(P.chamber_L - P.clearance_end);
        P.hwall_narrow(k)  = offset ;
        P.hwall_wide(k)    = offset+ slope*P.chamber_L;
    end
end

% Terminal-stage endpoint heights. The floor stage spans floor_L rather
% than L and is bounded by the chamber shell rather than by a further
% plate, so its endpoint heights carry no downstream taper contribution.
% The terminal entries are therefore adjusted once, after the stage loop.
% (These use P.clearance_end for the same reason as the stage loop above:
%  the taper contribution must be removed at the SAME horizontal stations
%  the loop used, or the subtraction will not cancel.)
P.hgap_narrow(end)  = P.hgap_narrow(end)  - tan(P.theta)*P.clearance_end;
P.hgap_wide(end)    = P.hgap_wide(end)    - tan(P.theta)*(P.chamber_L - P.clearance_end);
P.hwall_narrow(end) = P.hwall_narrow(end);   % offset only: no taper term to remove at x=0
P.hwall_wide(end)   = P.hwall_wide(end)   - tan(P.theta)*P.chamber_L;
% ---- Physical envelope (REPORTING ONLY -- no dynamic consequence) ----
% The three P.stack_* scalars below are the external build height of the
% assembly, printed by print_summary. They enter NO balance equation:
% every heat/mass transfer path uses the per-gap heights P.hgap_* and
% P.hwall_* instead. They exist so the buildability of a given parameter
% set can be checked against the reported performance.
%
% bottom_middle: mid-length height of the bottom compartment, i.e. the
%   taper rise accumulated from the floor's narrow end to mid-span, plus
%   the compartment height P.hcomp -- the SAME height the simulation uses
%   for the gap above the floor (P.h_bottom_gap enters nothing here; see
%   the parameter-consistency note below). It is additive, so hcomp maps
%   one-for-one onto the reported stack height and the envelope can be
%   checked for buildability against the geometry that produced the
%   performance figures.
% top_moddle:  mid-length height of the top compartment. The
%   (tan(beta)+tan(theta)) sum carries BOTH tilt contributions because the
%   glass and the plate stack are inclined at different angles, so the
%   compartment between them opens at the sum of the two slopes.
%   (This local is spelled 'top_moddle' throughout; it is a reporting
%   variable only.)
%
% stack_mean sums the mid-span heights of the interior gaps (indices
% 2:end-1, since the top and bottom compartments are supplied separately
% by top_moddle/bottom_middle rather than being double counted) and adds
% the two end compartments. stack_front/stack_back then apply the +/- half
% of the glass tilt rise across chamber_L, giving the low and high edges
% of the enclosure respectively.
% PARAMETER CONSISTENCY. bottom_middle uses P.hcomp, NOT P.h_bottom_gap.
% The compartment above the floor is sized by P.hcomp exactly like every
% other inter-plate gap, and h_bottom_gap enters no thermal/fluid
% geometry at all. Using it here would make the REPORTED
% physical envelope (stack_mean -> stack_front/stack_back, printed by
% print_summary) depend on a different parameter from the one the model
% simulates, so the stack height quoted in a manuscript would not
% correspond to the geometry that produced the performance figures.
% Buildability of the stack height is itself a design constraint.
% NOTE: this is the same horizontal-station geometry as the taper loop, so
% it uses P.clearance_end for the same reason (see the derivation above).
% REPORTING ONLY -- this affects the printed stack height and enters no
% balance equation.
bottom_middle = (P.chamber_L/2-P.clearance_end)*tan(P.theta)+P.hcomp;
top_moddle    = (tan(P.beta)+tan(P.theta))*P.chamber_L/2 + P.h_top_gap ;

P.stack_mean = sum(P.hwall_narrow(2:end-1)+P.hwall_wide(2:end-1))/2+bottom_middle+top_moddle ;
P.stack_back = P.stack_mean +P.chamber_L*tan(P.beta)/2 ;
P.stack_front = P.stack_mean-P.chamber_L*tan(P.beta)/2 ;
% ---- Reference/diagnostic scalars only -- NOT used for Re/Nu/hc ----

h_mean  = (P.hgap_narrow+P.hgap_wide)/2 ;
P.Acomp = P.W*h_mean; 
P.Dh    = 2*P.W*h_mean./(P.W+h_mean) ;

P.fanArea = 0.25*pi*FC.fandia^2 ;

P.Awall_first     = (P.chamber_L + P.W) * ...
                    (P.chamber_L*tan(P.beta) + 2*P.h_top_gap)-P.fanArea;

% ---- Awall_vec: per-k wall surface area, using mean of hwall_narrow/hwall_wide ----
% Exact, not approximate: 2*(chamber_L+W)*h is linear in h, so evaluating at
% the mean h equals the true integral of the tapered wall area along x.
P.Awall_vec = zeros(P.Ns,1);
for k = 1:P.Ns
    h_wall_mean_k  = (P.hwall_narrow(k) + P.hwall_wide(k)) / 2;
    P.Awall_vec(k) = 2*(P.chamber_L + P.W) * h_wall_mean_k;
end
if P.Ns >= 2
    P.Awall_rest  = P.Awall_vec(2);   % scalar alias; identical across k>1 since offset=hcomp is constant
else
    P.Awall_rest  = P.Awall_vec(1);   % degenerate Np=0 case: the floor is the ONLY stage, there is no "k>1" gap to reference
end

P.Ag              = P.chamber_L*P.W/cos(P.beta);    % [m2] glass cover area
P.Ap              = P.L*P.W;                        % [m2] single (real) plate's active area

% APERTURE AND TARGET AREA
% The beam admitted through the glazing carries power tau_g*Ig*Ag, while
% the transmitted flux is applied to stage 1 as an areal flux over the
% plate footprint Ap = L*W. Because the glazing aperture
% Ag = chamber_L*W/cos(beta) exceeds Ap, a fraction of the admitted beam
% falls beyond the plate edge onto the chamber shell. P.A_spill
% quantifies this area and the corresponding power is reported
% separately by print_energy_balance so that the solar accounting closes
% on the aperture. Setting P.opt.route_solar_spill = true deposits this
% power on the gap-1 wall node rather than accounting for it alone.
P.A_spill         = max(P.Ag - P.Ap, 0);            % [m2] admitted-but-untargeted aperture area
P.spill_frac      = P.A_spill / P.Ag;               % [-]  fraction of admitted solar with no target
P.Ap_floor        = P.floor_L * P.W;                % [m2] floor stage's footprint area (different length than a plate: floor_L = L+clearance, not chamber_L -- see note above)

% ---- Ap_stage / L_stage / dx_stage: per-stage footprint area, streamwise
% length, and node spacing (Ns x 1). Stages 1..Np are the ordinary
% absorber plates (area Ap, length L); stage Ns is the floor, whose
% footprint/length differ (Ap_floor, floor_L). ----
P.Ap_stage = P.Ap * ones(P.Ns,1);      P.Ap_stage(P.Ns) = P.Ap_floor;
P.L_stage  = P.L  * ones(P.Ns,1);      P.L_stage(P.Ns)  = P.floor_L;
P.dx_stage = P.L_stage / P.Nx;         % [m] per-stage streamwise node spacing

% ---- Per-stage solid material properties (Ns x 1): plates 1..Np use the
% absorber-plate properties; stage Ns (the floor) uses its own, different
% optical/thermal properties. ----
P.alpha_solid_stage = P.alpha_p * ones(P.Ns,1);   % alpha_solid_stage(Ns) assigned further below, once P.alpha_base is defined
P.RF_solid_stage    = P.RF_p    * ones(P.Ns,1);   % RF_solid_stage(Ns) IS used -- it sets q_refl at the floor in solar_split(). It is
                                                  % overwritten with P.RF_base further below. The quantity with no consumer at the
                                                  % terminal stage is IT_next, not RF (see the k==Ns branch in rhs()).
P.rho_solid_stage   = P.rho_p   * ones(P.Ns,1);
P.Cp_solid_stage    = P.Cp_p    * ones(P.Ns,1);
P.t_solid_stage     = P.t_p     * ones(P.Ns,1);
P.k_solid_stage     = P.k_p     * ones(P.Ns,1);
% (rho_base/Cp_base/t_base/k_base assigned into the Ns-th slot further
% below, once those floor-material properties are defined.)

% ---- V_gap_vec: per-k vapor-gap volume, using mean of hgap_narrow/hgap_wide ----
% Exact (not approximate), same reasoning as Awall_vec: volume = A_mean*W*L_gap
% is linear in h, so evaluating at the mean h equals the true integral along x.
% L_gap differs for the top compartment (chamber_L, under the glass -- a
% genuinely different length owing to the glass/plate-stack tilt-angle
% mismatch), the bottom compartment (floor_L, above the floor,
% matching the floor's own footprint length -- no such tilt mismatch),
% and the ordinary inter-plate compartments (L, plate-to-plate span).
P.V_gap_vec = zeros(P.Ns,1);
for k = 1:P.Ns
    if k == 1
        L_gap_k = P.chamber_L;
    elseif k == P.Ns
        L_gap_k = P.floor_L;
    else
        L_gap_k = P.L;
    end
    h_gap_mean_k   = (P.hgap_narrow(k) + P.hgap_wide(k)) / 2;
    P.V_gap_vec(k) = h_gap_mean_k * P.W * L_gap_k;
end

% ----------------------------- Derived Flow Parameters -----------------------------
P.mdot_a          = P.rho_air_fan_in * P.Qfan;                           % [kg/s] total moist-air mass flow rate from the fan 
% Still inlet humidity ratio -- CLOSED LOOP, SET BY THE COIL.
%   Loop air leaves the AWG SATURATED at T_coil, then passes through the
%   electric air heater, which adds only sensible heat. Sensible heating
%   conserves the humidity ratio, so the humidity ratio entering the
%   cascade is the coil saturation value and the still inlet RH is the
%   derived quantity (it falls as the air is heated at constant w):
%
%       w_in = w_sat(T_coil) = 0.622*Psat(T_coil)/(Patm - Psat(T_coil))
%
%   This is what propagates into the PDE solve: the evaporative driving
%   force in every gap is (Psat(Tw) - Pv), and Pv is built from this w.
%   The humidity ratio is coil-derived rather than ambient-derived, so
%   results are NOT comparable with once-through runs at the same
%   RH_amb.
%
%   Because T_coil is a FIXED setpoint this evaluates once, here, with
%   no state dependence -- the recycle closes at build time.
P.wv_fan_in       = humidity_ratio_from_RH(1.0, P.T_coil, P);            % [kg_vapor/kg_dry-air] saturated at the coil
P.wv_amb          = humidity_ratio_from_RH(P.RH_amb, P.Ta, P);           % [kg/kg] ambient humidity ratio (diagnostics only -- no ambient air enters the loop)
P.phi_fan_in      = P.wv_fan_in / humidity_ratio_from_RH(1.0, P.Tfan_in, P);  % [-] resulting RH at the still inlet (derived, not imposed)
P.Pv_fan_in       = P.wv_fan_in * 101325 / (0.622 + P.wv_fan_in);        % [Pa] vapor partial pressure entering the cascade (Patm as in vapor_partial_pressure) 
P.mdot_da         = P.mdot_a / (1 + P.wv_fan_in);                        % [kg/s] DRY-air mass flow 
P.mdot_vapor_in   = P.mdot_a - P.mdot_da;                                % [kg/s] vapor mass flow entering with the fan stream 

P.Gamma_feed      = (P.Vfeed*P.rhow_in/1000) / ...
                    (P.W*P.t_operating*3600);
P.mdot_feed_water = (P.Vfeed*P.rhow_in/1000 - P.Vfeed*P.TDSfeed/1000) ...
                    / (P.t_operating*3600);   % [kg/s] WATER in the feed,
                    % used by feed_preheat_T to get the brine flow from
                    % the instantaneous water balance.
                   

% ---- Feed preheat HX terminals ------------------------------------
% P.Tfeed_in  : seawater as SUPPLIED, i.e. the COLD side of the HX. This
%               is the genuine boundary condition of the plant.
% P.Tfeed     : seawater as SUPPLIED. IDENTICAL to P.Tfeed_in and it
%               STAYS that way -- it is never overwritten with the
%               solved hot-side value after the solve, because every
%               feed MASS expression (P.rhow_in, Gamma_feed, the
%               recovery denominator) is defined at the supply state.
%               The solved cascade inlet temperature lives in
%               results.feedpre.Tfeed_out / P.Tfeed_cascade_in.
%               With tfeed_dynamic it is a time-varying OUTPUT computed
%               in feed_preheat_T(); the value set here is only a seed
%               for the initial condition and for pre-solve property
%               evaluation. Seeding it AT the inlet temperature is
%               deliberate and physically right: at t = 0 no brine has
%               been raised, the evaporator has no duty, and the feed
%               genuinely does enter cold at sunrise, warming through
%               the morning as the cascade develops.
P.Tfeed_in        = P.Tfeed;                                             % [K] seawater as supplied (cold side)

% ---- Design target consistency ------------------------------------
% Placed here, not at the parameter declaration, because both guards
% read P.Tfeed_in and it does not exist until the line above.
if P.Tfeed_target <= P.Tfeed_in
    error('ZLDD:FeedTargetLow', ...
        ['Tfeed_target = %.1f K must exceed the supplied seawater ' ...
         'temperature %.1f K.'], P.Tfeed_target, P.Tfeed_in);
end
if P.Tfeed_target > P.bl.T_sat_vap - P.dT_pinch_HX
    warning('ZLDD:FeedTargetAbovePinch', ...
        ['Tfeed_target = %.1f K exceeds the pinch ceiling %.1f K; the ' ...
         'pinch will bind instead and the target is inoperative.'], ...
         P.Tfeed_target, P.bl.T_sat_vap - P.dT_pinch_HX);
end
% Hot side: the temperature at which the film actually ENTERS the
% cascade. SEEDED HERE at the supply value so the field always exists,
% then overwritten after the solve when tfeed_dynamic is on. Consumers
% that need the HOT side read THIS, never P.Tfeed -- P.Tfeed is the
% supply state and every feed MASS expression depends on it staying so.
%
% The seed matters for ordering, not just for safety: compute_mass_balance
% runs inside extract_results, BEFORE the post-solve reconciliation, and
% reads its own copy via results.P (published at results.P = P). Without
% the seed that copy has no such field and the run dies mid-extraction.
P.Tfeed_cascade_in = P.Tfeed;                                            % [K] cascade inlet (hot side)
P.mdot_feed_total = (P.Vfeed*P.rhow_in/1000) / (P.t_operating*3600);     % [kg/s] TOTAL feed (water + salt)

P.mdot_salt_feed  = (P.Vfeed/1000 * P.TDSfeed) / ...
                    (P.t_operating*3600);
                    % [same mass units as TDSfeed, per s] total salt mass input rate,

% ----------------------------- Derived Thermal Parameters -----------------------------
P.R_wall_out      = P.t_wall/P.k_wall + ...
                    P.t_ins/P.k_ins + ...
                    1/P.hwind;

P.mCp_wall_areal  = P.rho_wall*P.Cp_wall*P.t_wall + ...
                    P.rho_ins*P.Cp_ins*P.t_ins;
                    % [J/m2-K] areal thermal mass of the wall+insulation composite
                    % NOTE: Twall(Ns), the wall strip around the gap above the floor,
                    % uses this SAME generic areal capacitance, via wall_rhs(), exactly
                    % like every other Twall(j). The wall strip is therefore treated as
                    % the same wall+insulation composite everywhere in the stack, with no
                    % separate structural-material-only capacitance anywhere.

% ----------------------------- Time Parameters -----------------------------
P.t_clock_offset  = P.t_day_start;   % [s] alias of t_day_start

% ----------------------------- Floor <-> wall strip: conductive contact, not convective film -----------------------------
% Contact area = perimeter * base depth (wall material wraps the 10 cm
% floor edge). This conduction path (floor bulk -> wrapping wall strip)
% is applied as a lumped
% (spatially-averaged) loss/gain term inside the floor stage's (k==Ns)
% own energy balance in rhs(), and correspondingly as a gain inside
% wall_rhs() for Twall(Ns). See rhs() and wall_rhs() for the exact terms.
P.Acontact_floor_wall = 2*(P.floor_L + P.W) * P.t_base;
                    % [m2] contact area between floor edge and wrapping wall strip,
                    
P.R_floor_wallstrip = P.t_wall / P.k_wall;   % [m2-K/W] conduction resistance through the wall's thickness

% ---------- Lateral (spreading) resistance inside the floor slab ----------
% WHY A SPREADING RESISTANCE IS REQUIRED. Carrying only the
% through-thickness conduction of the wall skin,
%     UA_fw = Acontact_floor_wall / (t_wall/k_wall),
% and smearing the resulting heat flow uniformly over the floor footprint
% Ap_floor implicitly asserts ZERO lateral resistance between the floor
% perimeter and the floor interior. That is not a small error:
%
%   - R_floor_wallstrip is PROPORTIONAL to t_wall, so making the wall
%     thinner makes the bridge STRONGER, not weaker. At t_wall = 2 mm the
%     path carries UA_fw = O(1e4) W/K against UA_floor_ground = O(0.3) W/K
%     -- a ratio of order 1e4. At that conductance the floor is a perfect
%     thermal short to Twall(Ns): the temperature difference collapses,
%     Q = UA*dT becomes numerically indeterminate, and the floor stage's
%     energy balance (hence its evaporation contribution) is an artifact.
%
%   - The floor temperature field is resolved in x only. Transverse
%     (W-direction) conduction from the perimeter inward is absent from
%     the discretisation entirely, so it must be supplied as a lumped
%     series resistance or it is simply missing.
%
% Model: the floor slab is a 2-D conductor of thickness t_base, held at
% the wall-strip temperature around its ENTIRE perimeter (the wall
% material wraps all four edges, per Acontact_floor_wall = 2*(L+W)*t_base)
% and exchanging heat with the film over its whole area. The
% area-averaged temperature rise of such a slab defines the spreading
% resistance exactly. Writing T = (q'''/k) u with -grad^2 u = 1, u = 0 on
% the perimeter, the double sine series gives the closed form
%
%   <u> = sum_{m,n odd} 64 / ( pi^6 m^2 n^2 [ (m/a)^2 + (n/b)^2 ] )
%   R_spread = <u> / ( k_base * t_base * a * b )                    [K/W]
%
% with a = floor_L, b = W. This is EXACT for the stated boundary
% condition and carries no aspect-ratio restriction. A 1-D two-edge
% surrogate, which would assume a long slab fed only from the two long
% edges, overstates R_spread whenever the floor is not long and narrow --
% i.e. whenever W/floor_L is of order one.
% For a square the series converges to <u>/(ab) = 0.035144, so
% R_spread -> 0.035144/(k*t), independent of size.
% ---------- Floor material properties & the floor <-> ground path ----------
% These sit AHEAD of the spreading-resistance construction below, which
% consumes P.k_base.
P.k_base     = 205.;     % [W/m-K] thermal conductivity of the floor material
P.alpha_base = 0.95;     % [-] floor solar absorptivity
P.R_floor_ground = P.t_base/P.k_base + P.t_ins/P.k_ins;   % [m2-K/W] floor bulk -> ground:
                                                          % conduction through the floor slab
                                                          % in series with the insulation layer

nser   = 1:2:201;                              % odd modes only; even modes vanish
[mm, nn] = ndgrid(nser, nser);
u_mean = sum(sum( 64 ./ (pi^6 .* mm.^2 .* nn.^2 .* ...
                  ((mm./P.floor_L).^2 + (nn./P.W).^2)) ));
P.R_spread_floor = u_mean / (P.k_base * P.t_base * P.floor_L * P.W);   % [K/W]

% Total floor->wallstrip conductance: contact conduction IN SERIES with
% lateral spreading. NOTE this is a UA [W/K], not an areal resistance --
% rhs() and every diagnostic must consume it as such.
P.UA_floor_wall = 1 / ( P.R_floor_wallstrip/P.Acontact_floor_wall + P.R_spread_floor );   % [W/K]
P.UA_floor_grnd = P.Ap_floor / P.R_floor_ground;                                          % [W/K]

% The series form carries no aspect-ratio restriction, so no W/floor_L
% guard is required. What DOES remain an idealisation is the
% perimeter Dirichlet condition: it assumes the wall strip holds the
% whole floor edge at Twall(Ns). That is reasonable while the contact
% resistance is small compared with R_spread_floor (verify from the
% console: the printed contact-only UA should be >> UA_floor_wall), and
% it is the conservative direction -- it maximises the bridge.

% FLOOR OUTWARD FACE -- NO RADIATIVE BRANCH. The model carries NO
% radiative exchange between the floor's outward face and the
% ground/sky: that boundary is purely conductive through
% R_floor_ground = t_base/k_base + t_ins/k_ins. Behind 50 mm of
% k = 0.03 W/m-K insulation the conductive path dominates, so this is a
% defensible simplification, but it is a STATED assumption. If the floor
% is ever exposed (uninsulated underside, elevated stand), the radiative
% branch must be added explicitly.

% FLOOR SOLAR REFLECTIVITY
% The absorbing floor terminates the beam path. Although no stage lies
% below it, its reflectivity remains an active optical property: it
% partitions the arriving beam between absorption and reflection back
% into the enclosure. The floor is treated as opaque, so absorptivity
% and reflectivity sum to unity.
P.RF_base         = 1 - P.alpha_base;   % [-] floor solar reflectivity (opaque closure)

% ---- Wire the floor's own material properties into the per-stage vectors
% (k_solid_stage needs P.k_base, only just defined above) ----
P.rho_solid_stage(P.Ns)   = P.rho_base;
P.Cp_solid_stage(P.Ns)    = P.Cp_base;
P.t_solid_stage(P.Ns)     = P.t_base;
P.k_solid_stage(P.Ns)     = P.k_base;
P.alpha_solid_stage(P.Ns) = P.alpha_base;
P.RF_solid_stage(P.Ns)    = P.RF_base;    % opaque closure: alpha + RF = 1

% ---- Per-stage optical closure validation ----
% Each stage must satisfy alpha + RF <= 1, the remainder being the
% transmitted fraction. Violation would imply a stage disposing of more
% radiant power than it receives.
for kk = 1:P.Ns
    tot_k = P.alpha_solid_stage(kk) + P.RF_solid_stage(kk);
    if tot_k > 1 + 1e-12
        error('ZLDD:OpticalClosure', ...
            ['Stage %d: alpha (%.3f) + RF (%.3f) = %.3f > 1. The stage ' ...
             'disposes of more solar than reaches it. Fix the optical ' ...
             'properties before running.'], kk, ...
             P.alpha_solid_stage(kk), P.RF_solid_stage(kk), tot_k);
    end
end

% ---- Alias for reporting/diagnostic functions that refer to
% "the bottom-gap wall area": this is gap Ns's wall area. ----
P.Awall_bottom_gap = P.Awall_vec(P.Ns);   % [m2]
% ----------------------------- Gap laminar Nusselt closure -----------------------------
% Fully developed laminar Nu for each vapour gap, set by the two-surface
% wall-flux ratio r = q_top/q_bottom via Nu = 5.385/(1 - 0.3461*r).
% See the header of nusselt_gap() for the derivation, the anchors it
% reproduces, and the quasi-steady assumption implied by fixing r.
% Values below follow from r evaluated on a converged solution:
%   gap 1  r = -0.095  (glazing above, radiating to sky)      -> 5.21
%   gap 2  r = -0.577  (plate 1 underside, coolest solid)     -> 4.49
%   gap >=3 r ~ +1     (symmetric: T_plate_above ~ T_film)    -> 7.541
% 7.541 is the both-walls-ISOTHERMAL value, chosen over the both-walls-
% uniform-FLUX value 8.235 because the plates are metal with resolved
% axial conduction and explicit surface temperature fields.
P.Nu_lam_default = 7.541;                            % interior symmetric gaps
P.Nu_lam_gap     = P.Nu_lam_default * ones(P.Ns,1);  % [-] one per gap
P.Nu_lam_gap(1)  = 5.214;
if P.Ns >= 2
    P.Nu_lam_gap(2) = 4.489;
end

% ----------------------------- Physical safeguards -----------------------------
P.delta_dryout       = 5e-6;    % [m] film-thickness threshold (5 micron) below which dry-out logic begins to engage/ramp
P.delta_dryout_full  = 1.5e-6;  % [m] film-thickness threshold below which dry-out is treated as complete/full 
P.C_saturation       = 317;     % [kg/m3 of SOLUTION] NaCl saturation, on the same basis as C = Ms/delta.
                                % BASIS. C is salt mass per unit volume of SOLUTION, not per unit volume of water.
                                % Two independent confirmations: (i) invert_holdup() sets delta = M/rho with rho the
                                % BRINE density, so delta is solution volume per unit area; (ii) water_density(T,C)
                                % forms S = C/rho and feeds it to a correlation whose salinity argument is kg salt
                                % per kg SOLUTION. Both fix the basis as solution-volumetric.
                                % VALUE. Saturated NaCl at 25 C is 35.98 g per 100 g water, i.e. a solution mass
                                % fraction w = 0.2646. With rho(298.15 K, C) from water_density() the fixed point
                                % C = w*rho(T,C) converges to 317.1 kg/m3, at rho = 1198.5 kg/m3 (which reproduces
                                % the tabulated saturated-brine density to within 0.04%).
                                % TEMPERATURE. On THIS basis the limit is nearly flat: 318 kg/m3 at 0 C rising to
                                % 328 kg/m3 at 100 C (3.4% over the full range), because rising solubility and
                                % falling brine density partly cancel. Holding it constant is therefore a better
                                % approximation here than it would be on a per-water basis, where the same range
                                % spans 356 -> 391 kg/m3 (9.7%). Using the 25 C value is conservative for a hot film.
                                % CAVEAT. water_density() is a seawater-range correlation (fitted to S <~ 0.12 kg/kg)
                                % evaluated at S ~ 0.265 -- a 2x extrapolation. It agrees with tabulated data at 25 C,
                                % but that is one point of agreement, not validation across the range.
                                % This threshold gates the saturation/dry-out logic and hence the terminal-stage
                                % behaviour, so it is a first-order parameter for any zero-liquid-discharge claim.
P.C_saturation_ramp_start = 310; % [same basis as C_saturation] reporting threshold only -- there is no ramp in the
                                % RHS; the sole consumer is a [WARN] line in the verification block. Set to 98% of
                                % C_saturation, i.e. the warning fires once the film is within 2% of saturation.

% ---- State layout bookkeeping ----
% Per-plate 2D fields (Np x Nx): M, Tw, Ms, Tp, u. M and Ms are the
% conserved mass/salt holdup states; u is a genuine momentum-PDE state,
% not an algebraic Nusselt closure (see the momentum equation in rhs).
% Lumped per-gap (Np): Tv, wv
% Scalar: Tg
% Per-stage 2D fields are (Ns x Nx): stage Ns is the floor, carrying the
% same M/Tw/Ms/Tp/u field set as every real absorber plate (see P.Ns
% above). Likewise Tv/wv/Twall span Ns gaps, gap Ns being the gap above
% the floor. The floor carries no separate scalar node pair -- its
% temperature is Tp(Ns,:) and its wall segment is Twall(Ns).
P.n_field  = P.Ns*P.Nx;
P.idx.Tg   = 1;
P.idx.Tv   = P.idx.Tg + 1                  : P.idx.Tg + P.Ns;
P.idx.wv   = P.idx.Tv(end) + 1             : P.idx.Tv(end) + P.Ns;
P.idx.Twall = P.idx.wv(end) + 1             : P.idx.wv(end) + P.Ns;
base = P.idx.Twall(end);
P.idx.M     = base + (1:P.n_field);            base = base + P.n_field;
P.idx.Tw    = base + (1:P.n_field);            base = base + P.n_field;
P.idx.Ms    = base + (1:P.n_field);            base = base + P.n_field;
P.idx.Tp    = base + (1:P.n_field);            base = base + P.n_field;
P.idx.u     = base + (1:P.n_field);            base = base + P.n_field;
P.nstate    = base;

% ---- Gap radiative-enclosure geometry ----
% Precompute per-gap view factors (film / top-surface / wall enclosure)
% ONCE here -- pure geometry, depends only on P.L, P.W, P.hcomp, P.theta,
% P.beta, P.h_top_gap, so it automatically updates whenever those are
% changed, without needing separate bookkeeping. See compute_gap_view_factors.
P = compute_gap_view_factors(P);

% ===================================================================
% DESIGN-PARAMETER DIAGNOSTIC NOTICES
% Each check below identifies a parameter whose value materially
% determines a reported performance figure, and reports the magnitude of
% its effect so that the choice can be assessed against the physical
% configuration being represented.
% ===================================================================
% ---- Air-path consistency (CLOSED LOOP) ----
% The coil sets the still inlet humidity ratio, w_in = w_sat(T_coil), so
% it is an upstream quantity again. Two conditions follow.
%
% HARD: the air heater must heat, not cool. Loop air leaves the coil at
% T_coil and the heater raises it to T_air_in, so T_air_in > T_coil is
% structural -- violating it would mean the "heater" removes energy,
% which no unit in this flowsheet can do.
if P.T_air_in <= P.T_coil
    error('ZLDD:AirHeaterInverted', ...
        ['T_air_in = %.1f K must exceed T_coil = %.1f K. Loop air ' ...
         'leaves the AWG at the coil temperature and the electric air ' ...
         'heater raises it to the still inlet setpoint; the stated ' ...
         'values would require that heater to REMOVE %.1f K.'], ...
         P.T_air_in, P.T_coil, P.T_coil - P.T_air_in);
end
%
% SOFT: the coil must lie below the still EXHAUST dew point or the AWG
% recovers nothing from the air stream. That cannot be evaluated from
% the inlet state and is checked after the solve in
% check_baseline_validity. A necessary condition IS available up front:
% the exhaust is the inlet air plus whatever the cascade adds, so the
% exhaust dew point is at least
% T_coil. A coil close to the still inlet temperature therefore leaves
% little humidification room.
if P.T_coil >= P.T_air_in - 2
    warning('ZLDD:CoilHigh', ...
        ['T_coil = %.1f K sits within 2 K of the still inlet air ' ...
         'temperature (%.1f K). The loop air enters nearly saturated, ' ...
         'so the evaporative driving force in the cascade will be ' ...
         'small. Lower T_coil or raise T_air_in.'], ...
         P.T_coil, P.T_air_in);
end
if P.T_coil < 273.15
    warning('ZLDD:CoilFrosting', ...
        ['T_coil = %.1f K is below the freezing point. A practical ' ...
         'evaporator coil frosts at this temperature and the ' ...
         'humidity ratio would be limited by ice rather than by the ' ...
         'saturation line used here.'], P.T_coil);
end

if P.opt.warn_fan_preheat && abs(P.T_air_in - P.T_coil) > 1
    Q_fan_th_preview = P.mdot_da * moist_air_cp(P.T_air_in, P.wv_fan_in) ...
                       * (P.T_air_in - P.T_coil);
    warning('ZLDD:FanThermalInput', ...
        ['Still inlet air enters %.1f K above the AWG coil (%.1f K), ' ...
         'carrying %.0f W. The datum is the COIL, not ambient: the loop ' ...
         'is closed, so this air was not drawn from outside. The heat ' ...
         'is supplied by the ELECTRIC AIR HEATER and is a genuine ' ...
         'external input, appearing in the system energy balance.'], ...
         P.T_air_in - P.T_coil, P.T_coil, Q_fan_th_preview);
end

% Thermal-bridge check, evaluated on the series conductance (contact +
% lateral spreading). The contact-only value is reported alongside it so
% the weight of the spreading resistance is visible in the console record
% rather than buried in this file.
UA_floor_wall_contact_only = P.Acontact_floor_wall / P.R_floor_wallstrip;
if P.opt.warn_thermal_bridge && P.UA_floor_wall > 20*P.UA_floor_grnd
    warning('ZLDD:FloorWallThermalBridge', ...
        ['Floor->wallstrip conductance (%.1f W/K, incl. spreading) ' ...
         'exceeds floor->ground (%.2f W/K) by %.0fx. Contact-only value ' ...
         'is %.0f W/K (spreading reduces it by a factor of %.0f). ' ...
         'The wall-strip path carries no insulation in series while the ' ...
         'ground path carries %.0f mm of it; the skin is t_wall = %.1f mm ' ...
         'at k_wall = %.0f W/m-K.'], ...
         P.UA_floor_wall, P.UA_floor_grnd, P.UA_floor_wall/P.UA_floor_grnd, ...
         UA_floor_wall_contact_only, ...
         UA_floor_wall_contact_only/P.UA_floor_wall, ...
         1000*P.t_ins, 1000*P.t_wall, P.k_wall);
end

if P.opt.warn_plate_optics && P.alpha_p < 0.5 && P.alpha_p + P.RF_p < 0.5
    warning('ZLDD:PlateTransmissive', ...
        ['alpha_p = %.3f with RF_p = %.3f implies plate transmissivity ' ...
         'tau_p = %.3f -- i.e. the cascade plates are modelled as ' ...
         'SEMI-TRANSPARENT, passing %.0f%% of the solar to the stage below. ' ...
         'Combined with rho_p = %.0f kg/m3 and k_p = %.2f W/m-K (polymer ' ...
         'range) this is self-consistent for acrylic plates over an ' ...
         'absorbing floor, but is NOT an absorber-plate cascade. If the ' ...
         'plates are coated metal, alpha_p/RF_p/rho_p/Cp_p/k_p must ALL ' ...
         'be changed together.'], ...
         P.alpha_p, P.RF_p, 1-P.alpha_p-P.RF_p, 100*(1-P.alpha_p-P.RF_p), ...
         P.rho_p, P.k_p);
end

if P.opt.warn_aperture_spill && P.spill_frac > 0.02
    warning('ZLDD:ApertureSpill', ...
        ['%.1f%% of the solar admitted through the glass aperture (Ag = ' ...
         '%.3f m2) lands outside plate 1 (Ap = %.3f m2) and has no target ' ...
         'stage. It is now tracked explicitly (see the energy report). ' ...
         'Set P.opt.route_solar_spill = true to deposit it on the gap-1 ' ...
         'wall node instead of only accounting for it.'], ...
         100*P.spill_frac, P.Ag, P.Ap);
end

end 



function P = compute_gap_view_factors(P)
% Computes the gray-diffuse view factors for the 3-surface radiative
% enclosure (film / top-surface / wall) formed by each gap in the
% cascade. This is a purely geometric calculation -- it depends only on
% fixed dimensions and tilt angles, not on temperature -- and is
% therefore called once from build_parameters(), not on every rhs()
% evaluation.
%
% ENCLOSURE DEFINITION (per gap k = 1..Np):
%   "f"   = film k                           (bottom)  area = P.Ap
%   "top" = plate (k-1), or the glass if k=1 (top)      area = P.Ap for
%                                                        k>1, P.Ag for k=1
%   "w"   = lumped side wall k                (sides)   area = Awall_vec(k)
%
% METHOD: for each gap, the vertical separation between the film and the
% top surface is treated as varying linearly along the flow direction,
% from gap_narrow at the feed end (x=0) to gap_wide at the fan/opening
% end (x=L_gap). The view factor for this tapered geometry is obtained
% from tapered_rect_view_factor(), then the remaining five view factors
% of the enclosure are obtained algebraically from reciprocity and
% summation, without further approximation.
%
% ASSUMPTIONS (taper geometry differs between the two branches below):
%   k = 1 (top compartment, under the glass): the taper is driven by the
%   difference between the glass tilt angle (P.beta, e.g. 21 deg) and
%   plate 1's tilt angle (P.theta, e.g. 3 deg), with the extraction fan
%   located at the far (x=L_gap) end. The gap is assumed to be narrowest
%   at the feed end, equal to P.h_top_gap, and to widen linearly toward
%   the fan end. The length scale used here is P.chamber_L (the glass's
%   own footprint length, P.L + P.wall_clearance), matching the length
%   already implicit in P.Ag = P.chamber_L*P.W/cos(P.beta) -- using P.L
%   here instead would assume a smaller top surface than P.Ag actually
%   represents. This taper model should be checked against the physical
%   glass mounting if higher precision is required.
%
%   k > 1 (between stacked plates in the folded cascade): the taper is
%   driven by the alternating tilt (+theta / -theta) of consecutive
%   plates, so the relative tilt between any two adjacent plates is
%   2*theta regardless of position in the stack. The gap is assumed to
%   taper linearly over the plate length P.L, narrow at the feed end and
%   wide at the opening end, consistent with the clearance needed for
%   water and vapor to pass to the next plate.
%
%   Both tapers are converted to a view factor via a 2D (per-unit-width)
%   crossed-strings correction applied to Hottel's closed-form flat
%   parallel-rectangle result (evaluated at the mean gap height). This
%   multiplicative combination is an ENGINEERING APPROXIMATION, exact
%   only in the limit W/L -> infinity, where the 3D problem reduces to
%   the 2D case by construction. It has not been validated here against
%   a numerical (quadrature or Monte Carlo) reference; such a comparison
%   is recommended before quoting the radiative results as high-accuracy,
%   particularly if W is not much larger than L or the mean gap height
%   for the geometry in use.
%
% INPUT:
%   P : parameter structure; must already contain P.Np, P.L, P.W,
%       P.theta, P.beta, P.hcomp, P.h_top_gap, P.chamber_L, P.Ag, P.Ap,
%       and the fields consumed by wall_area_per_gap().
%
% OUTPUT (all added to P, one entry per gap k = 1..Np):
%   P.Ffp, P.Ffw               : film's view factor to top, to wall  [-]
%   P.Fpf, P.Fpw               : top's view factor to film, to wall  [-]
%   P.Fwf, P.Fwp, P.Fww        : wall's view factor to film, to top,
%                                 to itself (the wall ring is concave
%                                 and can see itself)                [-]
%   P.Atop_vec                 : top-surface area used for gap k     [m^2]

Awall_vec = wall_area_per_gap(P);
Ns = P.Ns;  W = P.W;
% No single shared length scale here: the k=1 (glass, chamber_L), k=Ns
% (floor's own footprint, floor_L), and 1<k<Ns (plate-to-plate, length L)
% gaps all have different physical top/film lengths, so each branch below
% sets its own L_gap.

P.Ffp = zeros(Ns,1);  P.Ffw = zeros(Ns,1);
P.Fpf = zeros(Ns,1);  P.Fpw = zeros(Ns,1);
P.Fwf = zeros(Ns,1);  P.Fwp = zeros(Ns,1);  P.Fww = zeros(Ns,1);
P.Atop_vec = zeros(Ns,1);

for k = 1:Ns
    if k == 1
        % ---- Top compartment: glass (tilt beta) over plate 1 (tilt theta) ----
        L_gap = P.chamber_L;                 % glass's own footprint length
        Atop  = P.Ag;
    elseif k == Ns
        % ---- Bottom compartment: plate Np's underside over the floor's
        % own film (the floor's footprint uses floor_L = L+clearance, NOT
        % chamber_L -- unlike the glass, the floor has no tilt-angle
        % mismatch to project through). The "top" surface here is still
        % an ordinary plate (Plate Np), so Atop = Ap as usual. ----
        L_gap = P.floor_L;
        Atop  = P.Ap;
    else
        % ---- Inter-plate gap: alternating +theta / -theta stack ----
        L_gap = P.L;                          % plate-to-plate span
        Atop  = P.Ap;
    end

    % ---- Taper read directly from the shared P.hgap_narrow/P.hgap_wide
    % fields (built once in build_parameters) rather than re-deriving the
    % narrow/wide formula inline, so this loop and film_gap_coeffs consume
    % one definition of the geometry.
    gap_narrow = P.hgap_narrow(k);
    gap_wide   = P.hgap_wide(k);
    
    Af = P.Ap_stage(k);   % film area for THIS gap's lower surface (Ap for k=1..Np, Ap_floor for k=Ns)
    Aw = Awall_vec(k);

    Ffp = tapered_rect_view_factor(L_gap, W, gap_narrow, gap_wide, k);

    % ---- Remaining view factors, obtained from summation and reciprocity ----
    Ffw = 1 - Ffp;                   % film is flat: sees only "top" and "wall"
    Fpf = (Af   * Ffp) / Atop;       % reciprocity; Af may differ from Atop (e.g. k=1, k=Ns)
    Fpw = 1 - Fpf;                   % top is flat: sees only "film" and "wall"
    Fwf = (Af   * Ffw) / Aw;         % reciprocity
    Fwp = (Atop * Fpw) / Aw;         % reciprocity
    Fww = 1 - Fwf - Fwp;             % wall is concave and sees itself

    P.Ffp(k) = Ffp;  P.Ffw(k) = Ffw;
    P.Fpf(k) = Fpf;  P.Fpw(k) = Fpw;
    P.Fwf(k) = Fwf;  P.Fwp(k) = Fwp;  P.Fww(k) = Fww;
    P.Atop_vec(k) = Atop;
end

end % compute_gap_view_factors


function Ffp = tapered_rect_view_factor(L, W, gap_narrow, gap_wide, k)
% Approximates the view factor between two directly opposed L-by-W
% rectangles whose separation varies linearly along L, from gap_narrow
% at one end to gap_wide at the other.
%
% METHOD: Hottel's closed-form flat parallel-rectangle result is
% evaluated at the mean gap height, then scaled by the ratio of two
% per-unit-width (2D) crossed-strings view factors -- one for the
% tapered geometry, one for the equivalent flat geometry at the same
% mean height. The formula collapses exactly to the standard flat
% parallel-plate result when gap_narrow = gap_wide: the two crossed
% diagonals are then equal, so the tapered and flat 2D view factors
% coincide and the scaling ratio is unity.
%
% APPROXIMATION: this multiplicative combination of a 3D flat baseline
% with a 2D taper correction is exact only as W/L -> infinity, the limit
% in which the 3D and 2D problems coincide. For finite W it is an
% interpolation, not a derived identity, and has not been checked here
% against a numerical reference. Its accuracy should be treated with
% more caution as W approaches L or the mean gap height in magnitude.
%
% INPUTS:
%   L, W                  : rectangle length and width [m]
%   gap_narrow, gap_wide  : separation at the two ends of the taper [m]
%   k (optional)          : gap index, used only to identify the source
%                            of the diagnostic warning below; omit or
%                            pass NaN if not applicable.
%
% OUTPUT:
%   Ffp : view factor from the "film" surface to the "top" surface [-],
%         clipped to the physical range [0, 1].

if nargin < 5, k = NaN; end

h_avg = (gap_narrow + gap_wide)/2;

% ---- 2D (per-unit-width) taper correction, via Hottel's crossed strings ----
AD = sqrt(L^2 + gap_wide^2);
BC = sqrt(L^2 + gap_narrow^2);
F2D_taper = (AD + BC - gap_narrow - gap_wide) / (2*L);
F2D_flat  = (sqrt(L^2 + h_avg^2) - h_avg) / L;

% ---- 3D flat parallel-rectangle baseline (Hottel's chart formula), at the mean gap ----
X = L/h_avg;  Y = W/h_avg;
F3D_flat = (2/(pi*X*Y)) * ( ...
      log( sqrt((1+X^2)*(1+Y^2)/(1+X^2+Y^2)) ) ...
    + X*sqrt(1+Y^2)*atan(X/sqrt(1+Y^2)) ...
    + Y*sqrt(1+X^2)*atan(Y/sqrt(1+X^2)) ...
    - X*atan(X) - Y*atan(Y) );

% ---- combine: scale the 3D baseline by the 2D taper-to-flat ratio ----
Ffp = F3D_flat * (F2D_taper / F2D_flat);

% A raw result outside [0,1] indicates the approximation is breaking
% down for this geometry (e.g. an aggressive taper relative to L or W);
% this is reported rather than silently clamped away.
if Ffp < -1e-6 || Ffp > 1 + 1e-6
    warning('tapered_rect_view_factor:outOfBounds', ...
        ['Gap %d: raw Ffp = %.4f is outside the physical range [0,1] ' ...
         'before clipping (L=%.3f, W=%.3f, gap_narrow=%.4f, gap_wide=%.4f). ' ...
         'The taper approximation may not be valid for this geometry.'], ...
        k, Ffp, L, W, gap_narrow, gap_wide);
end

Ffp = min(max(Ffp, 0), 1);
end 




function [Qf_gap, Qtop_gap, Qwall_gap] = radiosity_all_gaps(S, P)
% Solves the 3-surface gray-diffuse radiosity network (film / top /
% wall) for every gap in the cascade, using the view factors precomputed
% by compute_gap_view_factors(). The vapor phase is treated as
% non-participating in radiation, consistent with every other radiative
% term in this model; vapor still exchanges heat with all three
% surfaces, but only convectively, through the hc_wv and wall_rhs terms
% defined elsewhere.
%
% SIGN CONVENTION: Qf_gap(k), Qtop_gap(k), and Qwall_gap(k) are each the
% net radiative heat LOSS from that surface (positive = net emitter).
% By construction, Qf_gap(k) + Qtop_gap(k) + Qwall_gap(k) = 0 for every
% gap, to machine precision -- this follows directly from the
% reciprocity relations already enforced in compute_gap_view_factors(),
% and can be used as a hard assertion in testing rather than as an
% approximate check.
%
% ASSUMPTION -- spatial lumping: the film and top-surface temperatures
% entering the radiosity solution are each a single spatial mean over
% the Nx axial nodes (see Tf and Ttop below), even though both fields
% are Nx-resolved everywhere else in the model (conduction, evaporation,
% and solar absorption all act node-by-node). The resulting single
% Q-value per gap is then applied uniformly across all Nx nodes by the
% calling code. Consequently, any axial temperature gradient within a
% gap is not visible to the radiative exchange. This is a genuine
% simplification of the model, not only a convenience, and bears on the
% axial convergence behaviour of the film fields.
%
% CROSS-REFERENCE: Qtop_gap(k) is later used, as Qtop_gap(k+1), as
% plate k's own radiative loss term in its energy balance -- i.e. plate
% k acts as the top surface of the enclosure for gap (k+1). This index
% shift (gap k corresponds to plate (k-1) as its top surface) is the
% convention the plate energy balance assumes.
%
% INPUTS:
%   S : state structure; must contain Tw (film temperature field,
%       Np x Nx), Tp (plate temperature field, Np x Nx), Tg (glass
%       temperature, scalar), and Twall (per-gap wall temperature,
%       Np x 1)
%   P : parameter structure; must contain the view factors and areas
%       produced by compute_gap_view_factors(), plus P.eps_w, P.eps_g,
%       P.eps_p, P.eps_wall, and P.sigma
%
% OUTPUTS (one entry per gap k = 1..Np):
%   Qf_gap(k)    : net radiative loss from film k                  [W]
%   Qtop_gap(k)  : net radiative loss from gap k's top surface      [W]
%                   (the glass if k=1, otherwise plate (k-1))
%   Qwall_gap(k) : net radiative loss from wall k                  [W]

Ns = P.Ns;
Awall_vec = wall_area_per_gap(P);

Qf_gap    = zeros(Ns,1);
Qtop_gap  = zeros(Ns,1);
Qwall_gap = zeros(Ns,1);

for k = 1:Ns
    Tf = mean(S.Tw(k,:));   % spatial mean over Nx -- see lumping assumption above (k=Ns: the floor's own film)

    if k == 1
        Ttop    = S.Tg;
        eps_top = P.eps_g;
        Atop    = P.Ag;
    else
        % k = 2..Ns: top surface is always the underlying REAL plate
        % (k-1), including for k=Ns where the top surface of the bottom
        % gap is Plate Np's underside -- an ordinary plate, so eps_p/Ap
        % apply to it as to any other plate.
        Ttop    = mean(S.Tp(k-1,:));   % spatial mean over Nx -- see lumping assumption above
        eps_top = P.eps_p;
        Atop    = P.Ap;
    end

    [Qf_gap(k), Qtop_gap(k), Qwall_gap(k)] = solve_radiosity3( ...
        Tf, Ttop, S.Twall(k), ...
        P.eps_w, eps_top, P.eps_wall, ...
        P.Ffp(k), P.Ffw(k), P.Fpf(k), P.Fpw(k), P.Fwf(k), P.Fwp(k), P.Fww(k), ...
        P.Ap_stage(k), Atop, Awall_vec(k), P.sigma);
end
end 


function [Qf, Qp, Qw] = solve_radiosity3(Tf, Tp, Tw, eps_f, eps_p, eps_w, ...
                                          Ffp, Ffw, Fpf, Fpw, Fwf, Fwp, Fww, ...
                                          Af, Ap, Aw, sigma)
% Solves the gray-diffuse radiosity network for a closed 3-surface
% enclosure. Surfaces: f (film), p (top: the plate above, or the glass
% for the topmost gap), w (lumped side wall).
%
% GOVERNING EQUATIONS: for a gray-diffuse enclosure, the radiosity of
% each surface i satisfies
%     J_i = eps_i*Eb_i + (1 - eps_i) * sum_j F_ij * J_j
% which is rearranged below into the linear system A*J = b. The
% self-view terms F_ff and F_pp are zero, since the film and top
% surfaces are both flat and are simply omitted from those two rows;
% F_ww is retained on the diagonal of the wall row, since the wall ring
% is concave and can see itself. The term is required for the energy
% balance below to close; it vanishes only for a convex enclosure.
%
% SIGN CONVENTION: Q_i > 0 means surface i is a net emitter, losing
% radiative energy on balance to the other two surfaces.
%
% INPUTS:
%   Tf, Tp, Tw          : surface temperatures [K]
%   eps_f, eps_p, eps_w : surface emissivities [-]
%   Ffp, Ffw, Fpf, Fpw,
%   Fwf, Fwp, Fww       : view factors between the three surfaces [-];
%                          must already satisfy reciprocity (e.g.
%                          Af*Ffp = Ap*Fpf), or the exact conservation
%                          property below will not hold
%   Af, Ap, Aw          : surface areas [m^2]
%   sigma               : Stefan-Boltzmann constant [W/m^2/K^4]
%
% OUTPUTS:
%   Qf, Qp, Qw : net radiative heat loss from each surface [W]. By
%                construction, Qf + Qp + Qw = 0 to machine precision:
%                each pairwise term (A_i*F_ij - A_j*F_ji) vanishes
%                exactly under reciprocity, so this identity is a valid
%                hard check for unit testing, not an approximate one.

Ebf = sigma*Tf^4;  Ebp = sigma*Tp^4;  Ebw = sigma*Tw^4;

% ---- Radiosity linear system: A*J = b, solved for J = [Jf; Jp; Jw] ----
Amat = [ 1,                    -(1-eps_f)*Ffp,     -(1-eps_f)*Ffw   ;
        -(1-eps_p)*Fpf,         1,                 -(1-eps_p)*Fpw   ;
        -(1-eps_w)*Fwf,        -(1-eps_w)*Fwp,      1-(1-eps_w)*Fww ];  % F_ww retained: wall self-view
bvec = [ eps_f*Ebf ; eps_p*Ebp ; eps_w*Ebw ];

J = Amat \ bvec;

% ---- Net radiative loss per surface: Q_i = A_i * sum_{j != i} F_ij*(J_i - J_j) ----
Qf = Af*( Ffp*(J(1)-J(2)) + Ffw*(J(1)-J(3)) );
Qp = Ap*( Fpf*(J(2)-J(1)) + Fpw*(J(2)-J(3)) );
Qw = Aw*( Fwf*(J(3)-J(1)) + Fwp*(J(3)-J(2)) );
end 


function Jpat = estimate_jacobian_pattern(Y0, P)
% Estimates the sparsity pattern of the system Jacobian by finite-
% difference probing, so that ode15s can factor a sparse Jacobian
% instead of a dense n-by-n matrix during integration.
%
% METHOD: each state component Y0(i) is perturbed individually by a
% small forward-difference step, and the right-hand side rhs(0,Y,P) is
% re-evaluated. Any output component whose value changes by more than a
% mixed absolute/relative tolerance is recorded as being coupled to
% state i. The resulting (row, column) pairs define the nonzero
% structure of the Jacobian, which is returned as a sparse 0/1 matrix.
%
% This is a one-time, temperature-independent structural probe -- it is
% called once before integration begins, not inside the ODE right-hand
% side itself, so its computational cost (n evaluations of rhs) is
% incurred only once per simulation run.
%
% INPUTS:
%   Y0 : initial state vector [n x 1], used only as the point at which
%        the perturbation is evaluated
%   P  : parameter structure passed through unchanged to rhs()
%
% OUTPUT:
%   Jpat : sparse n-by-n matrix of 0s and 1s (returned as double, per
%          the format expected by odeset's 'JPattern' option), where
%          Jpat(j,i) = 1 indicates that state i affects output j.
%          The diagonal is always set to 1, regardless of what the
%          perturbation test detects, since every state is assumed to
%          affect its own derivative and this guards against a
%          self-coupling term being missed due to numerical noise near
%          the detection threshold.

n  = numel(Y0);
f0 = rhs(0, Y0, P);

rows = [];
cols = [];

for i = 1:n
    % Perturbation step: a fixed floor of 1e-6, or 1e-6 relative to the
    % state's own magnitude, whichever is larger. This avoids a step of
    % exactly zero when Y0(i) = 0.
    dy = max(1e-6, 1e-6*abs(Y0(i)));

    Yp = Y0;
    Yp(i) = Yp(i) + dy;
    fp = rhs(0, Yp, P);

    % Mixed absolute/relative detection threshold, applied elementwise:
    % a change is considered significant if it exceeds 1e-10 in
    % absolute terms, or 1e-10 relative to the unperturbed output
    % magnitude, whichever is larger.
    changed = find(abs(fp - f0) > 1e-10*max(1, abs(f0)));

    if ~isempty(changed)
        rows = [rows; changed(:)]; %#ok<AGROW>
        cols = [cols; i*ones(numel(changed),1)]; %#ok<AGROW>
    end
end

Jpat = sparse(rows, cols, true(numel(rows),1), n, n);
Jpat = Jpat | speye(n);   % guarantee the diagonal is populated

% ---- FEED PREHEAT RECYCLE: FORCED, NOT PROBED ------------------------
% The probe above is a finite difference taken AT Y0, and Y0 is the
% t = 0 state: nothing has evaporated, the brine is essentially the
% whole feed, and the implied vapour is large enough that
% feed_preheat_T() returns a PINCHED value (and, inside the startup
% window, possibly a CAPPED one). Both branches are flat -- their
% derivative with respect to every state is exactly zero -- so the probe
% sees no coupling and silently omits it.
%
% The coupling is nonetheless real for most of the run: once the cascade
% develops, the feed temperature responds to evaporation everywhere in
% the cascade. A Jacobian pattern that omits it would force ode15s to
% take very small steps, or fail to converge, precisely when the recycle
% becomes active.
%
% The entries are therefore asserted structurally rather than measured.
% Every M, Ms, Tw state feeds the cascade evaporation total, which sets
% the feed temperature, which enters stage 1's inlet BC and hence the
% d_M, d_Ms, d_Tw derivatives of stage 1. A few false nonzeros cost
% almost nothing; a missing true nonzero costs convergence.
if isfield(P,'tfeed_dynamic') && P.tfeed_dynamic
    src = [P.idx.M(:); P.idx.Ms(:); P.idx.Tw(:)];        % everything driving evaporation
    st1 = reshape(1:P.n_field, P.Ns, P.Nx);
    st1 = st1(1,:);                                       % stage 1 nodes only
    tgt = [P.idx.M(st1).'; P.idx.Ms(st1).'; P.idx.Tw(st1).'];
    Jpat(tgt, src) = true;
end

Jpat = double(Jpat);      % odeset('JPattern', ...) expects a double sparse matrix

end 



function status = progress_report(t, ~, flag, t_total) 
% ode15s OutputFcn: prints a percent-complete progress line to the
% console during integration, throttled to at most one update every 2
% seconds of wall-clock time (not simulation time), so that solver steps
% taken in quick succession do not flood the console.
%
% INTERFACE: this function implements MATLAB's ode15s OutputFcn
% signature, (t, y, flag), and is registered via
%     odeset('OutputFcn', @(t,y,flag) progress_report(t,y,flag,P.t_sim))
% t_total is supplied through the anonymous wrapper above rather than
% being part of the standard OutputFcn signature, since ode15s itself
% only ever calls the function with (t, y, flag).
%
% ode15s calls this function with three possible values of flag:
%   'init'  : once, before integration begins
%   ''      : (empty) after each accepted step or block of steps, with
%              t possibly a vector of several time points reached since
%              the previous call
%   'done'  : once, after integration finishes
%
% The output y is required by the OutputFcn interface but is not used
% here, since only elapsed simulation time is needed to report progress
% (suppressed from the "unused input" linter warning via the #ok
% directive above).
%
% status is always returned as 0, which tells ode15s to continue
% integration; returning 1 would request early termination.
%
% INPUTS:
%   t       : current solver time, or a vector of times if several
%              steps were taken since the last call [s]
%   y       : solution state at time(s) t (unused)
%   flag    : 'init', '' (empty), or 'done' -- see above
%   t_total : total simulation duration, supplied via the anonymous
%              wrapper at the odeset() call site, used only to compute
%              the percent-complete figure [s]
%
% OUTPUT:
%   status : always 0 (continue integration)

status = 0;
persistent last_wall_clock

if strcmp(flag, 'init')
    last_wall_clock = tic;
    fprintf('  solving:      0.0%% (t = 0 / %g s)\n', t_total);
elseif isempty(flag)
    tcur = t(end);   % most recent time point, if several were reported together
    if isempty(last_wall_clock) || toc(last_wall_clock) > 2
        pct = 100*tcur/t_total;
        fprintf('  solving: %6.1f%% (t = %8.1f / %g s)\n', pct, tcur, t_total);
        last_wall_clock = tic;
    end
elseif strcmp(flag, 'done')
    fprintf('  solving: 100.0%% -- integration complete\n');
end

end % progress_report

    

function idxNN = find_nonneg_indices(P)
% Returns the indices, within the flat ODE state vector, of the state
% components that must be constrained non-negative via ode15s' built-in
% 'NonNegative' option: vapor humidity ratio, areal mass holdup, areal
% salt holdup, and film velocity.
%
% Film velocity (P.idx.u) is included because the model represents
% down-slope flow only -- the momentum equation has no mechanism for
% reverse flow, so a negative value would be non-physical rather than a
% valid transient state.
%
% INPUT:
%   P : parameter structure; must contain the index map P.idx, as built
%       elsewhere when assembling the state vector layout
%
% OUTPUT:
%   idxNN : row vector of state-vector indices to be passed to
%           odeset('NonNegative', idxNN)

idxNN = [P.idx.wv, P.idx.M, P.idx.Ms, P.idx.u];
end


function [delta, C, rho, C_raw] = invert_holdup(M, Ms, T, P)
% Recovers film thickness (delta), salt concentration (C), and density
% (rho) from the two conserved areal holdup variables -- areal mass
% holdup M = rho*delta and areal salt holdup Ms = C*delta -- together
% with the local film temperature T.
%
% METHOD: a fixed-point (Picard) iteration, needed because rho depends
% on C through water_density(T,C), while C = Ms/delta and delta = M/rho
% both depend on rho -- the three quantities are mutually coupled and
% cannot be solved for directly in closed form. water_density() runs its
% own inner Picard loop, so this is an outer iteration wrapped around an
% inner one, called at every accepted and rejected ODE step.
%
% ASSUMPTION -- two different floors on delta, for two different
% purposes:
%   - When computing C, delta is floored at P.delta_dryout_full, a
%     physically meaningful minimum "wet" film thickness, so that C does
%     not blow up unphysically as the film approaches dry-out.
%   - The delta value iterated and returned is floored only at 1e-12, a
%     numerical safeguard against division by zero, not a physical limit.
%   The returned delta is therefore allowed below the dry-out threshold
%   (a genuinely thin or dry film) while C is protected from the
%   resulting near-zero denominator.
%
% C and rho are recomputed after the loop because inside it they are
% formed from delta at the START of the iteration. On exit via max_iter
% rather than convergence, those values would correspond to the
% second-to-last delta. Recomputing guarantees the returned set exactly
% satisfies C = Ms/delta and rho = water_density(T,C) either way.
%
% LIMITATION: exhausting max_iter without meeting tol returns a stale
% result silently, with no diagnostic.
%
% INPUTS:
%   M  : areal mass holdup, rho*delta [kg/m^2]
%   Ms : areal salt holdup, C*delta   [kg/m^2]
%   T  : local film temperature        [K]
%   P  : parameter structure; must contain P.delta_dryout_full and
%        P.C_saturation
%
% OUTPUTS:
%   delta : film thickness            [m]
%   C     : salt concentration        [kg/m^3], capped at P.C_saturation
%   rho   : film density              [kg/m^3]

delta = max(M,0) / 1000;   % initial guess: M / (approx. density of fresh water)
max_iter = 40; tol = 1e-10;

for iter = 1:max_iter
    delta_prev = delta;

    C = Ms ./ max(delta, P.delta_dryout_full);   % floored to avoid blow-up near dry-out
    C = min(C, P.C_saturation);                   % physical saturation cap

    rho   = water_density(T, C);
    delta = M ./ rho;
    delta = max(delta, 1e-12);                    % numerical-only floor, not physical

    if max(abs(delta - delta_prev)./max(delta_prev,1e-12), [], 'all') < tol
        break
    end
end

% Recompute C and rho from the FINAL delta, so the returned triple is
% self-consistent even if the loop above exited via max_iter (see NOTE).
C_raw = Ms ./ max(delta, P.delta_dryout_full);   % UNCLAMPED
C     = min(C_raw, P.C_saturation);
rho   = water_density(T, C);

% C vs C_raw -- WHICH TO USE WHERE.
%   C     (clamped): property evaluation ONLY -- water_density, water_cp,
%         water_viscosity, psat_saline. The clamp is there because the
%         property correlations are not defined above saturation, so
%         feeding them C_raw would extrapolate them into nonsense.
%   C_raw (unclamped): every DIAGNOSTIC, screen, drift metric and mass
%         balance. Clamping is a property-domain guard, not physics; the
%         model has no crystallisation closure, so a run that exceeds
%         P.C_saturation is INVALID and must be visibly reported as such.
%         Reporting the clamped value instead censors the signal at the
%         ceiling: C would read exactly P.C_saturation no matter how far
%         the true Ms/delta had gone, so any convergence or drift metric
%         built on it measures how long the clamp was active rather than
%         what the salt field did. Note that Ms itself is a conserved
%         state whose PDE has no sink -- salt is NOT destroyed by the
%         clamp, only hidden from the reporting path.
end 



function Y0 = build_initial_state(P)
% Assembles the initial condition vector Y0 for the ODE integration,
% consistent with the state-vector layout defined by P.idx.
%
% Every field set here is an INITIAL CONDITION ONLY. None of these
% values are enforced as constraints during integration -- the ODE is
% free to evolve away from any of them once solving begins.
%
% METHOD, by state group:
%   Tg, Twall, Tbase, Twall_bottom : all start at ambient temperature,
%       P.Ta -- a simple cold-start assumption for every solid/glass
%       thermal mass in the system.
%   Tv  : every vapor gap starts at the fan inlet temperature, P.Tfan_in
%       -- a uniform initial guess; the true steady-state profile
%       (which develops a gradient along the cascade as vapor picks up
%       evaporate from each gap) is left for the ODE to develop.
%   wv  : seeded near equilibrium with the warm feed film rather than
%       from dry ambient air, so the evaporative driving force
%       dP = Psat_w(Tfeed) - Pv(wv0) starts close to its true
%       early-transient value instead of an artificially large one (see
%       ASSUMPTION below).
%   Tw, Tp : film starts at the feed temperature P.Tfeed (uniformly
%       across all Np plates and Nx nodes); the plate starts at ambient
%       P.Ta, reflecting that it is a separate thermal mass with no
%       reason to start warm.
%   M, Ms : areal mass and salt holdup, built from an initial film
%       thickness delta0 (see ASSUMPTION below) and the feed density and
%       concentration, via M0 = rho0.*delta0 and Ms0 = C0.*delta0 --
%       consistent with the same holdup definitions used everywhere
%       else in the model (see invert_holdup).
%   u   : seeded on the Nusselt laminar-falling-film closure (see
%       ASSUMPTION below).
%
% ASSUMPTION -- initial humidity (wv): frac_init in [0,1] sets how close
% to saturation, relative to the feed's own saturation pressure, the
% initial vapor space is taken to be; 0.7 is a reasonable default, but
% should be tuned against measured or expected enclosure humidity if
% such data is available.
%
% ASSUMPTION -- initial film thickness and velocity (delta0, u0): delta0
% is computed by compute_delta0(P) as a SINGLE SCALAR value -- the
% Nusselt-closure film thickness at the feed's own temperature, salinity,
% and areal flow rate (P.Gamma_feed) -- and is then broadcast uniformly
% across every one of the Np*Nx film nodes. This is a deliberate
% simplification: rather than attempting to guess a spatially varying
% initial profile, the entire cascade is initialized as if it were
% everywhere equal to the feed's own inlet condition, and the ODE is
% left to develop the true spatial profile during integration. This
% initial delta0 is DISTINCT from the spatially and temporally resolved
% inflow film thickness (delta_in_k) computed elsewhere during
% integration from the true, evolving inlet holdup -- delta0 is used
% only here, once, to build Y0.
%
% u0 is likewise seeded on the Nusselt closure,
%     u0 = rho0*g*sin(theta)*delta0^2 / (3*mu0),
% which is EXACTLY the same closure used for the feed's inflow velocity
% u_in_k throughout rhs() (see the matching expression at the plate
% boundary condition). This is a deliberate consistency, not a
% coincidence: seeding u at t=0 with the same closure used at the inlet
% boundary avoids forcing the momentum PDE to relax away a spurious
% cold-start velocity transient that would otherwise appear if u0 were
% chosen independently of that closure.
%
% INPUT:
%   P : parameter structure; must contain P.nstate, P.idx, P.Ta,
%       P.Tfan_in, P.Np, P.Nx, P.Tfeed, P.TDSfeed, P.g, P.theta, and
%       everything consumed by compute_delta0(), psat_saline(), and
%       water_density()/water_viscosity()
%
% OUTPUT:
%   Y0 : full initial state vector, P.nstate x 1

Y0 = zeros(P.nstate,1);

Y0(P.idx.Tg) = P.Ta;
Y0(P.idx.Tv) = P.Tfan_in * ones(P.Ns,1);

% ---- Initial vapor-gap humidity: seeded near equilibrium with the feed film ----
frac_init   = 0.7;
Psat_init   = psat_saline(P.Tfeed, P.TDSfeed, P);
Pv_init     = frac_init * Psat_init;
Patm_init   = 101325;   % matches the constant used in humidity_ratio_from_RH/vapor_partial_pressure
wv0_equil   = 0.622*Pv_init/(Patm_init - Pv_init);

Y0(P.idx.wv) = wv0_equil * ones(P.Ns,1);

Y0(P.idx.Twall)        = P.Ta * ones(P.Ns,1);   % wall segments start at ambient temperature (including Twall(Ns), the strip beside the floor)

% At t = 0 no brine has been raised and the HX has no duty, so the feed
% genuinely enters at the supplied seawater temperature. P.Tfeed is the
% cold seed here, consistent with that condition.
delta0_val = compute_delta0(P, P.Tfeed);   % scalar Nusselt-closure feed thickness, broadcast below
Tw0    = P.Tfeed   * ones(P.Ns, P.Nx);
C0     = P.TDSfeed * ones(P.Ns, P.Nx);
Tp0    = P.Ta      * ones(P.Ns, P.Nx);   % Tp0(Ns,:) is the floor's initial temperature field (starts at ambient, like every plate)
rho0   = water_density(Tw0, C0);
mu0    = water_viscosity(Tw0, C0);

M0  = rho0 .* delta0_val;   % [Ns x Nx] areal mass holdup [kg/m2]
Ms0 = C0   .* delta0_val;   % [Ns x Nx] areal salt holdup [kg/m2]

u0 = rho0 .* P.g .* sin(P.theta) .* delta0_val.^2 ./ (3.*mu0);

Y0(P.idx.M)  = M0(:);
Y0(P.idx.Tw) = Tw0(:);
Y0(P.idx.Ms) = Ms0(:);
Y0(P.idx.Tp) = Tp0(:);
Y0(P.idx.u)  = u0(:);

end 


function [mdot_local_vec, evap_per_plate_prelim] = estimate_local_mdot_vec(S, P)
% One-Picard-sweep estimate of the LOCAL total moist-air mass flow in
% each gap, needed INSIDE rhs()'s plate loop where the true mevap
% distribution isn't known yet (forward dependency: gap k's flow needs
% evaporation from plates k..Np, but the k-loop only has plates 1..k-1
% computed at iteration k). Approach: evaluate a preliminary mevap field
% with the local flow held at the constant fan-supplied value, then
% re-evaluate the local flow from that estimate. This is a single
% fixed-point sweep rather than a converged nonlinear solve, which is
% justified because sum(mevap)/mdot_da is
% O(1-10%) at most here, so a second-order iteration error is negligible
% relative to other model uncertainties, and ode15s's implicit corrector
% re-evaluates rhs() multiple times per accepted step anyway, which
% further damps any residual lag.
%
% The local duct mass flow in gap k comprises the fan-supplied dry air
% and vapour plus the evaporation accumulated from every stage below it.
% Because that evaporation depends on the mass transfer coefficient,
% which in turn depends on the local mass flow, the relation is implicit
% and is closed by Picard iteration. The converged vector is used by
% every downstream consumer -- the film transfer coefficients, the
% vapour-gap balance, the wall balance and the glazing balance -- so
% that conjugate surfaces share identical coefficients.
%
% SECOND OUTPUT. evap_per_plate_prelim is the per-stage evaporation from
% the same Picard sweep. rhs() needs a cascade-wide evaporation total to
% evaluate the feed preheat recycle at k == 1, but the k-loop has only
% filled stages 1..k-1 by then -- the identical forward dependency this
% routine already exists to resolve. Returning the estimate here rather
% than recomputing it in rhs() keeps ONE evaporation estimate in the
% model rather than two that could disagree without detection.
mdot_flat = P.mdot_da + P.mdot_vapor_in;
evap_per_plate_prelim = zeros(P.Ns,1);
mdot_local_vec = mdot_flat * ones(P.Ns,1);

% Cache the delta/C inversion: it does not change between sweeps, and it
% is the single most expensive call in this routine.
C_all = zeros(P.Ns, P.Nx);
for k = 1:P.Ns
    [~, C_all(k,:)] = invert_holdup(S.M(k,:), S.Ms(k,:), S.Tw(k,:), P);
end

n_sweeps = max(1, round(P.opt.picard_iters));
for it = 1:n_sweeps
    mevap_prelim = zeros(P.Ns, P.Nx);
    for k = 1:P.Ns
        [~, hm_prelim] = film_gap_coeffs(k, S, P, mdot_local_vec(k));
        Tw_k   = S.Tw(k,:);
        Psat_w = psat_saline(Tw_k, C_all(k,:), P);
        Pv_k   = vapor_partial_pressure(S.wv(k), S.Tv(k), P);
        dP     = evap_driving_dP(Psat_w, Pv_k, P);
        mevap_prelim(k,:) = hm_prelim .* dP .* P.Mwater ./ (P.Rgas .* Tw_k);
    end
    evap_per_plate_prelim = sum(mevap_prelim,2).*P.dx_stage*P.W;   % Ns x 1

    % Air travels Ns -> 1 (fan draws in at the bottom, gap Ns): gap k's
    % local flow picks up evaporation from stage k and every stage below.
    mdot_new = mdot_flat + flipud(cumsum(flipud(evap_per_plate_prelim)));

    if max(abs(mdot_new - mdot_local_vec)./max(mdot_local_vec,eps)) < 1e-10
        mdot_local_vec = mdot_new;
        break
    end
    mdot_local_vec = mdot_new;
end
end % estimate_local_mdot_vec


function dP = evap_driving_dP(Psat_w, Pv, P)
% Evaporation driving potential.
%
% The interfacial driving force is the difference between the saturation
% pressure at the film surface and the vapour partial pressure in the
% adjacent gap. Three treatments are provided:
%
%   'smooth' : dP = 0.5*(x + sqrt(x^2 + eps^2)) - 0.5*eps, a C-infinity
%              approximation to the positive part that agrees with
%              max(x,0) to O(eps). Restricting the model to evaporation
%              while retaining differentiability is important for the
%              implicit integrator, whose numerically-evaluated Jacobian
%              is otherwise discontinuous at any node where the film and
%              gap approach equilibrium.
%   'twoway' : the signed driving force is retained, so nodes at which
%              the gap is supersaturated with respect to the film
%              produce condensation with a corresponding latent gain to
%              the film. In a counter-current cascade the cool upper
%              films are the expected condensation sites, and this
%              closure permits the associated internal heat recovery.
%   'hard'   : max(x,0), evaluated without smoothing.
x = Psat_w - Pv;
switch lower(P.opt.evap_clamp)
    case 'smooth'
        e  = P.opt.evap_clamp_eps;
        dP = 0.5*(x + sqrt(x.^2 + e^2)) - 0.5*e;   % -0.5*e so dP(0) == 0 exactly
        dP = max(dP, 0);                            % guard against round-off sign flips
    case 'twoway'
        dP = x;
    case 'hard'
        dP = max(x, 0);
    otherwise
        error('P.opt.evap_clamp must be ''smooth'', ''twoway'' or ''hard''.');
end
end % evap_driving_dP


function [q_w, q_p, IT_next, q_refl] = solar_split(delta, IT, k, P)
% Radiative partition of the beam across one stage.
%
% The falling film attenuates the beam according to Beer-Lambert, so the
% fraction absorbed within the film is 1 - exp(-kappa_w*delta), where
% kappa_w is the spectrally-averaged absorption coefficient of the
% liquid. For film thicknesses of order 1e-4 m this fraction is small,
% and the majority of the beam reaches the solid surface beneath, where
% it is partitioned into absorption, reflection and transmission
% according to that stage's optical properties.
%
% All four streams are returned so that the partition closes exactly:
%     q_w + q_p + q_refl + IT_next == IT
% The reflected component leaves the stage and is not tracked further,
% but is returned explicitly so that it appears in the solar accounting
% rather than being absorbed into the transmitted term.
switch lower(P.opt.solar_film_model)
    case 'beer'
        a_w = 1 - exp(-P.kappa_w .* delta);          % Beer-Lambert, correct
    case 'lumped'
        a_w = P.alpha_w .* (1 - exp(-P.kappa_w .* delta));
    otherwise
        error('P.opt.solar_film_model must be ''beer'' or ''lumped''.');
end

q_w     = a_w .* IT;                                  % absorbed in the film
IT_solid = (1 - a_w) .* IT;                           % reaches the solid surface
q_p     = P.alpha_solid_stage(k) .* IT_solid;         % absorbed by plate/floor
q_refl  = P.RF_solid_stage(k)    .* IT_solid;         % reflected, leaves the stage
IT_next = (1 - P.alpha_solid_stage(k) - P.RF_solid_stage(k)) .* IT_solid;
IT_next = max(IT_next, 0);
end % solar_split

function [Tfeed_out, info] = feed_preheat_T(~, evap_per_plate, P, t)
% FEED_PREHEAT_T
% Still FEED temperature evaluated from the CURRENT cascade state,
% inside the RHS.
%
% NOTE ON THE SIGNATURE: the cascade state S is not read. The salt flow
% leaving the cascade is the conserved feed value, not a local
% terminal-node holdup ratio, so nothing in this function depends on the
% local film state; the argument is kept for a uniform call signature.
%
% The feed preheat HX is a two-stream exchanger: evaporator vapour
% condenses, incoming SEAWATER is warmed. The vapour is raised from the
% cascade brine, so the feed temperature is not an independent input --
% it is whatever the brine at this instant can support:
%
%   brine      = feed water - cascade evaporation        [kg/s]
%   mdot_vap   = brine*(1 - w_salt/w_salt_target)        [kg/s]
%   Q_avail    = Q_desuperheat + mdot_vap*h_fg           [W]
%   T_feed_out = T_feed_in + Q_avail / (mdot_feed*Cp_f)
%
% Evaluating this here rather than in an outer loop closes the recycle
% WITHIN the integration: the feed temperature tracks the cascade as it
% develops instead of being held at a single assumed value for ten
% hours. Cost is a few flops per RHS call and extra Jacobian coupling.
%
% THE RECYCLE IS NEGATIVE FEEDBACK, hence stable:
%   hotter feed -> more cascade evaporation -> less brine
%              -> less vapour -> less HX duty -> cooler feed
% It is nonetheless a genuine algebraic loop, and
% estimate_jacobian_pattern() must carry the resulting coupling between
% every stage and stage 1's inlet.
%
% CONDENSING TEMPERATURE IS P.bl.T_sat_vap, NOT P.bl.T_boil. The vapour
% is pure steam: it desuperheats by BPE degrees and then condenses
% isothermally at Tsat(P_evap). Using T_boil here would credit the
% exchanger with the boiling-point elevation, which physically never
% reaches it. See the evaporator block in build_parameters().
%
% The brine flow is taken from the instantaneous mass balance rather
% than rho*u*delta at the outlet node: the balance form is cheaper and
% smoother, which matters inside a stiff solve.
%
% The two are NOT identical. They differ by the film storage rate, and
% preheater_heater_awg() is fed the hydrodynamic outlet
% (results.brine_out.mdot) while this function uses the balance. The gap
% closes only as the film holdups reach steady state; it is reported as
% info.mdot_brine so the two can be compared directly rather than
% assumed equal.

T_sat  = P.bl.T_sat_vap;    % [K] CONDENSING temperature -- no BPE here

% ---- Brine leaving the cascade: WATER and SALT on the SAME basis ----
% Both components are taken from the global balance. Only water is
% removed by evaporation, so the salt flow leaving the cascade IS the
% salt flow entering it -- a conserved constant, not a state-dependent
% quantity. Imposing it here makes salt conservation exact by
% construction and puts this stream on the same basis as the heater in
% preheater_heater_awg(), which already closes on conserved salt.
%
% w_salt is then the MIXED-CUP outlet fraction implied by that balance,
% rather than the areal holdup ratio Ms/M at the terminal node. The node
% ratio is a local, transient quantity: pairing it with a global water
% flow mixes bases and leaves the implied salt flow unconstrained.
%
% ASSUMPTION: no net salt accumulation in the films at steady state.
% This is precisely what the terminal-stage salt-field convergence check
% tests, so a run that passes that check satisfies it.
mdot_brine_water = max(P.mdot_feed_water - sum(evap_per_plate), 0);    % [kg/s] water
mdot_salt_in     = max(P.mdot_feed_total - P.mdot_feed_water, 0);      % [kg/s] salt, conserved
mdot_brine       = mdot_brine_water + mdot_salt_in;                    % [kg/s] total
w_salt           = min(mdot_salt_in / max(mdot_brine, eps), 0.999);    % [-] mixed-cup

% ---- Vapour the electric heater raises from that brine ----
mdot_vap = max(mdot_brine * (1 - w_salt/P.bl.w_salt_target), 0);       % [kg/s]

% ---- Heat available to the air ----
% Superheat carried by the vapour is exactly BPE (T_vap_heater is the
% ELEMENT temperature; the steam leaves at the liquor temperature
% T_boil = T_sat + BPE). Desuperheating is ~0.6% of the duty but is
% carried explicitly so the accounting is closed.
Cp_sup   = 1860 + 0.12*(0.5*(P.bl.T_boil + T_sat) - 273.15);
Q_desup  = mdot_vap * Cp_sup * max(P.bl.T_boil - T_sat, 0);            % [W]
hfg_sat  = latent_heat(T_sat, 0, P);                                   % [J/kg]

% FULL condensation. With the exchanger on the FEED side the condensing
% fraction is not a free design choice: how much condenses is set by the
% feed's heat capacity rate against the available latent duty, and the
% vapour the feed cannot absorb is limited by the PINCH rather than by a
% specified fraction.
Q_avail  = Q_desup + mdot_vap * hfg_sat;                               % [W]

% ---- Feed-side heat capacity rate ----
% Evaluated at the FEED INLET temperature and salinity. Using the inlet
% keeps this term state-independent; cp of seawater varies by under 1%
% across the rise, far below the uncertainty in the vapour flow.
Cp_f     = water_cp(P.Tfeed_in, P.TDSfeed, P);                         % [J/kg-K]
Cden     = P.mdot_feed_total * Cp_f;                                   % [W/K]

% ---- CEILING AND TARGET ---------------------------------------------
% Two distinct limits, and they are NOT the same thing:
%
%   T_pinch  = T_sat - dT_pinch_HX   THERMODYNAMIC. Heat cannot flow from
%                                    condensing steam into a stream hotter
%                                    than itself. Inviolable.
%   T_target = P.Tfeed_target        DESIGN. Where we CHOOSE to run, set
%                                    below the ceiling so the film does not
%                                    enter the cascade near flash.
%
% The ceiling always wins. Without it the feed temperature can be driven
% past the condensing steam temperature, which is thermodynamically
% impossible for this exchanger, leaves the liquid range and carries the
% property correlations outside their domain. The startup guard holds the
% approach down for the first t_preheat_guard seconds, so the clamp is
% what enforces the second-law limit once the guard releases.
T_pinch  = P.bl.T_sat_vap - P.dT_pinch_HX;
T_target = min(P.Tfeed_target, T_pinch);
T_cap    = T_pinch;                          % retained: reports read this
dT_pinch = T_pinch - P.Tfeed_in;
Q_target = Cden * max(T_target - P.Tfeed_in, 0);

% ---- THREE REGIMES, RESOLVED SMOOTHLY -------------------------------
% supply-limited : Q_avail < Q_target. ALL vapour condenses (f_cond = 1)
%                  and the feed lands SHORT of target. Not an error --
%                  the feed is an inlet BC, so a cooler feed is simply a
%                  cooler cascade inlet. Reached when the cascade
%                  evaporates hard: small brine -> small vapour.
% target-limited : Q_avail > Q_target, target below the ceiling. Partial
%                  condensation; the surplus leaves for the air preheater.
% pinch-limited  : the ceiling binds before the target does.
%
% The switch is a min(), which is C0 and would cost ode15s steps every
% time the trajectory crosses it. Smoothed on the DUTY (log-sum-exp, the
% same device used for the event channel) so the RHS stays C1.
dQ = 0.02*max(Q_target, eps);                % smoothing width, ~2% of duty
Q_absorbed = Q_target - dQ*log(1 + exp((Q_target - Q_avail)/dQ));
Q_absorbed = max(min(Q_absorbed, Q_avail), 0);

Q_surplus        = max(Q_avail - Q_absorbed, 0);                % [W]
f_cond           = min(Q_absorbed / max(Q_avail, eps), 1);      % [-]
mdot_vap_cond    = f_cond * mdot_vap;                           % [kg/s]
mdot_vap_surplus = max(mdot_vap - mdot_vap_cond, 0);            % [kg/s]

dT_raw    = Q_avail / max(Cden, eps);        % unconstrained; REPORTING ONLY
Tfeed_out = P.Tfeed_in + Q_absorbed/max(Cden, eps);

% ---- STARTUP GUARD: NUMERICAL ONLY, APPLIED AFTER THE PHYSICS -------
% Softens the first minutes, when nothing has evaporated and nearly the
% whole feed reports as brine. It is a SOLVER aid and it MUST expire: a
% guard still clamping at t_end is silently setting the answer.
% With the target branch active this will rarely engage at all.
if nargin < 4 || isempty(t)
    in_startup = false;                      % no clock -> post-solve call
else
    in_startup = (t < P.t_preheat_guard);
end
capped = in_startup && (Tfeed_out - P.Tfeed_in > P.dT_preheat_cap);
if capped
    Tfeed_out  = P.Tfeed_in + P.dT_preheat_cap;
    Q_absorbed = Cden * P.dT_preheat_cap;
    Q_surplus  = max(Q_avail - Q_absorbed, 0);
    f_cond     = min(Q_absorbed / max(Q_avail, eps), 1);
    mdot_vap_cond    = f_cond * mdot_vap;
    mdot_vap_surplus = max(mdot_vap - mdot_vap_cond, 0);
end

dT            = Tfeed_out - P.Tfeed_in;
pinched       = Tfeed_out >= T_pinch - 1e-9;
dT_over_pinch = Tfeed_out - T_cap;           % <= 0 whenever the cap is active

if     Q_avail <= Q_target, regime = 'supply';
elseif pinched,             regime = 'pinch';
else,                       regime = 'target';
end

info = struct('mdot_brine',mdot_brine,'mdot_brine_water',mdot_brine_water, ...
              'mdot_salt',mdot_salt_in, ...
              'w_salt',w_salt,'mdot_vap',mdot_vap, ...
              'Q_avail',Q_avail,'Q_absorbed',Q_absorbed,'dT',dT, ...
              'pinched',pinched, ...
              'dT_raw',dT_raw,'dT_pinch',dT_pinch, ...
              'T_cap',T_cap,'Cden',Cden, ...
              'dT_over_pinch',dT_over_pinch, ...
              'capped',capped, ...
              'Q_surplus',Q_surplus, ...
              'f_cond',f_cond, ...
              'mdot_vap_cond',mdot_vap_cond, ...
              'mdot_vap_surplus',mdot_vap_surplus, ...
              'T_target',T_target, ...
              'regime',regime);
end


function dY = rhs(t, Y, P)
% Master right-hand-side: unpack -> invert holdups -> compute all
% couplings -> pack. Inlet BCs for M, Ms, Tw are imposed ONLY here,
% via the conservative flux divergence (M, Ms) and a matching upwind
% derivative using the true inflow temperature (Tw). No separate
% penalty/relaxation BC exists anywhere else in the model. The momentum
% equation for u follows the SAME philosophy: its inlet condition is
% imposed via an upwind derivative using the true inflow velocity
% u_in_k (already required by, and computed for, the M/Ms flux terms
% below) -- no separate Dirichlet-penalty BC is introduced for u.

S = unpack_state(Y, P);

[Ig, IT1] = solar_irradiance(t, P);

d_M     = zeros(P.Ns, P.Nx);
d_Tw    = zeros(P.Ns, P.Nx);
d_Ms    = zeros(P.Ns, P.Nx);
d_Tp    = zeros(P.Ns, P.Nx);
d_u     = zeros(P.Ns, P.Nx);
u_film  = zeros(P.Ns, P.Nx);
Tw_top  = zeros(P.Ns,1);
mevap   = zeros(P.Ns, P.Nx);

IT = IT1;

% ---- Enclosure radiation for every gap (1..Ns, incl. gap Ns above the
% floor) ----
% Solve the 3-surface (film / top-surface / wall) radiosity network for
% every gap ONCE here -- "single point of truth" philosophy also used
% for the inlet BCs below. This is what actually drives film<->top-
% surface radiative exchange (rather than a simple pairwise view-
% factor-1 closure), and it also supplies the wall<->film/top-surface
% radiative pathway consumed by wall_rhs().
[Qf_gap, Qtop_gap, Qwall_gap] = radiosity_all_gaps(S, P);

% ---- Local (evap-accumulating) duct mass flow ----
% Converged local duct mass flow, used by every downstream consumer
% (film transfer coefficients, vapour gap, wall, glazing) so that
% conjugate surfaces share identical coefficients.
% evap_per_plate_now is the same sweep's per-stage evaporation. It is
% consumed ONLY by the feed preheat recycle at k == 1 below, which needs
% a cascade total before the k-loop has produced one.
[mdot_local_vec, evap_per_plate_now] = estimate_local_mdot_vec(S, P);

% ---- Conjugate interface fluxes ----
% Convective exchanges are evaluated once, from the nodal fields, and
% passed as absolute powers to the vapour-gap balance. Evaluating each
% side of an interface independently would require the spatially-lumped
% form mean(hc)*(mean(Tw) - Tv)*A, which differs from the nodal integral
% sum(hc(x)*(Tw(x) - Tv))*dx*W by the covariance mean(hc'*Tw'). Since hc
% varies monotonically along the tapering gap while the film temperature
% varies monotonically along x, that covariance is systematic rather
% than random, and the two forms would not agree. Sharing a single
% nodal flux makes each conjugate pair conservative by construction.
Q_film_to_gap = zeros(P.Ns,1);   % [W] film k -> gap k, convective
Q_gap_to_top  = zeros(P.Ns,1);   % [W] gap k -> its top surface (plate k-1); entry 1 filled by vapor_gap_rhs (glass)
Q_solar_refl  = 0;               % [W] solar reflected off solid surfaces, leaves the model

% ---- Floor stage's own energy balance total (computed inside the loop
% at k==Ns, but also needed by wall_rhs() for Twall(Ns)'s gain term --
% single point of truth, same philosophy as above). ----
Q_floor_to_wallstrip = 0;   % [W] populated when k == P.Ns below

for k = 1:P.Ns
    % k = 1..Np : ordinary absorber plates.
    % k = Ns    : the floor, carrying EXACTLY the same film/solid
    %             PDE structure as every other plate, just with its own
    %             optical/thermal/geometric properties (P.*_stage(Ns))
    %             and its own back-side boundary condition (conduction to
    %             the wall-strip and to ground, instead of a gap above).

    M_k     = S.M(k,:);
    Ms_k    = S.Ms(k,:);
    Tw_k    = S.Tw(k,:);
    Tp_k    = S.Tp(k,:);
    dx_k    = P.dx_stage(k);

    % ---- Invert conserved holdups for delta, C (exact-mass approach) ----
    [delta_k, C_k, rho_w] = invert_holdup(M_k, Ms_k, Tw_k, P);
    delta_k = max(delta_k, 1e-8);

    Tv_k    = S.Tv(k);
    wv_k    = S.wv(k);

    % Film velocity u_k is a momentum-PDE state (see the momentum
    % equation below), not an algebraic Nusselt closure. mu_w is needed
    % both for the momentum equation's wall-shear term and for the
    % temperature/mass-flux terms further down.
    mu_w    = water_viscosity(Tw_k, C_k);
    u_k     = max(S.u(k,:), 0);        % clamp: down-slope flow only
    u_film(k,:) = u_k;

    % ---- Local gap flow / heat-transfer coefficients ----
    [hc_wv, hm_wv] = film_gap_coeffs(k, S, P, mdot_local_vec(k));

    % ---- Evaporation flux and vapor pressures ----
    Psat_w = psat_saline(Tw_k, C_k, P);
    Pv_k   = vapor_partial_pressure(wv_k, Tv_k, P);
    dP     = evap_driving_dP(Psat_w, Pv_k, P);
    mevap_k = hm_wv .* dP .* P.Mwater ./ (P.Rgas .* Tw_k);
    mevap(k,:) = mevap_k;

    % ---- Solar absorption terms ----
    % Same Beer-Lambert-in-the-film split for every stage, including the
    % floor (k==Ns): the water layer sitting on the floor absorbs first,
    % and whatever solar energy survives the film is then absorbed by
    % the floor itself, per P.alpha_solid_stage(Ns) = P.alpha_base.
    % Radiative partition across this stage; the reflected stream is
    % returned explicitly for the solar accounting.
    [q_solar_w, q_solar_p, IT_next, q_solar_refl_k] = solar_split(delta_k, IT, k, P);
    Q_solar_refl = Q_solar_refl + sum(q_solar_refl_k)*dx_k*P.W;
    % (IT_next is unused after k==Ns: the floor is the terminal stage, so
    % whatever isn't absorbed there is simply not tracked further, exactly
    % as the untracked reflected fraction (1-alpha_p-RF_p) already was for
    % every interior plate.)

    % ---- Plate <-> film conduction ----
    hp_w = P.k_water ./ delta_k;

    % ---- Radiative coupling: FILM k <-> top-surface/wall enclosure of gap k ----
    % Resolved via the full 3-surface (film/top/wall) radiosity network
    % (radiosity_all_gaps, solved once above) rather than a pairwise F = 1
    % closure, so the wall participates in the exchange. Qf_gap(k) = net
    % radiative loss FROM film k [W], distributed areally across the Nx
    % nodes, with the sign convention positive = film losing energy.
    q_p_w    = hp_w .* (Tp_k - Tw_k);
    q_conv_wv= hc_wv .* (Tw_k - Tv_k);
    q_evap_wv= hm_wv .* latent_heat(Tw_k,C_k,P) .* dP .* P.Mwater ./ (P.Rgas.*Tw_k);
    q_rad_k  = (Qf_gap(k)/P.Ap_stage(k)) * ones(1, P.Nx);

    % Nodal integral of the film-to-gap convective exchange. This value
    % is passed unchanged to the vapour-gap balance, so the conjugate
    % pair conserves energy to machine precision.
    Q_film_to_gap(k) = sum(q_conv_wv) * dx_k * P.W;

    % ---- Inlet (upstream) conditions for stage k, in M/Ms/Tw/u terms ----
    % THIS is the single point of truth for Variable_in(k) = Variable_out(k-1).
    % Plate 1 uses the external feed condition; every later stage (2..Ns,
    % INCLUDING the floor at k==Ns) reads the PREVIOUS stage's outlet
    % (node Nx) state directly from S -- so continuity is exact by
    % construction (see print_continuity_check), and the floor's inlet is
    % simply Plate Np's outlet, entering near the bottom as intended.
    if k == 1
        % FEED PREHEAT RECYCLE. The film does not enter at the supplied
        % seawater temperature: it enters at the feed preheat HX outlet,
        % which is driven by vapour raised from THIS cascade's own
        % brine. It is therefore evaluated from the current state, and
        % the recycle closes inside the integration -- no outer
        % iteration. See feed_preheat_T().
        if P.tfeed_dynamic
            Tw_in_k = feed_preheat_T(S, evap_per_plate_now, P, t);
        else
            Tw_in_k = P.Tfeed;
        end
        % delta0 MUST be evaluated at Tw_in_k, not at P.Tfeed. The
        % conserved inlet quantity is the areal mass flow Gamma_feed;
        % holding delta at the cold-seed value while the film enters
        % preheated inflates the inlet flux by the viscosity ratio.
        [M_in_k, Ms_in_k] = feed_inlet_state(P, Tw_in_k);
        % SALT ON A CONSERVED MASS-FRACTION BASIS, NOT ON TDS.
        % P.TDSfeed is a VOLUMETRIC concentration [kg/m3] defined at the
        % SUPPLY state. Writing Ms_in = TDSfeed*delta0 injects salt at
        % that volumetric concentration into a film that has thermally
        % expanded to rho(Tw_in) < rhow_in, so the salt mass flux comes
        % out as TDSfeed*Gamma/rho(Tw_in) instead of TDSfeed*Gamma/rhow_in
        % -- salt is MANUFACTURED at the inlet in exact proportion to the
        % preheat-driven density change, and that is the entire source
        % of the plate-1 salt residual.
        %
        % The conserved inlet quantity is the salt MASS FRACTION of the
        % supplied seawater, w_f = TDSfeed/rhow_in, which is a property
        % of the feed and does not change when the stream is heated.
        % M_in_k and Ms_in_k come back self-consistent: invert_holdup
        % will recover exactly the rho, C and delta they were built
        % from, so u_in gives rho*u*delta = Gamma_feed identically.
    else
        M_in_k     = S.M(k-1,end);
        Ms_in_k    = S.Ms(k-1,end);
        Tw_in_k    = S.Tw(k-1,end);
    end
    [delta_in_k, C_in_k, rho_in_k] = invert_holdup(M_in_k, Ms_in_k, Tw_in_k, P);
    mu_in_k  = water_viscosity(Tw_in_k, C_in_k);

    % ---- Inflow velocity: INTER-STAGE MASS/SALT CONTINUITY ----
    %
    % INVARIANT: the inlet flux of stage k is the SOLVED outlet flux of
    % stage k-1. The velocity is not re-derived from the Nusselt closure.
    %
    % The mass leaving stage k-1 is rho*u_solved*delta*W, where u_solved
    % is the MOMENTUM PDE state S.u(k-1,end). The algebraic closure
    %       u_Nu = rho*g*sin(theta)*delta_in^2 / (3*mu)
    % agrees with it only to the accuracy with which the momentum solver
    % reproduces the Nusselt profile, so evaluating the inlet velocity
    % from that closure would create or destroy mass and salt at every
    % stage hand-off. Such an imbalance is invisible per-stage -- each
    % stage still balances against its OWN declared inlet flux -- but
    % accumulates across the cascade and appears only in the cascade
    % totals, far above the machine-precision invariants the model
    % otherwise satisfies.
    %
    % Physically the film is continuous across a plate edge, so velocity
    % is continuous at the hand-off and the upstream solved velocity is
    % the correct inflow condition; using it makes the hand-off
    % conservative BY CONSTRUCTION rather than to within the closure
    % error.
    %
    % Stage 1 is the exception: there is no upstream solver state, so the
    % Nusselt closure at the feed condition is the inlet BC.
    if k == 1
        u_in_k = rho_in_k * P.g * sin(P.theta) * delta_in_k^2 / (3*mu_in_k);
    else
        u_in_k = max(S.u(k-1,end), 0);   % clamp mirrors the down-slope clamp on u_k
    end

    flux_in_M  = M_in_k  * u_in_k;    % [kg/m-s]  mass flux entering stage k
    flux_in_Ms = Ms_in_k * u_in_k;    % [kg/m-s]  salt flux entering stage k

    % ---- Advective (upwind) derivative in x for TEMPERATURE ----
    % Node 1 uses the TRUE inflow temperature Tw_in_k (not a Neumann
    % extrapolation), consistent with the conservative M/Ms inlet flux.
    dTw_dx = upwind_deriv_1field(Tw_k, dx_k, Tw_in_k);

    Cp_w   = water_cp(Tw_k, C_k, P);
    kw     = P.k_water;

    % ---- EXACT conservative mass & salt PDEs ----
    d_M(k,:)  = -ddx_flux(M_k .*u_k, dx_k, flux_in_M)  - mevap_k;
    d_Ms(k,:) = -ddx_flux(Ms_k.*u_k, dx_k, flux_in_Ms);           % no evaporation sink: salt doesn't evaporate

    % ---- Temperature PDE ----
    Tw_xx = second_diff(Tw_k, dx_k);
    d_Tw(k,:) = (1./(rho_w.*Cp_w.*delta_k)) .* ( ...
                  kw.*delta_k.*Tw_xx ...
                  + q_solar_w + q_p_w - q_conv_wv - q_evap_wv - q_rad_k ) ...
                - u_k .* dTw_dx;

    
    du_dx     = upwind_deriv_1field(u_k, dx_k, u_in_k);
    ddelta_dx = upwind_deriv_1field(delta_k, dx_k, delta_in_k);
    u_xx      = second_diff(u_k, dx_k);

    delta_k_visc = max(delta_k, P.delta_dryout_full);
    tau_b        = 3.*mu_w.*u_k ./ delta_k_visc;
    nu_w         = mu_w ./ rho_w;

    % Depth-averaged streamwise momentum. The advective term carries the
    % Boussinesq momentum correction factor beta = mean(u^2)/mean(u)^2,
    % which equals 6/5 for the parabolic velocity profile consistent
    % with the wall shear closure tau_b = 3*mu*u/delta used below.
    %
    % Two limitations of this formulation should be noted. First, the
    % wall shear closure presumes a stress-free interface and therefore
    % neglects the shear exerted by the gap airflow on the film surface.
    % Second, the streamwise viscous term nu*u_xx is retained for
    % numerical regularisation; its magnitude relative to the other
    % terms is O(delta^2/L^2) and it does not represent a physically
    % significant transport mechanism at these aspect ratios.
    d_u(k,:) = -P.opt.momentum_beta.*u_k.*du_dx ...
               + P.g*sin(P.theta) ...
               - P.g*cos(P.theta).*ddelta_dx ...
               - tau_b./(rho_w.*delta_k) ...
               + nu_w.*u_xx;

    % ---- Absorber-plate / floor-slab PDE ----
    Tp_xx = second_diff(Tp_k, dx_k);

    if k < P.Ns
        Tv_below   = S.Tv(k+1);
        % Streamwise frame reversal across the serpentine stack.
        %
        % Consecutive plates are inclined in opposite senses and carry
        % their films in opposite directions, so each plate's local
        % streamwise coordinate runs counter to its neighbour's: node i
        % of plate k lies at the opposite end of the chamber from node i
        % of plate k+1. The transfer coefficients for gap k+1 are indexed
        % in plate (k+1)'s frame, since that gap's taper is defined
        % there, and must therefore be index-reversed before being
        % applied against the plate-k temperature field.
        %
        % The inter-plate gap height varies from
        % hcomp + 2*tan(theta)*clearance to
        % hcomp + 2*tan(theta)*(chamber_L - clearance), so hc = Nu*k_a/Dh
        % varies by a comparable factor along the stage and the pairing
        % of the two fields materially affects the nodal distribution.
        [hc_wv_above, ~] = film_gap_coeffs(k+1, S, P, mdot_local_vec(k+1));
        if P.opt.alternating_plates
            hc_wv_above = fliplr(hc_wv_above);
        end
        q_conv_vp  = hc_wv_above .* (Tv_below - Tp_k);

        % Plate k forms the upper surface of gap k+1. The nodal integral
        % is taken over plate k's own area and node spacing, and the same
        % value is used by the gap balance, so the two control volumes
        % agree on both the flux and the interfacial area.
        Q_gap_to_top(k+1) = sum(q_conv_vp) * dx_k * P.W;

        % ---- Radiative gain: plate k is the TOP SURFACE of gap (k+1) ----
        % (mirrors the film/top/wall roles used for q_rad_k above),
        % resolved via that gap's own 3-surface radiosity network. Sign
        % flip: Qtop_gap(k+1) = net LOSS from plate k in that role
        % (positive = plate k radiating energy away); the GAIN added to
        % plate k's own energy balance is therefore its negative. This
        % branch covers k==Np as well: Plate Np's back side faces gap Ns,
        % the active gap above the floor, exactly as every other interior
        % plate faces the gap below it.
        q_rad_gain = -(Qtop_gap(k+1)/P.Ap_stage(k)) * ones(1, P.Nx);
    else
        % ---- Floor stage (k == Ns): back side conducts to the wrapping
        % wall-strip and to the ground, instead of facing another gap.
        % Both loss paths are evaluated from the floor's SPATIAL MEAN
        % temperature and applied as a uniform per-node areal loss,
        % consistent with the enclosure radiosity terms elsewhere in this
        % model, which are likewise spatially lumped per gap.
        Tp_bar = mean(Tp_k);
        Q_floor_to_ground    = (Tp_bar - P.Tground) / P.R_floor_ground * P.Ap_floor;              % [W]
        Q_floor_to_wallstrip = (Tp_bar - S.Twall(P.Ns)) * P.UA_floor_wall;  % [W]
        q_ambient_loss_areal = (Q_floor_to_ground + Q_floor_to_wallstrip) / P.Ap_floor;            % [W/m2]

        q_conv_vp  = zeros(1, P.Nx);
        q_rad_gain = -q_ambient_loss_areal * ones(1, P.Nx);
    end

    d_Tp(k,:) = (1./(P.rho_solid_stage(k).*P.Cp_solid_stage(k).*P.t_solid_stage(k))) .* ( ...
                  P.k_solid_stage(k).*P.t_solid_stage(k).*Tp_xx + q_solar_p - q_p_w + q_conv_vp + q_rad_gain );

    Tw_top(k) = mean(Tw_k);

    % Beam hand-off between stages of differing footprint.
    %
    % The transmitted quantity IT_next is an areal flux, but consecutive
    % stages do not share a footprint: the cascade plates span Ap = L*W
    % whereas the floor spans Ap_floor = floor_L*W. The conserved
    % quantity across the hand-off is power, so the flux is rescaled by
    % the area ratio: the power IT_next*Ap_stage(k) leaving stage k is
    % redistributed over Ap_stage(k+1) on arrival.
    if k < P.Ns
        % The beam reaches stage k+1 at the same physical location at
        % which it left stage k, which corresponds to the reversed node
        % index under the serpentine arrangement.
        if P.opt.alternating_plates
            IT_next = fliplr(IT_next);
        end
        IT = IT_next * (P.Ap_stage(k) / P.Ap_stage(k+1));
    else
        IT = IT_next;
    end

end

% ---- Vapor-space energy & moisture balances ----
[d_Tv, d_wv, Q_gap_to_top(1)] = vapor_gap_rhs(S, P, mevap, ...
                                    Q_film_to_gap, Q_gap_to_top, mdot_local_vec);

% ---- Wall-node energy balance (per gap j): vapor convects in, ambient
% conducts/convects out, and the wall exchanges radiatively with the
% film+top-surface enclosure (see radiosity_all_gaps / wall_rhs).
% Twall(Ns), the wrapping wall strip around the floor, ALSO gains the
% floor's own wall-strip conduction (Q_floor_to_wallstrip, computed
% inside the k==Ns branch above) -- passed in for energy-balance
% closure, single point of truth. ----
% Optional deposition of the aperture spill on the gap-1 wall node.
Q_spill_to_wall = 0;
if P.opt.route_solar_spill
    Q_spill_to_wall = P.opt.alpha_wall_solar * P.tau_g * Ig * P.A_spill;   % [W]
end
d_Twall = wall_rhs(S, P, Qwall_gap, mdot_local_vec, Q_floor_to_wallstrip, Q_spill_to_wall);

% ---- Glass cover energy balance ----
% The convective glazing/gap-1 exchange is supplied by vapor_gap_rhs as
% an absolute power, so both control volumes use a single value.
d_Tg = glass_rhs(S, P, Ig, Qtop_gap(1), Q_gap_to_top(1));

dY = pack_state(d_Tg, d_Tv, d_wv, d_Twall, d_M, d_Tw, d_Ms, d_Tp, d_u, P);

end % rhs
    
function S = unpack_state(Y, P)
% Splits the flat ODE state vector Y into a struct S of named,
% correctly-shaped fields using the index map built in P.idx. Every
% per-gap field (Tv, wv, Twall) spans Ns gaps, and every per-stage field
% (M, Tw, Ms, Tp, u) spans Ns stages -- stage/gap Ns is the floor, which
% carries no separate scalar node pair (see P.Ns notes in
% build_parameters).
S.Tg    = Y(P.idx.Tg);
S.Tv    = Y(P.idx.Tv);
S.wv    = Y(P.idx.wv);
S.Twall = Y(P.idx.Twall);

S.M     = reshape(Y(P.idx.M),     P.Ns, P.Nx);
S.Tw    = reshape(Y(P.idx.Tw),    P.Ns, P.Nx);
S.Ms    = reshape(Y(P.idx.Ms),    P.Ns, P.Nx);
S.Tp    = reshape(Y(P.idx.Tp),    P.Ns, P.Nx);   % S.Tp(P.Ns,:) is the floor's own temperature field
S.u     = reshape(Y(P.idx.u),     P.Ns, P.Nx);

end

function dY = pack_state(d_Tg, d_Tv, d_wv, d_Twall, d_M, d_Tw, d_Ms, d_Tp, d_u, P)
% Reassembles the individual derivative fields back into the flat ODE
% derivative vector dY, in the same P.idx layout unpack_state() reads.
dY = zeros(P.nstate,1);
dY(P.idx.Tg)    = d_Tg;
dY(P.idx.Tv)    = d_Tv;
dY(P.idx.wv)    = d_wv;
dY(P.idx.Twall) = d_Twall;
dY(P.idx.M)     = d_M(:);
dY(P.idx.Tw)    = d_Tw(:);
dY(P.idx.Ms)    = d_Ms(:);
dY(P.idx.Tp)    = d_Tp(:);
dY(P.idx.u)     = d_u(:);
end

function dTw_dx = upwind_deriv_1field(Tw_k, dx, Tw_in)
% First-order upwind (backward) spatial derivative along a single plate
% (streamwise direction), with the true inflow value Tw_in supplying the
% boundary condition at node 1 instead of a one-sided/extrapolated stencil.
dTw_dx = zeros(size(Tw_k));
dTw_dx(2:end) = (Tw_k(2:end) - Tw_k(1:end-1))/dx;
dTw_dx(1) = (Tw_k(1) - Tw_in)/dx;
end

function d = ddx_flux(flux, dx, flux_in)
% Upwind divergence of an already-computed flux field (flux = state .* u),
% used for the conservative mass/salt PDEs. flux_in is the true inflow
% flux crossing the node-1 boundary.
d = zeros(size(flux));
d(2:end) = (flux(2:end) - flux(1:end-1))/dx;
d(1) = (flux(1) - flux_in)/dx;
end

function d2 = second_diff(f, dx)
% Second spatial derivative via a centered 3-point stencil in the
% interior, and a one-sided (Neumann/zero-gradient-consistent) stencil
% at both ends of the 1-D node array f.
n = numel(f);
d2 = zeros(1,n);
d2(2:n-1) = (f(3:n) - 2*f(2:n-1) + f(1:n-2))/dx^2;
d2(1)   = (f(2)-f(1))/dx^2;
d2(n)   = (f(n-1)-f(n))/dx^2;
end

function delta0 = compute_delta0(P, T_in)
% Feed film thickness at the top of plate 1, from the Nusselt
% laminar-falling-film closure at the feed's own temperature/salinity
% and the areal feed flow rate P.Gamma_feed.
%
% T_IN IS REQUIRED WHENEVER THE FEED IS PREHEATED. The conserved
% quantity at the inlet is the areal MASS FLOW Gamma_feed, not the film
% thickness. Thickness and velocity are related through
%
%     u = rho*g*sin(theta)*delta^2 / (3*mu),   Gamma = rho*u*delta
%
% so delta must be solved at the SAME temperature the inflow velocity is
% evaluated at, or the product rho*u*delta will not equal Gamma_feed.
%
% This is material whenever the feed is preheated. Holding delta at its
% P.Tfeed (cold seed) value while the film enters at the solved preheat
% temperature inflates the inlet mass flux by exactly the viscosity
% ratio, mu(300 K)/mu(352 K) = 2.38, so plate 1 would manufacture water
% by that factor. Evaluating at the actual inlet temperature makes
% rho*u*delta = Gamma_feed identically, at any preheat level.
if nargin < 2 || isempty(T_in), T_in = P.Tfeed; end
rho0 = water_density(T_in, P.TDSfeed);
mu0  = water_viscosity(T_in, P.TDSfeed);
Gamma = P.Gamma_feed;
delta0 = (3*mu0*Gamma/(rho0^2 * P.g * sin(P.theta)))^(1/3);
end




function [M_in, Ms_in, delta0, rho0, C0] = feed_inlet_state(P, T_in)
% FEED_INLET_STATE  Stage-1 film inlet holdups, CONSISTENT with what
% invert_holdup() will recover from them.
%
% REQUIREMENT. The inlet must satisfy rho0*u(delta0)*delta0 = Gamma_feed
% exactly, where u is the Nusselt closure
%       u_in = rho*g*sin(theta)*delta^2/(3*mu).
% Building the inlet from P.TDSfeed directly does not achieve this: once
% the feed is preheated the film expands, so the same salt mass fraction
% sits at a lower volumetric concentration, and the C that invert_holdup()
% recovers is not P.TDSfeed. The rho and delta the closure was built on
% then disagree with the recovered pair. The error is second order in the
% density change but enters as delta^2, and appears as spurious mass
% creation at plate 1.
%
% METHOD. The conserved feed properties are the areal mass flow
% Gamma_feed and the salt MASS FRACTION w_f = TDSfeed/rhow_in, both
% defined at the supply state. The inlet volumetric concentration is the
% fixed point of
%       C = w_f * rho(T_in, C)
% which converges in a handful of iterations because rho depends only
% weakly on C.
w_f = P.TDSfeed / P.rhow_in;             % [-] conserved salt mass fraction

C0 = P.TDSfeed;                          % initial guess
for it = 1:50
    rho0  = water_density(T_in, C0);
    C_new = w_f * rho0;
    if abs(C_new - C0) < 1e-12*max(C0,1), C0 = C_new; break; end
    C0 = C_new;
end
rho0 = water_density(T_in, C0);
mu0  = water_viscosity(T_in, C0);

% Nusselt closure at the SELF-CONSISTENT state, so that
% rho0*u(delta0)*delta0 == Gamma_feed identically.
delta0 = (3*mu0*P.Gamma_feed/(rho0^2 * P.g * sin(P.theta)))^(1/3);

M_in  = rho0 * delta0;                   % [kg/m2] water+salt areal holdup
Ms_in = C0   * delta0;                   % [kg/m2] salt areal holdup  (= w_f*M_in)
end


function [hc_wv, hm_wv] = film_gap_coeffs(k, S, P, mdot_local_k)
% mdot_local_k [kg/s]: the ACTUAL local total moist-air mass flow through
% gap k -- i.e. mdot_da + mdot_vapor_in + (evaporation accumulated from
% every plate at or below k, in the direction of air travel, Np -> 1).
% REQUIRED argument (no fallback to a constant P.mdot_da): the duct flow
% picks up evaporated vapor at every plate, so treating it as constant
% ignores that accumulation and biases Re/Nu/hc/hm low, worst at gap 1
% where the most vapor has accumulated -- see estimate_local_mdot_vec().
%
% Compartment cross-section TAPERS along x (wedge-shaped gap between
% inclined plates, or under the tilted glass for k==1) -- see P.hgap_narrow
% / P.hgap_wide built in build_parameters(). Air velocity, Re, Dh, Nu, and
% hence hc/hm are all computed LOCALLY at each of the P.Nx streamwise
% nodes from the true local gap height, rather than from a single
% mean-area/mean-Dh value. Mass flow mdot_local_k is conserved along x;
% it is the cross-sectional area (and therefore velocity) that varies.

Tv_k  = S.Tv(k);
Tw_k = mean(S.Tw(k,:));
rho_a = air_density(Tv_k);
mu_a  = air_viscosity(Tv_k);
k_a   = air_conductivity(Tv_k);
Cp_a  = air_cp(Tv_k);
Pr_a  = mu_a * Cp_a / k_a;

% ---- local gap height, linearly interpolated feed end -> opening end ----
% NOTE: gaps 1..Np use P.L uniformly here, including gap 1, whereas
% compute_gap_view_factors takes gap 1's taper length as chamber_L. The
% two lengths differ by the end clearances, so the transfer coefficients
% and the view factors for gap 1 rest on slightly different taper spans;
% this is a disclosed inconsistency in the closure. Only gap Ns (above
% the floor's own film, footprint length floor_L =
% L+clearance -- no cos(theta) projection like the glass has) needs its
% own case.
if k == P.Ns
    L_k = P.floor_L;
else
    L_k = P.L;
end
x = linspace(0, L_k, P.Nx);                 % 1 x Nx streamwise coordinate
h_local  = P.hgap_narrow(k) + (P.hgap_wide(k) - P.hgap_narrow(k)) .* (x / L_k);   % 1 x Nx
A_local  = h_local * P.W;                                                          % 1 x Nx
Dh_local = 2*P.W*h_local ./ (P.W + h_local);                                       % 1 x Nx

% ---- local velocity, Re, Nu, hc, hm (all 1 x Nx) ----
ugap = mdot_local_k ./ (rho_a * A_local);            % 1 x Nx
Re   = rho_a .* ugap .* Dh_local / mu_a;              % 1 x Nx

% nusselt_gap has a scalar if/elseif branch on Re -- vectorize via arrayfun
Nu =  nusselt_gap(Re, Pr_a, Tv_k, Tw_k, P, k);   % 1 x Nx

hc_wv = Nu .* k_a ./ Dh_local;                        % 1 x Nx, evaluated node by node
hm_wv = hc_wv ./ (rho_a * Cp_a) * P.Le^(-2/3);        % 1 x Nx
end % film_gap_coeffs



function [d_Tv, d_wv, Q_gap1_to_glass] = vapor_gap_rhs(S, P, mevap, Q_film_to_gap, Q_gap_to_top, mdot_local_vec)
% Vapor-gap sensible-energy and moisture balances, one lumped node per gap.
%
% Convective exchanges with the film below and the plate above are
% received as absolute powers evaluated nodally by the stage loop in
% rhs(), so that each conjugate pair is conservative by construction.
%
% ENTHALPY DATUM. The latent heat h_fg does not appear explicitly in the
% three vapour streams below. On a liquid-at-Tref datum the latent
% contributions are
%     mdot_v_in*h_fg + mevap*h_fg - mdot_v_out*h_fg
% and since mdot_v_out = mdot_v_in + mevap by the mass balance enforced
% here, they cancel identically. The film correspondingly loses only the
% latent term q_evap_wv; the sensible enthalpy of the departing liquid
% is accounted for automatically by the non-conservative form of the
% film temperature equation, in which the mass-loss term multiplies the
% film temperature.
%
% This cancellation relies on all vapour generated within a gap leaving
% with the exhaust stream. If condensation onto a surface is permitted
% (P.opt.evap_clamp = 'twoway'), vapour leaves the gas phase to a
% surface, mdot_v_out ceases to equal mdot_v_in + mevap, and h_fg must
% appear explicitly in each stream. The check below identifies that
% configuration.

d_Tv = zeros(P.Ns,1);
d_wv = zeros(P.Ns,1);

evap_per_plate = sum(mevap,2).*P.dx_stage*P.W;

V_gap     = P.V_gap_vec;
Awall_vec = wall_area_per_gap(P);

% ---- Gap 1's top surface is the glass, which the k-loop cannot supply
% (it only walks stages). Computed here and returned so glass_rhs uses
% the identical value. ----
hc_gap1 = mean(film_gap_coeffs(1, S, P, mdot_local_vec(1)));
Q_gap1_to_glass = hc_gap1 * (S.Tv(1) - S.Tg) * P.Ag;
Q_gap_to_top(1) = Q_gap1_to_glass;

for j = 1:P.Ns
    if j < P.Ns
        Tin  = S.Tv(j+1);
        Min  = P.mdot_a + sum(evap_per_plate(j+1:end));
        win  = S.wv(j+1);
        Cp_in = moist_air_cp(S.Tv(j+1), S.wv(j+1));
    else
        % j == Ns: the fan delivers loop air into the bottom gap.
        %
        % CLOSED LOOP. Both inlet properties are FIXED, and neither
        % depends on the cascade state:
        %   Tin  = the electric air heater's outlet SETPOINT. The heater
        %          is specified by outlet temperature, so this is an
        %          input and its duty is the output (computed post-solve
        %          in preheater_heater_awg).
        %   win  = w_sat(T_coil). Loop air leaves the AWG saturated at
        %          the coil and the heater adds only sensible heat, so
        %          the humidity ratio is unchanged from the coil state.
        %
        % The vapour recycle acts on the FEED (see feed_preheat_T, called
        % at k == 1 in rhs). The air side carries no recycle at all, so
        % this inlet temperature is not a dynamic quantity.
        Tin  = P.T_air_in;
        Min  = P.mdot_a;
        win  = P.wv_fan_in;
        Cp_in = moist_air_cp(Tin, P.wv_fan_in);
    end

    Mout = P.mdot_a + sum(evap_per_plate(j:end));

    Tv_j  = S.Tv(j);
    wv_j  = S.wv(j);
    rho_j = moist_air_density(Tv_j, wv_j);
    Cp_j  = moist_air_cp(Tv_j, wv_j);

    Mv_j   = rho_j * V_gap(j);
    Cp_v_j = Cp_j;

    % ---- Convective exchanges: taken directly from the nodal integrals
    % computed in rhs(), rather than recomputed from lumped means ----
    if P.opt.interface_consistent
        q_conv_in = Q_film_to_gap(j);    % film j -> gap j  (positive = gap gains)
        q_conv_vp = Q_gap_to_top(j);     % gap j -> top     (positive = gap loses)
    else
        % Spatially-lumped alternative, retained for sensitivity
        % assessment. The lumped areas and mean-field evaluation do not
        % reproduce the nodal integrals used on the film and plate
        % sides, so this path is not conservative across interfaces.
        hc_local  = mean(film_gap_coeffs(j, S, P, mdot_local_vec(j)));
        q_conv_in = hc_local*(mean(S.Tw(j,:)) - S.Tv(j))*P.W*P.L_stage(j);
        if j == 1
            Tpartner_above = S.Tg;
        else
            Tpartner_above = mean(S.Tp(j-1,:));
        end
        q_conv_vp = hc_local*(S.Tv(j) - Tpartner_above)*P.W*P.L_stage(j);
    end

    hc_wall_j = mean(film_gap_coeffs(j, S, P, mdot_local_vec(j)));
    q_wall_j  = hc_wall_j * (S.Tv(j) - S.Twall(j)) * Awall_vec(j);

    % Sensible cp of the freshly evaporated vapour, evaluated at its
    % source temperature (the film surface).
    Cp_vapor_j = 1860 + 0.12*(mean(S.Tw(j,:)) - 273.15);

    d_Tv(j) = ( Min*Cp_in*(Tin-P.Tref) + q_conv_in ...
                + evap_per_plate(j)*Cp_vapor_j*(mean(S.Tw(j,:))-P.Tref) ...
                - q_conv_vp - q_wall_j - Mout*Cp_v_j*(S.Tv(j)-P.Tref) ) / (Mv_j*Cp_v_j);

    % ---- Moisture balance ----
    % The gap is treated as a well-mixed control volume of dry-air mass
    % rho_da*V, so the humidity ratio obeys
    %     V*rho_da*dw/dt = mdot_da*(w_in - w) + mdot_evap
    % whose relaxation time is V*rho_da/mdot_da. The alternative 'relax'
    % closure imposes the same steady state through a first-order lag
    % with a prescribed time constant, and is retained only to assess
    % the sensitivity of transient results to that assumption; the two
    % differ whenever the gap volume or the fan duty is varied.
    switch lower(P.opt.humidity_balance)
        case 'volume'
            rho_da_j = rho_j / (1 + wv_j);            % dry-air density in the gap
            m_da_j   = max(rho_da_j * V_gap(j), eps); % dry-air mass held in the gap
            d_wv(j)  = ( P.mdot_da*(win - wv_j) + evap_per_plate(j) ) / m_da_j;
        case 'relax'
            w_steady = win + evap_per_plate(j) / P.mdot_da;
            d_wv(j)  = (w_steady - S.wv(j)) / P.opt.tau_w_relax;
        otherwise
            error('P.opt.humidity_balance must be ''volume'' or ''relax''.');
    end
end
end % vapor_gap_rhs

function Awall = wall_area_per_gap(P)
% Side-wall area exposed to each vapor gap (gap 1, under the glass, has
% a different footprint than the repeating inter-plate gaps).
Awall = P.Awall_vec;
Awall(1) = P.Awall_first;
end % wall_area_per_gap

function d_Twall = wall_rhs(S, P, Qwall_gap, mdot_local_vec, Q_floor_to_wallstrip, Q_spill_to_wall)
% Q_floor_to_wallstrip [W]: the floor stage's own conductive loss to its
% wrapping wall strip (computed once in rhs(), single point of truth),
% positive = floor losing heat / Twall(Ns) gaining it.
% Q_spill_to_wall [W]: the portion of the admitted beam falling beyond
% the plate edge, deposited on the gap-1 wall node. Zero unless
% P.opt.route_solar_spill is enabled.
if nargin < 5, Q_floor_to_wallstrip = 0; end
if nargin < 6, Q_spill_to_wall      = 0; end

hc_vec = zeros(P.Ns,1);
for j = 1:P.Ns
    hc_j = film_gap_coeffs(j, S, P, mdot_local_vec(j));
    hc_vec(j) = mean(hc_j);
end

q_in_areal  = hc_vec .* (S.Tv - S.Twall);
q_out_areal = (S.Twall - P.Ta) ./ P.R_wall_out;

Awall_vec = wall_area_per_gap(P);

% ---- Wall radiative exchange ----
% In addition to convective exchange with the vapor (q_in_areal) and
% conduction to ambient (q_out_areal) / axially (q_cond), the wall also
% exchanges radiatively with the film+top-surface enclosure. Qwall_gap(j)
% is the net radiative LOSS from wall j to that enclosure (from
% radiosity_all_gaps, solved once in rhs()); it is converted to areal and
% subtracted here, using the same sign convention as q_out_areal (both
% are loss terms).
q_rad_areal = Qwall_gap ./ Awall_vec;

% ---- axial conduction along the continuous wall sheet ----
% The wall is discretised into one node per gap, with node j spanning a
% height L_cond(j). Conduction is therefore a property of the LINK
% between two adjacent nodes, not of either node individually: the
% thermal resistance of the link between j and j+1 is the series
% combination of the half-spans on either side,
%
%     R_link(j) = [ L_cond(j)/2 + L_cond(j+1)/2 ] / (k_wall * A_cond)
%
% and BOTH nodes must use that same value, so that the heat leaving one
% node equals the heat entering its neighbour.
%
% Assigning each node its own resistance would make the link
% non-conservative whenever adjacent spans differ. That matters here
% because the spans are strongly non-uniform: the top compartment is
% several times taller than an inter-plate gap, and the floor gap
% differs again, so the mismatch is concentrated at the two ends of the
% stack -- precisely where the floor-to-wallstrip contact injects the
% largest single heat flow in the device.
A_cond   = 2*(P.chamber_L + P.W) * P.t_wall;
L_cond   = (P.hgap_narrow + P.hgap_wide) / 2;
R_link   = (L_cond(1:end-1)/2 + L_cond(2:end)/2) ./ (P.k_wall * A_cond);   % [K/W], Ns-1 links

q_cond = zeros(P.Ns,1);   % absolute W, net INTO each node
for j = 1:(P.Ns-1)
    q_link = (S.Twall(j) - S.Twall(j+1)) / R_link(j);   % [W] j -> j+1
    q_cond(j)   = q_cond(j)   - q_link;
    q_cond(j+1) = q_cond(j+1) + q_link;
end
% ---- Floor's own wall-strip conduction: a GAIN to Twall(Ns), the wall
% wrapping the floor's edge (this contact path is unique to the floor
% stage -- no other Twall(j) receives an analogous term). ----
q_cond(P.Ns) = q_cond(P.Ns) + Q_floor_to_wallstrip;

% Aperture spill deposited on the gap-1 wall node, if enabled.
q_cond(1) = q_cond(1) + Q_spill_to_wall;

d_Twall = (q_in_areal - q_out_areal - q_rad_areal + q_cond./Awall_vec) ./ P.mCp_wall_areal;
end

function d_Tg = glass_rhs(S, P, Ig, Q_rad_glass, Q_conv_gap1_to_glass)
% Glass cover lumped energy balance.
%
% The convective exchange with gap 1 is received as an absolute power
% from vapor_gap_rhs, so the glazing and the gap agree exactly on the
% heat crossing between them.
%
% UNITS: P.mCp_g is an AREAL heat capacity [J/m^2-K] -- d_Tg is never
% multiplied by P.Ag -- so every term here must be areal [W/m^2]. The
% two incoming absolute-watt quantities are divided by P.Ag accordingly.
Tg = S.Tg;

q_rad_wg  = -Q_rad_glass         / P.Ag;   % gain = -(net radiative loss from glass)
q_conv_vg =  Q_conv_gap1_to_glass / P.Ag;  % gain from gap 1 (positive = glass heating)

hc_gamb  = P.hwind;
Tsky     = 0.0552*P.Ta^1.5;
hr_gsky  = P.eps_g*P.sigma*(Tg^2+Tsky^2)*(Tg+Tsky);

q_solar_g   = P.alpha_g*Ig;
q_conv_gamb = hc_gamb*(Tg - P.Ta);
q_rad_gsky  = hr_gsky*(Tg - Tsky);

d_Tg = (1/P.mCp_g) * (q_solar_g + q_rad_wg + q_conv_vg - q_conv_gamb - q_rad_gsky);

end % glass_rhs

function hr = radiative_coeff(T1, T2, eps1, eps2, P)
% Linearized gray-diffuse radiative heat-transfer coefficient between two
% parallel surfaces at T1, T2 with emissivities eps1, eps2, such that
% q = hr*(T1-T2). Assumes surfaces with view factor 1 (infinite parallel
% plates); used wherever the model treats a pair of surfaces as simple
% two-body radiative exchange rather than the full 3-surface enclosure
% (see solve_radiosity3 for the latter).
hr = P.sigma.*(T1.^2+T2.^2).*(T1+T2) ./ (1/eps1 + 1/eps2 - 1);
end

function [Ig, IT1] = solar_irradiance(t, P)
% Plane-of-array irradiance forcing.
%
% TWO MODES, selected by P.opt.constant_irradiance:
%
%   false (default) -- Ig(t) is interpolated from the measured PVGIS
%       record. The RHS is then NON-AUTONOMOUS: no state variable can
%       reach a true steady state, because every one tracks the diurnal
%       forcing. All steady-state diagnostics must use the quasi-steady
%       (tracking-ratio) criterion in that case.
%
%   true -- Ig = P.G_const for all t. Since every OTHER forcing in this
%       model is already a constant scalar (P.Ta, P.Vwind, P.Tground,
%       P.Tfan_in, P.Gamma_feed at fixed P.TDSfeed), fixing Ig makes the
%       system AUTONOMOUS,
%             dy/dt = f(y; P),
%       so y(t) -> y_inf and mfw(t) -> mfw_ss. The reported daily figure
%       then ceases to be an extrapolation of a windowed mean and becomes
%       an exact evaluation of a conditional design-point statement:
%       "under sustained irradiance G, the cascade produces at this rate."
%
%       This is a WEAKER claim than a diurnal yield prediction and must
%       be reported as such. Because mfw is a nonlinear functional of G,
%             mfw(<G>) ~= <mfw(G(t))>,
%       so the constant-G daily total is NOT the daily total obtained
%       from the measured record. Quantify that gap once (see
%       P.G_const guidance in build_parameters) rather than assuming it
%       is negligible.
if P.opt.constant_irradiance
    Ig  = P.G_const;
else
    t_clock = mod(P.t_clock_offset + t, 86400);
    Ig = interp1(P.t_irr_data, P.Gi_irr_data, t_clock, 'pchip', 0);
    Ig = max(Ig, 0);
end
IT1 = P.tau_g * Ig;
end

function w = humidity_ratio_from_RH(RH, T, ~)
% Converts relative humidity RH (0-1) at temperature T [K] to humidity
% ratio w [kg water vapor / kg dry air], at fixed atmospheric pressure.
psat = psat_pure(T);
Pv   = RH*psat;
Patm = 101325;
w = 0.622*Pv/(Patm-Pv);
end

function Pv = vapor_partial_pressure(w, ~, ~)
% Inverse of humidity_ratio_from_RH: recovers water-vapor partial
% pressure Pv [Pa] from humidity ratio w at fixed atmospheric pressure.
Patm = 101325;
Pv = w*Patm/(0.622+w);
end

function psat = psat_pure(T)
% Saturation vapor pressure of pure water [Pa] as a function of
% temperature T [K] (empirical correlation).
psat = exp(23.1964 - 3816.44./(T-46.13));
end

function psat = psat_saline(T, C, ~)
% Saturation vapor pressure of SALINE water [Pa]: pure-water psat_pure(T)
% reduced by a salinity correction factor, given areal salt
% concentration C [kg/m^3] converted internally to salinity S [g/kg].
rho = water_density(T, C);
S   = 1000.*C./rho;
psat = psat_pure(T) .* (1 - 0.000537.*S - 1.0278e-6.*S.^2);
end

function rho = water_density(T, C)
% Density of saline water [kg/m^3] as a function of temperature T [K]
% and salt concentration C [kg/m^3]. Iterates because the correlation is
% naturally expressed in salinity S = 1000*C/rho (kg salt per m^3
% solution vs. g salt per kg solution), which itself depends on rho.
t = T - 273.15;
rho0 = 999.9 + 2.034e-2*t - 6.162e-3*t.^2 + ...
       2.261e-5*t.^3 - 4.657e-8*t.^4;

rho = rho0;

for k = 1:50
    S = C./rho;

    rho_new = rho0 + ...
        (802*S - 2.001*S.*t + 1.677e-2*S.*t.^2 ...
        -3.060e-5*S.*t.^3 -1.613e-5*S.^2.*t.^2);

    if abs(rho_new-rho) < 1e-3
        break;
    end

    rho = rho_new;
end
end

function mu = water_viscosity(T, C)
% Dynamic viscosity of saline water [Pa.s] as a function of temperature
% T [K] and salt concentration C [kg/m^3].
rho = water_density(T, C);
S_gkg = 1000.*C./rho;
S = S_gkg/1000;                     % convert g/kg -> kg/kg mass-fraction scale
t = T - 273.15;
mu_w = 4.2844e-5 + 1./(0.157*(t+64.993).^2 - 91.296);
A = 1.541 + 1.998e-2*t - 9.52e-5*t.^2;
B = 7.974 - 7.561e-2*t + 4.724e-4*t.^2;
mu = mu_w .* (1 + A.*S + B.*S.^2);
end

function Cp = water_cp(T, C,~)
% Specific heat of saline water [J/kg-K] as a function of temperature T
% [K] and salt concentration C [kg/m^3].
rho = water_density(T, C);
S = 1000.*C./rho;         % g/kg
A = 5.328 - 9.76e-2*S + 4.04e-4*S.^2;
B = -6.913e-3 + 7.351e-4*S - 3.15e-6*S.^2;
Cc = 9.6e-6 - 1.927e-6*S + 8.23e-9*S.^2;
D = -2.5e-9 + 1.666e-9*S - 7.125e-12*S.^2;
Cp = 1000*(A + B.*T + Cc.*T.^2 + D.*T.^3);   % kJ/kg·K -> J/kg·K to match your existing units
end

function Lv = latent_heat(T, C, ~)
% Pure-water latent heat (linear fit, valid over this model's ~30-45 C
% operating range) times the Sharqawy et al. (2010) salinity
% correction: hfg,sw = hfg,w * (1 - S/1000), S in g/kg. Without this
% factor q_evap/Q_evap are biased high, by an amount that grows with
% feed salinity.
% C is optional (areal salt concentration, kg/m^3); if omitted, falls
% back to pure water (S=0), so call sites without C remain valid.
Tc = T - 273.15;
Lv_w = 2.501e6 - 2361.*Tc;
if nargin < 2 || isempty(C)
    Lv = Lv_w;
    return
end
if nargin < 3
    rho = 1000;   % coarse fallback if P not supplied; avoids a hard error
else
    rho = water_density(T, C);
end
S = 1000.*C./rho;             % g/kg
Lv = Lv_w .* (1 - S./1000);
end

function rho = air_density(T)
% Ideal-gas density of dry air [kg/m^3] at atmospheric pressure and
% temperature T [K].
rho = 101325 ./ (287.05 .* T);
end

function Nu = nusselt_gap(Re, Pr, T_hot_side, T_cold_side, P, k)
% Nusselt correlation for the fan-driven duct flow inside a vapour gap.
%
% The turbulent branch uses the Dittus-Boelter correlation with a fixed
% Prandtl exponent of 0.4. Distinguishing the heating and cooling
% exponents (0.4 and 0.3) alters Nu by approximately 4% at Pr = 0.7,
% which is well within the scatter of the correlation itself, while
% introducing a discontinuity in dNu/dT at every node where the gap and
% surface temperatures cross. The surface-temperature arguments are
% retained in the signature for the 'switched' variant.
%
% The laminar-turbulent transition is represented by a smoothstep blend
% over 2300 <= Re <= 4000 rather than a discontinuous switch, since the
% jump between the fully-developed laminar value and the Dittus-Boelter
% value at Re = 2300 is approximately a factor of two and would
% otherwise appear as a discontinuity in the numerically-evaluated
% Jacobian.
%
% LAMINAR CONSTANT -- BOUNDARY CONDITION. Every gap here exchanges with
% TWO thermally active surfaces: the film below and the plate (or the
% glazing, for gap 1) above. For fully developed laminar flow between
% parallel plates the limiting Nusselt number is not a single number but
% a function of the wall-flux ratio
%       r = q_top / q_bottom        (both positive INTO the air)
% through the exact superposition result
%       Nu = 5.385 / (1 - 0.3461*r)
% which reproduces the standard anchors: r = 0 -> 5.385 (one surface
% adiabatic), r = +1 -> 8.235 (equal fluxes), r = -1 -> 4.000 (heated
% below, cooled above; also obtainable in closed form, since the profile
% is then linear and Tb = (T1+T2)/2).
%
% The single constant 5.385 quoted for ducts is the
% ONE-SURFACE-ADIABATIC value, i.e. the r = 0 special case, and does not
% describe this geometry.
%
% r is a MODEL OUTPUT, not an input. Evaluated on a converged solution
% (see diagnose_flux_ratio.m) it separates cleanly into two regimes:
%
%   gaps 3..Ns : r = 0.90 to 1.09, tightly bounded. Since r ~ 1 implies
%                T_plate_above ~ T_film_below with the air cooler than
%                both, this IS the symmetric two-surface problem, not an
%                approximation to it. Nu = 7.541 is used: the both-walls-
%                ISOTHERMAL value rather than 8.235 (both-walls uniform
%                FLUX), because the plates are metal with resolved axial
%                conduction and the model carries explicit surface
%                temperature fields, so the isothermal family is the apt
%                one. The 9% gap between 7.541 and 8.235 is the residual
%                boundary-condition uncertainty.
%
%   gaps 1..2  : r = -0.10 and -0.58 respectively -- heated from the film
%                below, cooled through the top (the glazing radiates to
%                sky; plate 1's underside is the coolest solid in the
%                stack). Nu = 5.21 and 4.49 from the relation above.
%
% ASSUMPTION AND ITS IMPLICATION: r is treated as a fixed per-gap
% constant, evaluated at end-of-run. r is driven by the solar temperature
% ladder and is therefore time-varying; holding it fixed is a quasi-steady
% approximation to the thermal boundary condition, not to the flow. It
% requires re-evaluation for a different stack geometry, plate count, or
% air path. Nu_lam is a first-order parameter for the evaporative flux
% (hm scales with it one-for-one through Chilton-Colburn) and belongs in
% the sensitivity screening with range [4.0, 8.4].
%
% REMAINING LIMITATIONS. (i) Dittus-Boelter is established for Re > 10000
% and L/D > 10, so its use from Re = 2300 upward is an extrapolation;
% Gnielinski would cover the transition band properly. (ii) The values
% above are FULLY DEVELOPED asymptotes. If the thermal entry length
% L_th ~ 0.05*Re*Pr*Dh is comparable to the gap length, the mean Nusselt
% number exceeds these constants and an entry correlation with them as
% asymptotes is required. Neither is addressed here.
if nargin < 5 || ~isfield(P,'opt') || ~isfield(P.opt,'nusselt_model')
    model = 'smooth';
else
    model = P.opt.nusselt_model;
end

% Per-gap laminar constant. k is optional: call sites that omit the gap
% index fall back to the value for the symmetric interior gaps, which are
% the large majority.
if nargin >= 6 && ~isempty(k) && isfield(P,'Nu_lam_gap') ...
        && k >= 1 && k <= numel(P.Nu_lam_gap)
    Nu_lam = P.Nu_lam_gap(k);
elseif isfield(P,'Nu_lam_default')
    Nu_lam = P.Nu_lam_default;
else
    Nu_lam = 7.541;
end

switch lower(model)
    case 'switched'
        Nu = zeros(size(Re));
        laminar = Re < 2300;
        Nu(laminar) = Nu_lam;
        if T_hot_side > T_cold_side
            Nu(~laminar) = 0.023*Re(~laminar).^0.8.*Pr.^0.3;
        else
            Nu(~laminar) = 0.023*Re(~laminar).^0.8.*Pr.^0.4;
        end

    case 'smooth'
        Re_lo = 2300;    % start of transition band
        Re_hi = 4000;    % end of transition band
        Re_s  = max(Re, eps);
        Nu_turb = 0.023*Re_s.^0.8.*Pr.^0.4;

        % Smoothstep blend: C^1 continuous at both band edges.
        xi = min(max((Re_s - Re_lo)/(Re_hi - Re_lo), 0), 1);
        w  = xi.^2 .* (3 - 2*xi);
        Nu = (1 - w).*Nu_lam + w.*Nu_turb;

    otherwise
        error('P.opt.nusselt_model must be ''smooth'' or ''switched''.');
end
end

function mu = air_viscosity(T)
% Dynamic viscosity of dry air [Pa.s] at temperature T [K] (Sutherland's law).
Tref = 273.15;
Sc   = 110.4;
mu   = 1.716e-5 .* (T./Tref).^1.5 .* (Tref + Sc) ./ (T + Sc);
end

function k = air_conductivity(T)
% Thermal conductivity of dry air [W/m-K] at temperature T [K] (Sutherland's law).
Tref = 273.15;
Sc   = 194;
k    = 0.0241 .* (T./Tref).^1.5 .* (Tref + Sc) ./ (T + Sc);
end

function Cp = air_cp(T)
% Specific heat of dry air [J/kg-K] at temperature T [K] (linear fit).
Cp = 1005 + 0.05 .* (T - 273.15);
end

function rho_m = moist_air_density(T, w)
% Density of moist air [kg/m^3] at temperature T [K] and humidity ratio
% w [kg water/kg dry air].
rho_da_T = 101325 ./ (287.05 .* T);
rho_m    = rho_da_T .* (1 + w) ./ (1 + 1.608 .* w);
end

function Cp_m = moist_air_cp(T, w)
% Mass-weighted specific heat of moist air [J/kg-K] (dry air + water
% vapor) at temperature T [K] and humidity ratio w.
Cp_da_T  = 1005 + 0.05 .* (T - 273.15);
Cp_vap_T = 1860 + 0.12 .* (T - 273.15);
Cp_m     = (Cp_da_T + w .* Cp_vap_T) ./ (1 + w);
end

function [value, isterminal, direction] = zldd_events(~, Y, P) 
S = unpack_state(Y, P);
delta_all = invert_holdup(S.M(:), S.Ms(:), S.Tw(:), P);
% One event channel per node (NOT min(delta_all)-P.delta_dryout).
% Collapsing all nodes to a single scalar via min() makes the effective
% event function non-smooth right at the crossing -- whichever node is
% "the minimum" can flip between solver steps, producing a kink exactly
% where odezero is trying to root-find. That kink is what causes the
% "odezero: an event disappeared (internal error)" failure mode. Giving each
% node its own smooth event value (and letting ode15s/odezero handle the
% multi-event bookkeeping) avoids the kink; isterminal=1 for all of them
% still stops the integration as soon as the FIRST node crosses.
value      = delta_all(:) - P.delta_dryout;
isterminal = ones(size(value));
direction  = -ones(size(value));
end



function results = extract_results(t, Y, P)
% Post-processes the raw ode15s trajectory (t, Y) into a fully expanded,
% human- and plot-friendly results struct: recovers delta/C at every
% saved time step, rebuilds per-plate inlet/outlet tables, integrates
% cumulative energy and distillate quantities, evaluates efficiency and
% GOR metrics, and runs the final-time mass/energy balance checks.
% This is the single place all downstream reporting/plotting functions
% pull their data from -- nothing below here touches Y or t directly.


nT = size(Y,1);
x  = (0.5:P.Nx-0.5)*P.dx_stage(1);   % streamwise node centers, plate-1 spacing (reporting/plotting reference axis; stage Ns uses its own dx_stage(P.Ns) internally wherever it matters)

results.t = t;
results.x = x;
% ---- Floor's own x-axis: stage Ns's node spacing is dx_stage(Ns), built
% from floor_L (= L+clearance), NOT plate-1's dx_stage(1) used for
% results.x above. Anything plotting a per-node profile for stage Ns
% against a streamwise coordinate needs this instead. ----
results.x_base = (0.5:P.Nx-0.5)*P.dx_stage(P.Ns);
results.P = P;

results.Tg = Y(:,P.idx.Tg);
results.Tv = Y(:,P.idx.Tv);
results.wv = Y(:,P.idx.wv);
results.Twall = Y(:,P.idx.Twall);
% NOTE: the floor carries no separate scalar temperature fields -- its
% temperature is results.Tp(P.Ns,:,:) (spatially resolved) and its wall
% segment is results.Twall(:,P.Ns), exactly like every other stage/gap.
% See build_parameters (P.Ns) and rhs().

results.Tw    = zeros(P.Ns, P.Nx, nT);
results.delta = zeros(P.Ns, P.Nx, nT);
results.C     = zeros(P.Ns, P.Nx, nT);
results.Ms    = zeros(P.Ns, P.Nx, nT);   % conserved salt holdup state [kg/m2]
results.C_raw = zeros(P.Ns, P.Nx, nT);   % unclamped Ms/delta [kg/m3 solution]
results.Tp    = zeros(P.Ns, P.Nx, nT);
results.u     = zeros(P.Ns, P.Nx, nT);
results.mevap = zeros(P.Ns, P.Nx, nT);
results.mevap_plate = zeros(nT, P.Ns);
results.I     = zeros(nT,1);

results.Tw_in_true    = zeros(nT, P.Ns);
results.delta_in_true = zeros(nT, P.Ns);
results.C_in_true     = zeros(nT, P.Ns);
results.rho_in_true   = zeros(nT, P.Ns);
results.u_in_true     = zeros(nT, P.Ns);
results.hfg_plate     = zeros(nT, P.Ns);
for i = 1:nT
    S = unpack_state(Y(i,:).', P);
    [S.delta, S.C, ~, S.C_raw] = invert_holdup(S.M, S.Ms, S.Tw, P);

    results.Tw(:,:,i)    = S.Tw;
    results.delta(:,:,i) = S.delta;
    results.C(:,:,i)     = S.C;
    % Ms is a SOLVER STATE and is stored verbatim. Downstream salt
    % accounting must read this, never C.*delta: C is clamped at
    % P.C_saturation, so reconstructing Ms from it manufactures a salt
    % deficit wherever the clamp is active. C_raw is the uncensored
    % concentration, for diagnostics and screens only.
    results.Ms(:,:,i)    = S.Ms;
    results.C_raw(:,:,i) = S.C_raw;
    results.Tp(:,:,i)    = S.Tp;

    [Ig, ~] = solar_irradiance(t(i), P);
    results.I(i) = Ig;

    % ---- Pass 1: preliminary evaporation using the uniform initial
    % duct-flow approximation, purely to build the local mass-flow
    % profile below. This mirrors estimate_local_mdot_vec() in rhs() --
    % same Picard-sweep justification (single sweep, not a converged
    % nonlinear solve; residual lag is second-order in sum(mevap)/mdot_da). ----
    mdot_flat = P.mdot_da + P.mdot_vapor_in;
    mevap_prelim_i = zeros(P.Ns, P.Nx);
    for ki_p = 1:P.Ns
        [~, hm_prelim] = film_gap_coeffs(ki_p, S, P, mdot_flat);
        Psat_w_p = psat_saline(S.Tw(ki_p,:), S.C(ki_p,:), P);
        Pv_k_p   = vapor_partial_pressure(S.wv(ki_p), S.Tv(ki_p), P);
        dP_p     = evap_driving_dP(Psat_w_p, Pv_k_p, P);
        mevap_prelim_i(ki_p,:) = hm_prelim .* dP_p .* P.Mwater ./ (P.Rgas .* S.Tw(ki_p,:));
    end
    evap_per_plate_prelim_i = sum(mevap_prelim_i,2).*P.dx_stage*P.W;   % Ns x 1
    mdot_local_vec_i = mdot_flat + flipud(cumsum(flipud(evap_per_plate_prelim_i)));

    % ---- Pass 2: final quantities using the converged local mass flow ----
for ki = 1:P.Ns
 
        Tw_k    = S.Tw(ki,:);
        C_k     = S.C(ki,:);
        Tv_k    = S.Tv(ki);
        wv_k    = S.wv(ki);

        % u is a genuine state -- read it directly rather than
        % recomputing it from the Nusselt algebraic closure.
        results.u(ki,:,i) = max(S.u(ki,:), 0);

        [~, hm_wv] = film_gap_coeffs(ki, S, P, mdot_local_vec_i(ki));
        Psat_w = psat_saline(Tw_k, C_k, P);
        Pv_k   = vapor_partial_pressure(wv_k, Tv_k, P);
        dP     = evap_driving_dP(Psat_w, Pv_k, P);
        mevap_k = hm_wv .* dP .* P.Mwater ./ (P.Rgas .* Tw_k);
        results.mevap(ki,:,i) = mevap_k;
        results.mevap_plate(i,ki) = sum(mevap_k) * P.dx_stage(ki) * P.W;
        results.hfg_plate(i,ki) = mean(latent_heat(S.Tw(ki,:), C_k, P));
        if ki == 1
            % Reproduce the RHS inlet BC exactly, INCLUDING the
            % time-varying preheat. Using P.Tfeed here would use the
            % cold seed (this loop runs before the post-solve feed
            % reconciliation) and would not describe the stream the
            % cascade actually received at this instant.
            if P.tfeed_dynamic
                Tw_in_k = feed_preheat_T(S, results.mevap_plate(i,:).', P, t(i));
            else
                Tw_in_k = P.Tfeed;
            end
            % Self-consistent inlet, IDENTICAL to the rhs construction.
            % Salt must enter on the conserved areal-holdup basis: the
            % volumetric form Ms = TDSfeed*delta0 is not invariant under
            % invert_holdup, so it would not reproduce the solver's own
            % inlet state.
            [M_in_k, Ms_in_k] = feed_inlet_state(P, Tw_in_k);
        else
            M_in_k  = S.M(ki-1,end);
            Ms_in_k = S.Ms(ki-1,end);
            Tw_in_k = S.Tw(ki-1,end);
        end
        [delta_in_k, C_in_k, rho_in_k] = invert_holdup(M_in_k, Ms_in_k, Tw_in_k, P);
        mu_in_k = water_viscosity(Tw_in_k, C_in_k);
        % Must mirror rhs() EXACTLY, or the reported inlet flux will not
        % be the one the solver actually imposed. Stage 1 uses the
        % Nusselt closure at the feed condition; every downstream stage
        % inherits the upstream solved velocity (see the conservation
        % note in rhs()).
        if ki == 1
            u_in_k = rho_in_k * P.g * sin(P.theta) * delta_in_k^2 / (3*mu_in_k);
        else
            u_in_k = max(S.u(ki-1,end), 0);
        end

        results.Tw_in_true(i,ki)    = Tw_in_k;
        results.delta_in_true(i,ki) = delta_in_k;
        results.C_in_true(i,ki)     = C_in_k;
        results.rho_in_true(i,ki)   = rho_in_k;
        results.u_in_true(i,ki)     = u_in_k;
end
end

results.mfw = sum(results.mevap_plate, 2);

% NOTE: reshape (not squeeze) is used here so the Ns dimension is
% preserved even when P.Ns == 1 (i.e. Np == 0, a degenerate case).
% squeeze() strips ALL singleton dimensions, so for Ns==1 it would
% collapse (1,1,nT) straight to (nT,1) instead of (1,nT), and the
% subsequent transpose would then produce a (1,nT) matrix instead of the
% intended (nT,1) -- silently breaking every downstream
% results.*_in/out(iE,k) index.
results.Tw_in     = reshape(results.Tw(:,1,:),     P.Ns, []).';
results.Tw_out    = reshape(results.Tw(:,end,:),   P.Ns, []).';
results.delta_in  = reshape(results.delta(:,1,:),  P.Ns, []).';
results.delta_out = reshape(results.delta(:,end,:),P.Ns, []).';
results.u_in      = reshape(results.u(:,1,:),      P.Ns, []).';
results.u_out     = reshape(results.u(:,end,:),    P.Ns, []).';
results.C_in      = reshape(results.C(:,1,:),      P.Ns, []).';
results.C_out     = reshape(results.C(:,end,:),    P.Ns, []).';
% Uncensored companion to C_out. C_out is clamped at P.C_saturation, so
% any convergence, drift or trajectory metric built on it saturates at
% the ceiling and reports how long the clamp was active rather than what
% the salt field did. Use C_out for property-consistent reporting and
% C_out_raw for every diagnostic that asks WHERE THE SALT WENT.
results.C_out_raw = reshape(results.C_raw(:,end,:), P.Ns, []).';

results.Tp_in     = reshape(results.Tp(:,1,:),   P.Ns, []).';
results.Tp_out    = reshape(results.Tp(:,end,:), P.Ns, []).';

results.Tv_out    = results.Tv;
results.Tv_in     = zeros(nT, P.Ns);
results.Tv_in(:,end)     =P.Tfan_in;    % fan draws in at gap Ns, the gap above the floor
if P.Ns >= 2
    results.Tv_in(:,1:end-1) = results.Tv(:,2:end);
end

results.I_layer          = zeros(nT, P.Ns);
results.q_solar_w_plate  = zeros(nT, P.Ns);
results.q_solar_p_plate  = zeros(nT, P.Ns);
results.q_solar_refl_plate = zeros(nT, P.Ns);   % beam reflected from solid surfaces
results.q_solar_terminal   = zeros(nT, 1);      % beam power leaving the terminal stage
for i = 1:nT
    [~, IT1_i] = solar_irradiance(t(i), P);
    IT_node = IT1_i * ones(1, P.Nx);
    for kj = 1:P.Ns
        delta_k = max(results.delta(kj,:,i), 1e-6);
        % Uses the same partition routine as rhs(), so this diagnostic
        % cannot diverge from the absorption applied by the solver.
        [q_w_node, q_p_node, IT_next_node, q_refl_node] = solar_split(delta_k, IT_node, kj, P);
        results.I_layer(i,kj)         = mean(IT_node);
        results.q_solar_w_plate(i,kj) = sum(q_w_node)    * P.dx_stage(kj) * P.W;
        results.q_solar_p_plate(i,kj) = sum(q_p_node)    * P.dx_stage(kj) * P.W;
        results.q_solar_refl_plate(i,kj) = sum(q_refl_node) * P.dx_stage(kj) * P.W;
        % area-scaled hand-off, as in rhs()
        if kj < P.Ns
            if P.opt.alternating_plates
                IT_next_node = fliplr(IT_next_node);
            end
            IT_node = IT_next_node * (P.Ap_stage(kj) / P.Ap_stage(kj+1));
        else
            IT_node = IT_next_node;
            results.q_solar_terminal(i) = sum(IT_node) * P.dx_stage(kj) * P.W;
        end
    end
end

results.mdot_air_in         = P.mdot_da * ones(nT,1);
results.mdot_humid_per_gap  = P.mdot_da + P.mdot_vapor_in + ...
    fliplr(cumsum(fliplr(results.mevap_plate), 2));
results.mdot_air_out        = results.mdot_humid_per_gap(:,1);

u_final = results.u(:,:,end);
u_final = max(u_final, 1e-9);
results.retention_time_plate = sum(P.dx_stage ./ u_final, 2);
results.retention_time_total = sum(results.retention_time_plate);
results.retention_time_cumulative = cumsum(results.retention_time_plate);

% ---- RECONCILE THE CASCADE INLET TEMPERATURE FIRST ----------------
% compute_mass_balance rebuilds the stage-1 inlet from scratch and must
% do it at the SAME temperature the RHS used, or its plate-1 row is a
% different stream from the one that was solved. The full reconciliation
% block further down runs too late for that, so the hot-side value is
% established here, before any consumer reads it.
if P.tfeed_dynamic
    S_end_pre = unpack_state(Y(end,:).', P);
    [Tfeed_pre, ~] = feed_preheat_T(S_end_pre, results.mevap_plate(end,:).', P, t(end));
    P.Tfeed_cascade_in         = Tfeed_pre;
    results.P.Tfeed_cascade_in = Tfeed_pre;
end

results.mass_balance = compute_mass_balance(results);

Nx = P.Nx;

hfg_ref = latent_heat(mean(results.Tw(:)), mean(results.C(:)), P);

results.total_evap_final   = sum(results.mevap_plate(end,:));

% ---- Cumulative distillate (needed below; also stored for later use) ----
results.cum_distillate = cumtrapz(t, results.mfw);

% ---- Initial-spike relaxation cutoff ----
% t=0 initial conditions (vapor space not yet in equilibrium with the
% warm feed film) cause a fast, non-physical evaporation spike before
% the vapor gap catches up. Detect it directly from mfw(t) rather than
% assuming a fixed window length. Search is restricted to the first
% t_relax_search_max seconds so this doesn't accidentally grab a later
% (e.g. solar-noon) peak/trough if irradiance is time-varying.
%
% IMPORTANT: the spike peak is NOT the right cutoff by itself. After the
% initial spike, mfw typically keeps DECLINING for a while (the vapor
% gap is still relaxing toward the film/plate temperatures) before it
% turns around and climbs smoothly under the real solar/film dynamics.
% That declining tail after the spike is still transient relaxation,
% not real production, so it must also be excluded. We therefore:
%   1) find the spike peak (as before), then
%   2) find the first point AFTER that peak where mfw stops decreasing
%      and starts increasing again (a slope-sign-change trough), and
%      use that time as the actual relaxation cutoff.
% NOTE: islocalmin() is NOT used here. islocalmin
% auto-derives a prominence threshold from the overall variation of
% the whole signal, so a genuine but SHALLOW post-spike dip (small
% relative to the full-day rise in mfw) gets
% filtered out as "noise" and islocalmin returns no candidate at all
% -- even though the dip is clearly visible and real in mfw(t). A
% direct sign-change-in-slope test has no such global-scale blind
% spot: it flags the first point where the trend locally reverses,
% regardless of how deep that reversal is.
% If no reversal is found before the search window ends (e.g. a very
% short run, or a spike-free run), fall back to the peak itself so
% behavior degrades gracefully rather than erroring.
% IMPORTANT: searching for the LARGEST value within a fixed window is
% NOT robust. If the long-term smooth production rise is steep enough,
% mfw(t) can climb past the original spike's height before the search
% window closes -- so "max within [0, t_relax_search_max]" ends up
% picking a point on the smooth rising branch near the window's edge,
% not the actual spike. The spike and trough are therefore detected by
% their SHAPE (a rise-then-fall, followed by a fall-then-rise), not by
% their magnitude relative to the rest of the run:
%   1) find the FIRST point where mfw stops rising and starts falling
%      (the spike peak, however tall or short it is), then
%   2) find the FIRST point after that where mfw stops falling and
%      starts rising again (the trough), and use that as t_relax.
%
% NOTE on robustness: reversals are confirmed by TOLERANCE (the value
% must move by more than tol_abs from the running extremum), not by
% requiring the new slope to persist for a fixed number of samples.
% A sample-count persistence check silently fails whenever the actual
% decline/rise only lasts 1-2 samples before flattening into a slow
% plateau (as commonly happens right after a fast startup spike) --
% the true reversal gets rejected simply because it didn't "hold" for
% N more samples, even though it clearly happened. Tracking the
% running max/min directly and requiring only a small tolerance move
% handles both a 1-sample-long reversal and a 50-sample-long one
% identically, and still ignores genuine step-to-step solver noise
% (tolerance is set relative to the local signal range).
t_relax_search_max = min(0.5*t(end), 3600);
n_search = find(t <= t_relax_search_max, 1, 'last');
mfw_win = results.mfw(1:n_search);
local_range = max(mfw_win) - min(mfw_win);
if local_range <= 0
    local_range = max(abs(mfw_win)) + eps;
end
tol_abs = 1e-3 * local_range;   % ignore moves smaller than 0.1% of the local signal range

% ---- Step 1: track running max; confirm peak once value drops by > tol_abs from it ----
idx_peak = 1;
has_spike = false;
for jj = 2:n_search
    if mfw_win(jj) > mfw_win(idx_peak)
        idx_peak = jj;              % still rising -- update running peak candidate
    elseif mfw_win(idx_peak) - mfw_win(jj) > tol_abs
        has_spike = true;           % genuine decline confirmed -- peak locked in
        break
    end
end
peak_rate = mfw_win(idx_peak);

% ---- Step 2: from the peak, track running min; confirm trough once value rises by > tol_abs from it ----

trough_rate = NaN;


if has_spike
    idx_trough = idx_peak;
    for jj = (idx_peak+1):n_search
        if mfw_win(jj) < mfw_win(idx_trough)
            idx_trough = jj;                 % still falling -- update running trough candidate
        elseif mfw_win(jj) - mfw_win(idx_trough) > tol_abs
            break                            % genuine rise confirmed -- trough locked in
        end
    end
    idx_relax   = idx_trough;
    trough_rate = mfw_win(idx_trough);
else
    idx_relax = 1;
end



t_relax = t(idx_relax);

results.spike_peak_rate    = peak_rate;
results.spike_peak_time    = t(idx_peak);
results.spike_trough_rate  = trough_rate;
results.spike_relax_time   = t_relax;   % trough time, i.e. the end of the startup relaxation

if t(end) - t_relax < 0.1*t(end)
    warning(['Less than 10%% of the run remains after the detected ' ...
             'post-spike relaxation time (t_relax = %.1f s, ' ...
             't_end = %.1f s). Extend t_sim before trusting the ' ...
             'reported production rate.'], t_relax, t(end));
end

% ---- Clean accumulated production over the post-relaxation window (diagnostic total) ----
win = t > t_relax;
results.window = win;
distillate_at_relax = interp1(t, results.cum_distillate, t_relax);
total_distillate_clean = results.cum_distillate(end) - distillate_at_relax;
duration_clean          = t(end) - t_relax;

% ---- Clean average production rate (post-relaxation window only) ----
mfw_avg_rate = trapz(t(win), results.mfw(win)) / duration_clean;

results.total_distillate_clean = total_distillate_clean;   % [kg], true integral, post-relaxation only (diagnostic)
results.duration_clean          = duration_clean;            % [s]
results.average_production = mfw_avg_rate;                    % [kg/s], post-relaxation average rate
results.mfw_hourly = mfw_avg_rate * 3600;

% ---- Daily production ------------------------------------------------
%
% The estimator BRANCHES on the irradiance forcing mode, because the two
% modes support fundamentally different claims.
%
% ================= CONSTANT-IRRADIANCE MODE (autonomous) =================
% With P.opt.constant_irradiance = true every forcing is constant, the
% RHS is autonomous, and the trajectory relaxes to a fixed point. The
% production rate is then genuinely CONSTANT, so
%       mfw_daily = mfw_ss * t_operating
% is not an extrapolation at all: it is the exact evaluation of a
% conditional design-point statement ("sustained at G, the cascade
% produces at this rate"). mfw_ss is taken as the mean over the FINAL
% stationary window, not over [t_relax, t_end] -- the latter still spans
% the whole settling transient and would bias the mean low.
%
% Stationarity is asymptotic, not automatic, so the relative drift of
% mfw across that window is computed and STORED
% (results.mfw_ss_drift_rel) and a warning is raised if it exceeds
% P.tol_stationary. Do not quote mfw_daily from a run that warns: it
% means t_sim was too short for the slowest mode. That mode is the film
% SALT field in the terminal (floor) stage, not any thermal node -- the
% floor receives the concentrate of every stage above it and has the
% longest settling time in the system.
%
% ================= DIURNAL MODE (non-autonomous, measured record) ========
% No state is ever stationary, so no steady rate exists to project. When
% the run already spans the full operating day, the only defensible
% figure is the TRUE INTEGRAL of production over the clean window.
% Multiplying a clean-window MEAN RATE by the full t_operating would
% re-attribute production onto the very interval the spike-exclusion
% logic deliberately discards and would overstate the result by the
% factor t_operating/duration_clean, which grows with the length of the
% discarded transient; that bias would propagate into
% water_recovery_pct, Edistill, and every efficiency and SEC figure.
% Rate extrapolation is retained ONLY for genuinely short runs, where it
% is a projection onto unsimulated time rather than a re-attribution.
sim_duration_hr = t(end) / 3600;

if P.opt.constant_irradiance
    % ---- Stationary window: final fraction of the run ----
    t_stat_start = t(end) - P.frac_stationary_window * (t(end) - t(1));
    win_stat     = t >= t_stat_start;
    if nnz(win_stat) < 3
        win_stat = true(size(t));   % degenerate short run: fall back to all points
    end
    mfw_ss = trapz(t(win_stat), results.mfw(win_stat)) / ...
             max(t(end) - min(t(win_stat)), eps);

    % Relative drift of mfw across the stationary window: least-squares
    % slope over the window, expressed as the fractional change the slope
    % would produce over the window itself. This is the quantitative
    % convergence criterion that should be quoted alongside the yield.
    tw_s   = t(win_stat) - min(t(win_stat));
    pmfw   = polyfit(tw_s, results.mfw(win_stat), 1);
    mfw_ss_drift_rel = abs(pmfw(1)) * (max(tw_s)) / max(abs(mfw_ss), eps);

    results.mfw_ss           = mfw_ss;              % [kg/s] steady-state production rate
    results.mfw_ss_drift_rel = mfw_ss_drift_rel;    % [-]    fractional drift over the window
    results.mfw_ss_window    = win_stat;

    if mfw_ss_drift_rel > P.tol_stationary
        warning('ZLDD:NotStationary', ...
            ['Distillate rate has NOT converged: relative drift over the ' ...
             'final %.0f%% of the run is %.2e (tolerance %.1e). mfw_daily ' ...
             'and every derived efficiency are unreliable. Increase ' ...
             'P.t_sim -- the limiting mode is the terminal-stage salt ' ...
             'field, not the thermal nodes.'], ...
             100*P.frac_stationary_window, mfw_ss_drift_rel, P.tol_stationary);
    end

    results.average_production = mfw_ss;            % override: steady rate, not transient-contaminated mean
    results.mfw_hourly         = mfw_ss * 3600;
    results.mfw_daily          = mfw_ss * P.t_operating * 3600;   % [kg] -- exact, not extrapolated
    results.mfw_daily_basis    = 'steady-state rate x t_operating (constant irradiance)';

else
    if sim_duration_hr >= P.t_operating
        % Full operating day simulated: use the true integral. Production
        % over [0, t_relax] is treated as exactly zero, which is the
        % assumption the spike-exclusion logic already makes.
        results.mfw_daily       = total_distillate_clean;
        results.mfw_daily_basis = 'true integral over clean window (diurnal)';
    else
        results.mfw_daily       = results.mfw_hourly * P.t_operating;
        results.mfw_daily_basis = 'rate extrapolation, run shorter than operating day (diurnal)';
    end
    results.mfw_ss           = NaN;   % no steady state exists under time-varying forcing
    results.mfw_ss_drift_rel = NaN;
end

results.water_recovery_pct =100* results.mfw_daily/(P.Vfeed*P.rhow_in/1000-P.Vfeed*P.TDSfeed/1000);

% ---- Esolar: must share the SAME temporal basis as mfw_daily ----------
% Otherwise every efficiency and SEC figure is biased by the ratio
% t_sim/t_operating. In constant-irradiance mode t_sim is deliberately
% DECOUPLED from t_operating: t_sim is a numerical convergence parameter
% (run it long enough to reach the fixed point), whereas t_operating is a
% reporting basis. They are distinct quantities and must not be
% conflated.
if P.opt.constant_irradiance
    % Analytic and unambiguous -- no quadrature needed.
    Esolar = P.G_const * P.Ag * P.t_operating * 3600;              % [J]
elseif sim_duration_hr >= 0.99*P.t_operating
    Esolar = trapz(t, results.I) * P.Ag;
else
    % Shorter-than-operating-day run: extrapolate consistently with
    % mfw_daily above, using the same average-rate-times-duration logic.
    Esolar = trapz(t, results.I) * P.Ag / t(end) * P.t_operating*3600;
end


%Edistill = trapz(t, sum(results.mevap_plate .* results.hfg_plate, 2)); 
Edistill = results.mfw_daily * hfg_ref;   % [J], energy of the daily distillate mass on whatever basis mfw_daily was formed (see results.mfw_daily_basis) -- Esolar above is formed on the SAME basis, so the ratio is dimensionally and temporally consistent
results.eta_first_effect = trapz(t, results.mevap_plate(:,1))*hfg_ref / max(Esolar,eps);
results.hfg_ref          = hfg_ref;
results.Esolar_incident  = Esolar;

Edistill_per_plate = trapz(t, results.mevap_plate) * hfg_ref;
Esolar_per_plate   = trapz(t, results.I_layer .* P.Ap);
results.Edistill_per_plate = Edistill_per_plate;
results.Esolar_per_plate   = Esolar_per_plate;
results.eta_plate_share    = Edistill_per_plate ./ max(Esolar, eps);
results.eta_plate_local    = Edistill_per_plate ./ max(Esolar_per_plate, eps);

results.E_solar_incident       = Esolar;
results.E_solar_glass          = trapz(t, P.alpha_g .* results.I .* P.Ag);
results.E_solar_films_total    = trapz(t, sum(results.q_solar_w_plate, 2));
results.E_solar_plates_total   = trapz(t, sum(results.q_solar_p_plate, 2));
results.E_solar_reflected_total = trapz(t, sum(results.q_solar_refl_plate, 2));
results.E_solar_absorbed_total = results.E_solar_glass ...
                               + results.E_solar_films_total ...
                               + results.E_solar_plates_total;

% ---- Aperture accounting ----
% E_solar_admitted is what actually passes the glass; E_solar_targeted is
% what the model applies to plate 1. The difference is the spill.
results.E_solar_admitted = trapz(t, P.tau_g .* results.I) * P.Ag;
results.E_solar_targeted = trapz(t, P.tau_g .* results.I) * P.Ap;
results.E_solar_spill    = results.E_solar_admitted - results.E_solar_targeted;
% Whatever survives the terminal stage (floor) is not tracked further.
% Evaluated directly at the terminal stage rather than inferred by
% difference, so the closure residual below is an independent check.
results.E_solar_terminal_transmitted = trapz(t, results.q_solar_terminal);
results.solar_closure_residual = results.E_solar_targeted ...
                               - results.E_solar_films_total ...
                               - results.E_solar_plates_total ...
                               - results.E_solar_reflected_total ...
                               - results.E_solar_terminal_transmitted;
results.E_distillate_total     = Edistill;

Tsky = 0.0552*P.Ta^1.5;
results.Tsky = Tsky * ones(nT,1);
results.q_glass_conv    = P.hwind * (results.Tg - P.Ta) * P.Ag;
hr_gsky_t = P.eps_g*P.sigma*(results.Tg.^2 + Tsky^2).*(results.Tg + Tsky);
results.q_glass_rad_sky = hr_gsky_t .* (results.Tg - Tsky) * P.Ag;
results.E_loss_glass_conv = trapz(t, results.q_glass_conv);
results.E_loss_glass_rad  = trapz(t, results.q_glass_rad_sky);

rho_brine_f = water_density(results.Tw(:,end,end), results.C(:,end,end));
Cp_brine_f  = water_cp(results.Tw(:,end,end),      results.C(:,end,end), P);
mdot_brine  = rho_brine_f(end) * results.u(end,end,end) * results.delta(end,end,end) * P.W;
results.E_loss_brine_sensible = mdot_brine * mean(Cp_brine_f) * ...
                                (results.Tw(end,end,end) - P.Tfeed) * t(end);

Ns = P.Ns;

results.q_evap_plate    = zeros(nT, Ns);
results.q_pw_plate      = zeros(nT, Ns);
results.q_conv_wv_plate = zeros(nT, Ns);
results.q_rad_wp_plate  = zeros(nT, Ns);
results.q_conv_vp_plate = zeros(nT, Ns);
results.q_rad_gain_plate= zeros(nT, Ns);
results.q_glass_rad_wg  = zeros(nT, 1);
results.q_rad_wall_plate= zeros(nT, Ns);   % wall's net radiative loss per gap [W]

results.hc_wv = zeros(nT, Ns);
results.hm_wv = zeros(nT, Ns);
results.hp_w  = zeros(nT, Ns);
results.hr_wp = zeros(nT, Ns);
results.hc_vg = zeros(nT, 1);

results.Re_gap = zeros(nT, Ns);
results.Nu_gap = zeros(nT, Ns);
results.Pr_gap = zeros(nT, Ns);
results.Sh_gap = zeros(nT, Ns);
results.Re_film= zeros(nT, Ns);

results.Pv_gap   = zeros(nT, Ns);
results.RH_gap   = zeros(nT, Ns);
results.Tdew_gap = zeros(nT, Ns);

results.Psat_w     = zeros(Ns, Nx, nT);
results.dP_driving = zeros(Ns, Nx, nT);
results.rho_w  = zeros(Ns, Nx, nT);
results.mu_w   = zeros(Ns, Nx, nT);
results.Cp_w   = zeros(Ns, Nx, nT);
results.Lv     = zeros(Ns, Nx, nT);
results.Sgkg   = zeros(Ns, Nx, nT);

for i = 1:nT
    S.Tg    = results.Tg(i);
    S.Tv    = results.Tv(i,:).';
    S.wv    = results.wv(i,:).';
    S.Tw    = results.Tw(:,:,i);
    S.delta = results.delta(:,:,i);
    S.C     = results.C(:,:,i);
    S.Tp    = results.Tp(:,:,i);

    % Solve the SAME 3-surface radiosity network used inside rhs() --
    % single point of truth, so this diagnostic cannot silently drift
    % from what the solver actually computed.
    [Qf_gap_i, Qtop_gap_i, Qwall_gap_i] = radiosity_all_gaps(S, P);

    results.q_glass_rad_wg(i) = -Qtop_gap_i(1);   % gain to glass = -(loss from glass)

    % ---- Exact local mass-flow vector for THIS saved time step ----
    % results.mevap_plate(i,:) is already fully populated by the earlier
    % pass above, so no lag/estimate is needed here (unlike inside rhs()'s
    % forward k-loop).
    mdot_local_vec_i = P.mdot_da + P.mdot_vapor_in + ...
        flipud(cumsum(flipud(results.mevap_plate(i,:).')));

    Tv1 = S.Tv(1);

    rho_a1 = air_density(Tv1); mu_a1 = air_viscosity(Tv1);
    k_a1   = air_conductivity(Tv1); Cp_a1 = air_cp(Tv1);
    Pr_a1  = mu_a1 * Cp_a1 / k_a1;

    % Local (tapered) gap-1 geometry, NOT the flat P.Acomp/P.Dh -- gap 1 is
    % under the glass and uses P.hgap_narrow(1)/P.hgap_wide(1), which differ
    % from the inter-plate compartments (see build_parameters). Kept as an
    % inline calc rather than calling film_gap_coeffs(1,...) because the
    % cold-side comparison here is S.Tg (glass), not mean(S.Tw(1,:)) as
    % film_gap_coeffs assumes -- gap 1's "cold side" is physically the
    % glass, not a plate film.
    x1 = linspace(0, P.chamber_L, P.Nx);
    h1_local  = P.hgap_narrow(1) + (P.hgap_wide(1) - P.hgap_narrow(1)) .* (x1 / P.chamber_L);
    A1_local  = h1_local * P.W;
    Dh1_local = 2*P.W*h1_local ./ (P.W + h1_local);

    ugap1 = mdot_local_vec_i(1) ./ (rho_a1 * A1_local);
    Re1   = rho_a1 .* ugap1 .* Dh1_local / mu_a1;

    Nu1 = nusselt_gap(Re1, Pr_a1, Tv1, S.Tg, P, 1);

    results.hc_vg(i) = mean(Nu1 .* k_a1 ./ Dh1_local);

    for k = 1:Ns
        delta_k = max(S.delta(k,:), 1e-6);
        Tw_k    = S.Tw(k,:); 
        C_k     = S.C(k,:); 
        Tp_k    = S.Tp(k,:);
        Tv_k    = S.Tv(k);   
        wv_k    = S.wv(k);
        dx_k    = P.dx_stage(k);

        rho_w = water_density(Tw_k, C_k);
        mu_w  = water_viscosity(Tw_k, C_k);
        Cp_w  = water_cp(Tw_k, C_k, P);
        Lv    = latent_heat(Tw_k, C_k, P);
        Sgkg  = 1000 .* C_k ./ rho_w;

        results.rho_w(k,:,i) = rho_w;
        results.mu_w(k,:,i)  = mu_w;
        results.Cp_w(k,:,i)  = Cp_w;
        results.Lv(k,:,i)    = Lv;
        results.Sgkg(k,:,i)  = Sgkg;

        Psat_w = psat_saline(Tw_k, C_k, P);
        Pv_k   = vapor_partial_pressure(wv_k, Tv_k, P);
        dP     = evap_driving_dP(Psat_w, Pv_k, P);
        results.Psat_w(k,:,i)     = Psat_w;
        results.dP_driving(k,:,i) = dP;

        [hc_wv, hm_wv] = film_gap_coeffs(k, S, P, mdot_local_vec_i(k));
        results.hc_wv(i,k) = mean(hc_wv);
        results.hm_wv(i,k) = mean(hm_wv);
        results.hp_w(i,k)  = mean(P.k_water ./ delta_k);

        % results.hr_wp is kept as an INFORMATIONAL, pairwise-equivalent
        % (view factor = 1) reference coefficient only -- the actual
        % energy balance is driven by Qf_gap_i(k)/Qtop_gap_i, from the
        % 3-surface network solved once above. Guarded divide avoids
        % blow-up when the two temperatures are nearly equal.
        if k == 1
            Tpartner_above = S.Tg * ones(1,Nx); eps_above = P.eps_g;
        else
            Tpartner_above = S.Tp(k-1,:);       eps_above = P.eps_p;   % k==Ns: partner is Plate Np, an ordinary plate -- eps_p is correct here too
        end
        hr_wp = radiative_coeff(Tw_k, Tpartner_above, P.eps_w, eps_above, P);
        results.hr_wp(i,k) = mean(hr_wp);

        results.q_evap_plate(i,k)    = sum(hm_wv .* Lv .* dP .* P.Mwater ./ (P.Rgas.*Tw_k)) * dx_k * P.W;
        results.q_pw_plate(i,k)      = sum((P.k_water./delta_k) .* (Tp_k - Tw_k)) * dx_k * P.W;
        results.q_conv_wv_plate(i,k) = sum(hc_wv .* (Tw_k - Tv_k)) * dx_k * P.W;

        % ---- Actual radiative terms used by the solver (from the shared
        % radiosity_all_gaps solve above). The pairwise hr_wp/hr_gain
        % closures computed just above feed only the informational
        % results.hr_wp field, not these terms. ----
        results.q_rad_wp_plate(i,k)   = Qf_gap_i(k);          % net loss from film k
        results.q_rad_wall_plate(i,k) = Qwall_gap_i(k);       % net loss from wall k

         if k < Ns
            Tv_below = S.Tv(k+1);
            [hc_wv_above, ~] = film_gap_coeffs(k+1, S, P, mdot_local_vec_i(k+1));
            if P.opt.alternating_plates
                hc_wv_above = fliplr(hc_wv_above);   % serpentine frame reversal
            end
            results.q_conv_vp_plate(i,k)  = sum(hc_wv_above .* (Tv_below - Tp_k)) * dx_k * P.W;
            % Plate k is the TOP SURFACE of gap (k+1); gain = -(loss from that role)
            results.q_rad_gain_plate(i,k) = -Qtop_gap_i(k+1);
        else
            % Floor stage (k == Ns): mirrors rhs()'s k==Ns branch exactly
            % (ground + wall-strip conduction, evaluated from the floor's
            % spatial mean temperature), so this diagnostic matches what
            % the solver actually used.
            Tp_bar = mean(Tp_k);
            Q_floor_to_ground_i    = (Tp_bar - P.Tground) / P.R_floor_ground * P.Ap_floor;
            Q_floor_to_wallstrip_i = (Tp_bar - S.Twall(P.Ns)) * P.UA_floor_wall;
            results.q_conv_vp_plate(i,k)  = 0;
            results.q_rad_gain_plate(i,k) = -(Q_floor_to_ground_i + Q_floor_to_wallstrip_i);
        end

        rho_a = air_density(Tv_k); mu_a = air_viscosity(Tv_k);
        k_a   = air_conductivity(Tv_k); Cp_a = air_cp(Tv_k);
        Pr    = mu_a * Cp_a / k_a;

        % Local (tapered) gap-k geometry, NOT the flat P.Acomp/P.Dh -- see
        % P.hgap_narrow(k)/P.hgap_wide(k) in build_parameters. Re_gap/Nu_gap
        % are reported here as the MEAN across the Nx streamwise nodes,
        % since velocity/Re/Nu vary along x with the taper;
        % a single scalar is a summary, not the true local value at any x.
        % k==Ns (the floor's own gap) uses floor_L = L+clearance (its own
        % footprint length, no cos(theta) projection like the glass has).
        if k == 1
            L_gap_k = P.chamber_L;
        elseif k == Ns
            L_gap_k = P.floor_L;
        else
            L_gap_k = P.L;
        end
        xk = linspace(0, L_gap_k, P.Nx);
        hk_local  = P.hgap_narrow(k) + (P.hgap_wide(k) - P.hgap_narrow(k)) .* (xk / L_gap_k);
        Ak_local  = hk_local * P.W;
        Dhk_local = 2*P.W*hk_local ./ (P.W + hk_local);

        ugap = P.mdot_da ./ (rho_a * Ak_local);
        Re_x = rho_a .* ugap .* Dhk_local / mu_a;
        Nu_x = nusselt_gap(Re_x, Pr, Tv_k, mean(Tw_k), P, k);

        results.Re_gap(i,k) = mean(Re_x);
        results.Nu_gap(i,k) = mean(Nu_x);
        results.Pr_gap(i,k) = Pr;
        results.Sh_gap(i,k) = mean(Nu_x) * (P.Le)^(1/3);
        Gamma = mean(rho_w) * mean(results.u(k,:,i)) * mean(delta_k);
        results.Re_film(i,k) = 4*Gamma / mean(mu_w);

        results.Pv_gap(i,k)   = Pv_k;
        results.RH_gap(i,k)   = Pv_k / max(psat_pure(Tv_k), eps);
        results.Tdew_gap(i,k) = 46.13 + 3816.44 / max(23.1964 - log(max(Pv_k,1e-6)), 1e-6);
    end
end
results.rho_in    = reshape(results.rho_w(:,1,:),   P.Ns, []).';
results.rho_out   = reshape(results.rho_w(:,end,:), P.Ns, []).';

results.Tw_final    = results.Tw(:,:,end);
results.Tp_final    = results.Tp(:,:,end);
results.delta_final = results.delta(:,:,end);
results.C_final     = results.C(:,:,end);
results.u_final     = results.u(:,:,end);
results.mevap_final = results.mevap(:,:,end);
results.Tv_final    = results.Tv(end,:).';
results.wv_final    = results.wv(end,:).';
results.Tg_final    = results.Tg(end);
% NOTE: "Tbase" is the floor's OWN Tp field (stage Ns), which is
% spatially resolved rather than a single lumped scalar. The scalar field
% name carries its spatial mean per time step, for the reporting code.
results.Tbase       = squeeze(mean(results.Tp(P.Ns,:,:), 2));   % [nT x 1]
results.Tbase_final = results.Tbase(end);

% ============ 1. Q_wall_total (all gaps, quasi-steady lumped through Twall) ============
Awall_vec = wall_area_per_gap(P);   % single source of truth: per-k tapered wall area, not a flat per-gap value

q_wall_areal = (results.Twall - P.Ta) ./ P.R_wall_out;      % [nT x Np], W/m^2
q_wall_total_t = q_wall_areal * Awall_vec;                  % [nT x 1], W (sum over gaps)

Q_wall_ss   = mean(q_wall_total_t(win));                    % window-mean rate, W
Q_wall_end  = q_wall_total_t(end);                          % instantaneous rate at t_final, W
results.Q_wall  = Q_wall_ss * P.t_operating * 3600;

% ============ 2. Q_feed_out  (concentrated brine leaving the floor stage, Ns) ============
mdot_brine_t = results.rho_out(:,end) .* results.u_out(:,end) .* ...
               results.delta_out(:,end) .* P.W;             % [nT x 1], kg/s
Cp_brine_t   = squeeze(results.Cp_w(end,end,:));            % already stored, [nT x 1]

H_brine_out_t = mdot_brine_t .* Cp_brine_t .* (results.Tw_out(:,end) - P.Ta);

Q_feed_out_ss  = mean(H_brine_out_t(win));                  % window-mean rate, W
Q_feed_out_end = H_brine_out_t(end);                        % instantaneous rate at t_final, W
results.Q_feed_out = Q_feed_out_ss * P.t_operating * 3600;  

% ---- Brine stream leaving the cascade, as a stream specification ----
% This is the feed to the crystalliser. It IS the outlet of
% the last plate: the same stream the PART D mass balance measures at
% x = L of stage Ns. There is no unit between them.
%
% BASIS. Reported at the FINAL INSTANT, not on the stationary window.
% Window-mean values are also formed below under _win names, for
% reference only, and must not be fed to the heater: mdot and w_salt are
% averaged INDEPENDENTLY over the window, and w_salt is a ratio of two
% quantities that are still drifting, so mean(C)/mean(rho) is not the
% salt fraction of any stream that ever existed. A heater inlet
% back-computed from such a pair as salt/w_salt matches the cascade
% outlet in neither flow rate nor salt fraction, so water and salt would
% be lost at the cascade -> heater handoff -- an inconsistency between two
% streams rather than a difference of averaging basis. Taking the final
% instant makes the S1w = S3 + S4 identity close to the storage term
% (verified and printed in PART B), and makes the heater sizing
% consistent with the audited outlet.
%
% The final instant is a legitimate operating point here: the
% quasi-steady check at the end of PART D confirms the floor stage and
% plate 1 both track the forcing without lag over the closing window,
% so there is no residual transient to average away.
iE_b = numel(results.t);

results.brine_out.mdot     = mdot_brine_t(iE_b);                         % [kg/s] total brine
results.brine_out.C        = results.C_out(iE_b,end);                    % [kg/m3] TDS
results.brine_out.rho      = results.rho_out(iE_b,end);                  % [kg/m3]
results.brine_out.T        = results.Tw_out(iE_b,end);                   % [K]
results.brine_out.w_salt   = results.brine_out.C / results.brine_out.rho;% [-] salt mass fraction
results.brine_out.CF       = results.brine_out.C / P.TDSfeed;            % [-] concentration factor
results.brine_out.mdot_salt = results.brine_out.mdot * results.brine_out.w_salt;   % [kg/s]
results.brine_out.mdot_water= results.brine_out.mdot - results.brine_out.mdot_salt;% [kg/s]
results.brine_out.basis     = 'final instant (t_end), plate Ns outlet face';

% Window-mean values, RETAINED FOR REFERENCE ONLY. Not fed to any unit.
results.brine_out.mdot_win   = mean(mdot_brine_t(win));                  % [kg/s]
results.brine_out.C_win      = mean(results.C_out(win,end));             % [kg/m3]
% UNCENSORED companions. C and C_win above are clamped at P.C_saturation
% and are the values consistent with the properties actually evaluated.
% The two below are the true Ms/delta, and are what the convergence and
% drift diagnostics must use -- a clamped signal cannot show how far the
% terminal salt field ran past the ceiling.
if isfield(results,'C_out_raw')
    results.brine_out.C_raw     = results.C_out_raw(iE_b,end);           % [kg/m3]
    results.brine_out.C_win_raw = mean(results.C_out_raw(win,end));      % [kg/m3]
    results.brine_out.CF_raw    = results.brine_out.C_raw / P.TDSfeed;   % [-]
end
results.brine_out.rho_win    = mean(results.rho_out(win,end));           % [kg/m3]
results.brine_out.T_win      = mean(results.Tw_out(win,end));            % [K]
results.brine_out.CF_win     = results.brine_out.C_win / P.TDSfeed;      % [-]
results.brine_out.mdot_end   = mdot_brine_t(iE_b);   % alias field name, same value

% Cross-check against the independently computed mass balance. Both
% measure the same stream; a mismatch means the outlet-face
% extrapolation and the node-centre field have diverged, which is a
% discretisation problem, not a reporting one.
if isfield(results,'mass_balance') && isfield(results.mass_balance,'total')
    % BOTH SIDES TOTAL BRINE. results.brine_out.mdot is the total stream
    % (water + salt); the PART D terminal must therefore be brine_total,
    % not brine_water. Comparing against brine_water would read the outlet
    % salt mass fraction out as a false mismatch of the same magnitude as
    % that fraction.
    mb_ref = results.mass_balance.total.brine_total;
    d_rel = abs(results.brine_out.mdot - mb_ref) / max(mb_ref, eps);
    results.brine_out.check_vs_MB = d_rel;
    if d_rel > 1e-3
        warning('brine_out:MBmismatch', ...
           ['Cascade outlet differs from the PART D mass balance by %.2e ' ...
            '(rel). Stream table and audit will not agree.'], d_rel);
    end
end

% ============ 3. Q_vapor_out (humid air exiting gap 1, fan-drawn) ============
mdot_evap_t    = sum(results.mevap_plate, 2);                  % [nT x 1], kg/s
mdot_vap_out_t = P.mdot_vapor_in + mdot_evap_t;                % [nT x 1], kg/s
Cp_vapor_out_t = 1860 + 0.12*(results.Tv(:,1) - 273.15);       % air_cp-vapor form, inlined
Cp_da_out_t    = 1005 + 0.05*(results.Tv(:,1) - 273.15);       % air_cp dry-air form

% Liquid sensible step, ambient datum -> film temperature, for the
% EVAPORATED stream only. A single-step form pricing that water as
% liquid@Ta -> vapour@Tv via h_fg(Tw) + Cp_v*(Tv - Ta) omits
% Cp_w*(Tw - Ta) and over-applies Cp_v across (Tw - Ta); the net
% omission is (Cp_w - Cp_v)*(Tw - Ta). Films run below Ta, so that
% term is negative and the ledger would over-charge the stream.
Tw_bar_t  = squeeze(mean(results.Tw,   2)).';                  % [nT x Ns]
Cpw_bar_t = squeeze(mean(results.Cp_w, 2)).';                  % [nT x Ns]

H_liq_step_t = sum( results.mevap_plate ...
                    .* (Cpw_bar_t - Cp_vapor_out_t) ...
                    .* (Tw_bar_t   - P.Ta), 2 );               % [nT x 1], W

% DATUM: ambient, matching the inlet term Q_fan_thermal and the latent
% reference used by Q_evaporation. See the datum note there.
H_air_out_t = P.mdot_da .* Cp_da_out_t .* (results.Tv(:,1) - P.Ta) ...
            + mdot_vap_out_t .* Cp_vapor_out_t .* (results.Tv(:,1) - P.Ta) ...
            + H_liq_step_t;

Q_vapor_out_ss  = mean(H_air_out_t(win));                   % window-mean rate, W
Q_vapor_out_end = H_air_out_t(end);                         % instantaneous rate at t_final, W
results.Q_vapor_out = Q_vapor_out_ss * P.t_operating * 3600;  



results.Y_final = Y(end,:).';   % needed for the rhs() re-eval above

results.total_evap_ts   = sum(results.mevap_plate, 2);
results.cum_distillate  = cumtrapz(t, results.mfw);
results.cum_solar_input = cumtrapz(t, results.I .* P.Ag);
results.GOR_ts          = (results.cum_distillate .* hfg_ref) ./ ...
                          max(results.cum_solar_input, eps);

results.TDS_out_per_plate  = squeeze(results.C(:,end,end));


[Qfan_elec, Qfan_therm] = fan_energy(P, t(end));
results.Q_feed = feed_thermal_energy(P);

results.Qfan_electrical_J = Qfan_elec;
results.Qfan_thermal_J    = Qfan_therm;

% ===================================================================
% PERFORMANCE METRICS
% -------------------------------------------------------------------
% Each metric is reported together with its numerator, denominator and
% integration limits in results.metrics, so that every quoted figure can
% be reconstructed from the tabulated powers.
%
% The window-consistent variant integrates production and solar input
% over identical limits, excluding the initial transient. The
% alternative variant extrapolates a windowed production rate to a full
% operating day while integrating the solar input over the whole
% simulated domain; it is reported for comparison but its numerator and
% denominator do not share a time base.
% ===================================================================
M = struct();

% ---- Window-consistent variant: same integration limits both sides ----
E_dist_win  = trapz(t(win), results.mfw(win)) * hfg_ref;              % [J]
E_solar_win = trapz(t(win), results.I(win) .* P.Ag);                  % [J]
M.window.t_start        = t(find(win,1));
M.window.t_end          = t(end);
M.window.E_distillate_J = E_dist_win;
M.window.E_solar_J      = E_solar_win;
M.window.GOR            = E_dist_win / max(E_solar_win, eps);

% ---- Mixed-basis variant, reported for comparison ----
M.mixed_basis.E_distillate_J = Edistill;   % mfw_daily * hfg_ref (extrapolated rate)
M.mixed_basis.E_solar_J      = Esolar;     % integrated over the full domain
M.mixed_basis.GOR            = Edistill / Esolar;
M.mixed_basis.note = ['numerator extrapolated from a windowed rate, ' ...
                      'denominator integrated over the full domain; ' ...
                      'time bases differ'];

M.hfg_ref_J_per_kg = hfg_ref;
M.timebase_note = ['the window variant uses identical integration ' ...
                   'limits for numerator and denominator'];

results.metrics = M;

% Headline GOR uses the window-consistent definition.
results.GOR             = M.window.GOR;
results.GOR_mixed_basis = M.mixed_basis.GOR;

% ---- Specific energy consumption ----
% The total specific energy consumption comprises the blower electrical
% work and the enthalpy carried in by the intake air above ambient. When
% the intake is preheated the latter dominates, so the two contributions
% are reported separately as well as in total.
vol_daily_m3 = results.mfw_daily/1000;
results.SEC_electrical_only = (Qfan_elec/3.6e6)               / max(vol_daily_m3, eps);
results.SEC_thermal_only    = (max(Qfan_therm,0)/3.6e6)       / max(vol_daily_m3, eps);
results.SEC_kWh_per_m3      = ((Qfan_elec + max(Qfan_therm,0))/3.6e6) / max(vol_daily_m3, eps);
results.metrics.SEC.electrical_kWh_per_m3 = results.SEC_electrical_only;
results.metrics.SEC.thermal_kWh_per_m3    = results.SEC_thermal_only;
results.metrics.SEC.total_kWh_per_m3      = results.SEC_kWh_per_m3;
results.metrics.SEC.volume_m3_per_day     = vol_daily_m3;

% ---- Glass-boundary losses, on the SAME basis as Esolar (full-domain
% integral when the sim already spans the full operating day; rate-
% extrapolated only for genuinely shorter runs). These are physical
% loss rates driven by Tg/Ta, not artifacts of the t=0 evaporation
% spike, so -- like Esolar -- there's no reason to exclude [0,t_relax]
% when the real simulated span already covers the operating day.
if sim_duration_hr >= 0.99*P.t_operating
    E_reflect_glass     = trapz(t, P.RF_g .* results.I) * P.Ag;
    E_loss_glass_conv_w = trapz(t, results.q_glass_conv);
    E_loss_glass_rad_w  = trapz(t, results.q_glass_rad_sky);
else
    dt_win = t(end) - t(find(win,1));
    E_reflect_glass     = trapz(t(win), P.RF_g .* results.I(win)) * P.Ag / dt_win * P.t_operating*3600;
    E_loss_glass_conv_w = trapz(t(win), results.q_glass_conv(win))        / dt_win * P.t_operating*3600;
    E_loss_glass_rad_w  = trapz(t(win), results.q_glass_rad_sky(win))     / dt_win * P.t_operating*3600;
end

results.E_reflect_glass     = E_reflect_glass;
results.E_loss_glass_conv_w = E_loss_glass_conv_w;
results.E_loss_glass_rad_w  = E_loss_glass_rad_w;

% ============ 4. Q_ground_loss (floor -> ground conduction) ============
Tbase_t = results.Tbase;   % [nT x 1]
q_ground_areal_t = (Tbase_t - P.Tground) ./ P.R_floor_ground;      % [nT x 1], W/m^2
Q_ground_t        = q_ground_areal_t .* P.Ap_floor;                 % [nT x 1], W

Q_ground_ss   = mean(Q_ground_t(win));
results.Q_ground_loss = Q_ground_ss * P.t_operating * 3600;

% ============ 5. Q_wallstrip_ambient (bottom-gap wall strip -> ambient) ============
Twallb_t = results.Twall(:,P.Ns);   % [nT x 1]
q_wallstrip_amb_areal_t = (Twallb_t - P.Ta) ./ P.R_wall_out;        % [nT x 1], W/m^2
Q_wallstrip_amb_t       = q_wallstrip_amb_areal_t .* P.Awall_bottom_gap;   % [nT x 1], W

Q_wallstrip_ss  = mean(Q_wallstrip_amb_t(win));
results.Q_wallstrip_loss = Q_wallstrip_ss * P.t_operating * 3600;


% Q_in must include the fan's THERMAL contribution (Qfan_therm), not
% only its electrical draw, and must subtract all three glass loss
% modes; omitting either breaks the closure.
results.Etotal = Esolar + Qfan_elec + Qfan_therm + results.Q_feed ...
               - min(E_reflect_glass,0) ...
               - min(E_loss_glass_conv_w,0) ...
               - min(E_loss_glass_rad_w,0) ...
               - min(results.Q_wall,0) ...
               - min(results.Q_feed_out,0) ...
               - min(results.Q_vapor_out,0) ...
               - min(results.Q_ground_loss,0) ...
               - min(results.Q_wallstrip_loss,0);

results.Etotal_thermal = Esolar+ Qfan_therm + results.Q_feed ...
                        - min(E_reflect_glass,0) ...
                        - min(E_loss_glass_conv_w,0) ...
                        - min(E_loss_glass_rad_w,0) ...
                        - min(results.Q_wall,0) ...
                        - min(results.Q_feed_out,0) ...
                        - min(results.Q_vapor_out,0) ...
                        - min(results.Q_ground_loss,0) ...
                        - min(results.Q_wallstrip_loss,0);

results.eta_thermal = Edistill /  results.Etotal_thermal* 100;

% BASIS NOTE -- READ BEFORE QUOTING THIS ALONGSIDE GOR.
% eta_thermal and steady_check.GOR_total_input are the SAME ratio,
% E_evap / E_in, on two different denominators, and they disagree by
% roughly a factor of two:
%   eta_thermal          : INTEGRATED over the window, denominator
%                          Etotal_thermal, which carries Qfan_therm --
%                          the AMBIENT-datum air enthalpy balance term.
%   GOR_total_input      : FINAL INSTANT, denominator solar + air lift +
%                          feed preheat on the T_coil datum.
% Both are valid on their own basis, but they are not interchangeable and must never be
% tabulated side by side without their bases. The denominator is exported
% here so the two can be reconciled in the report rather than guessed at.
results.Etotal_thermal_basis = ...
    'integrated over clean window; denominator carries ambient-datum Qfan_therm';
results.eta_thermal_basis    = results.Etotal_thermal_basis;


results.GOR_overall = Edistill / results.Etotal ;
results.eta_overall = results.GOR_overall * 100;


results.Edistill    =Edistill;
% Steady state energy balance results 
Y_end = Y(end,:) ;
y_end = Y(end,:).';
results.error = rhs(t(end),y_end,P) ;
S = unpack_state(Y_end, P);
[~, C_end, ~] = invert_holdup(S.M, S.Ms, S.Tw, P);

% ---- Q_solar (instantaneous, at final time) ----
[Ig, ~] = solar_irradiance(t(end), P);

% SOLAR INPUT CREDITED TO THE STEADY BALANCE
% Only the absorbed fraction of the incident beam enters the thermal
% balance. The remainder is disposed of by glazing reflection, by
% spillage beyond the plate footprint, and by reflection from the
% internal surfaces, none of which raises the temperature of any modelled
% control volume. The aperture and absorbed quantities are therefore
% reported separately, and the absorbed value is the one entering the
% efficiency denominator.
Q_solar_aperture = Ig * P.Ag;                            % [W] strikes the glass
Q_solar_absorbed = P.alpha_g*Ig*P.Ag ...                 % absorbed by the glass
                 + sum(results.q_solar_w_plate(end,:)) ...   % absorbed in the films
                 + sum(results.q_solar_p_plate(end,:));      % absorbed in plates/floor
Q_solar_unused   = Q_solar_aperture - Q_solar_absorbed;  % [W] never enters the model
Q_solar = Q_solar_absorbed;                              % [W] credited to the balance

% ---- Q_evap (instantaneous, summed across all plates) ----
Q_evap = 0;
mdot_evap_total = 0;
mdot_local_vec_end = P.mdot_da + P.mdot_vapor_in + ...
    flipud(cumsum(flipud(results.mevap_plate(end,:).')));
for k = 1:P.Ns
    Tw_k = S.Tw(k,:); C_k = C_end(k,:);
    Tv_k = S.Tv(k);   wv_k = S.wv(k);
    [~, hm_wv] = film_gap_coeffs(k, S, P, mdot_local_vec_end(k));
    Psat_w = psat_saline(Tw_k, C_k, P);
    Pv_k   = vapor_partial_pressure(wv_k, Tv_k, P);
    dP     = evap_driving_dP(Psat_w, Pv_k, P);
    mevap_k = hm_wv .* dP .* P.Mwater ./ (P.Rgas .* Tw_k);
    mdot_k  = sum(mevap_k) * P.dx_stage(k) * P.W;         % [kg/s], this stage
    hfg_k   = mean(latent_heat(Tw_k, C_k, P));            % LOCAL hfg, this stage
    Q_evap  = Q_evap + mdot_k * hfg_k;                    % [W]
    mdot_evap_total = mdot_evap_total + mdot_k;
end

% NOTE: the feed reconciliation below MUST precede every energy term
% that uses P.Tfeed. With tfeed_dynamic the RHS used a solved feed
% temperature, not the seed, and the feed thermal input, the GOR
% denominators and the still-only SEC all depend on it.
%
% The AIR inlet needs no reconciliation in this configuration: the
% electric air heater is specified by its outlet temperature, so
% P.T_air_in is what the RHS used, unchanged, at every step.

% =====================================================================
%  RECONCILE Tfeed WITH WHAT THE RHS ACTUALLY USED
% ---------------------------------------------------------------------
% With P.tfeed_dynamic the RHS evaluated the feed preheat HX outlet from
% the instantaneous brine at every step, so the value carried in P.Tfeed
% (the seed, = the supplied seawater temperature) is not what the
% cascade saw. Recompute it at the final state and publish that, so the
% HX block, the enthalpy accounting and the temperature report all
% describe the same stream.
if P.tfeed_dynamic
    S_end = unpack_state(Ysol(end,:).', P);
    mevap_end  = results.mevap_plate(end,:).';                  % [kg/s per stage]
    [Tfeed_end, tfinfo] = feed_preheat_T(S_end, mevap_end, P, t(end));
    % P.Tfeed IS NOT OVERWRITTEN HERE.
    % It is the SUPPLIED seawater temperature and it is what every feed
    % MASS expression is defined at (P.rhow_in, Gamma_feed, the recovery
    % denominator). Overwriting it with the preheated value would evaluate
    % the supplied seawater at the hot-side density, understating the
    % water fed by rho(Tfeed_in)/rho(Tfeed_out) and inflating every
    % recovery figure by the same factor. The solved cascade inlet
    % temperature is published as
    % results.feedpre.Tfeed_out and results.P.Tfeed_cascade_in; consumers
    % that need the hot-side state read those, not P.Tfeed.
    P.Tfeed_cascade_in         = Tfeed_end;
    results.P.Tfeed_cascade_in = Tfeed_end;
    results.feedpre.dynamic    = true;
    results.feedpre.Tfeed_out  = Tfeed_end;
    results.feedpre.Tfeed_in   = P.Tfeed_in;
    results.feedpre.dT         = Tfeed_end - P.Tfeed_in;
    results.feedpre.capped     = tfinfo.capped;
    results.feedpre.pinched    = tfinfo.pinched;
    results.feedpre.mdot_vap   = tfinfo.mdot_vap;
    results.feedpre.dT_raw     = tfinfo.dT_raw;
    results.feedpre.dT_pinch   = tfinfo.dT_pinch;
    results.feedpre.Cden       = tfinfo.Cden;
    results.feedpre.Q_avail    = tfinfo.Q_avail;
    results.feedpre.Q_absorbed = tfinfo.Q_absorbed;
    results.feedpre.T_cond     = P.bl.T_sat_vap;
    results.feedpre.T_cap      = P.bl.T_sat_vap - P.dT_pinch_HX;
    results.feedpre.seed       = P.Tfeed_in;
    % ---- Three-regime reporting (target / pinch / supply) ----
    % Which limit actually bound is a physical statement: a run that
    % crosses between regimes mid-solve has two qualitatively different
    % exchangers in one trajectory, and a single-number f_cond quoted
    % from it describes neither. P.Tfeed_regime is read by
    % preheater_heater_awg().
    results.feedpre.regime           = tfinfo.regime;
    results.feedpre.T_target         = tfinfo.T_target;
    results.feedpre.Q_surplus        = tfinfo.Q_surplus;
    results.feedpre.f_cond           = tfinfo.f_cond;
    results.feedpre.mdot_vap_surplus = tfinfo.mdot_vap_surplus;
    P.Tfeed_regime                   = tfinfo.regime;
    results.P.Tfeed_regime           = tfinfo.regime;
else
    results.feedpre.dynamic    = false;
    results.feedpre.Tfeed_out  = P.Tfeed;
    results.feedpre.Tfeed_in   = P.Tfeed_in;
end

% Air-side inlet record: fixed, for symmetry in the reporting blocks.
results.tfan.dynamic = false;
results.tfan.Tfan_in = P.T_air_in;
results.tfan.dT      = P.T_air_in - P.T_coil;   % rise across the air heater


% ---- Q_fan_thermal (instantaneous power, no integration needed) ----
% ---- Condenser reheat carried in by the return air ----
% Evaluated against an ambient datum for comparability with the loss
% terms. In the closed loop this flux originates at the AWG condenser
% and is therefore recirculated rather than externally supplied; it is
% classified accordingly in the accounting below.
% CONVENTION: evaluated per component (dry air and water vapour
% separately), identical to the construction of the outgoing stream
% H_air_out_t. This matters because moist_air_cp returns a specific heat
% per kilogram of MOIST air,
%
%     Cp_m = (Cp_da + w*Cp_v) / (1 + w)
%
% so it must be multiplied by the moist mass flow, not by the dry-air
% mass flow. Pairing Cp_m with mdot_da would understate the incoming
% enthalpy by a factor (1+w) relative to the per-component outgoing
% stream, leaving a spurious imbalance across the air loop.
w_fan_in     = P.wv_fan_in;
mdot_vap_in  = P.mdot_da * w_fan_in;                       % [kg/s]

% ---- DATUM: T_coil, NOT ambient. ----------------------------------
% In the CLOSED loop no air is taken from ambient, so an ambient datum
% measures nothing that physically happens. The loop air leaves the AWG
% at T_coil and the electric air heater raises it to T_air_in; that rise
% is the external energy actually injected into the air stream.
%
% An ambient datum measures (T_air_in - Ta), while fan_energy() -- the
% integrated report of the SAME stream -- uses (T_air_in - T_coil). The
% two spans are unrelated, so the two blocks disagree on one unit, and
% an ambient-datum value carried into the thermal-share diagnostic
% understates the air heater's share.
%
% NOTE ON THE STILL ENERGY BALANCE: this term is paired with the
% outgoing air enthalpy H_air_out_t, which is referenced to the SAME
% datum below. The datum cancels in (Q_in - Q_out), so the closure is
% insensitive to the choice, while the reported SHARES are not.
% DATUM FOR THE CLOSURE BLOCK: AMBIENT, and it must stay ambient.
%
% A T_coil datum for this term alone would match the reported "air
% heater" share to the unit's real duty but would break the still energy
% balance: the datum does NOT cancel between inlet and outlet here,
% because the outgoing air carries MORE vapour mass than the incoming air
% (that is the whole point of the cascade), and the latent term
% Q_evaporation is referenced to the film state. Shifting only the
% sensible datum leaves a closure residual of order 9% instead of 0.4%.
%
% The two purposes are therefore separated:
%   - THIS term, ambient datum, pairs with H_air_out_t and closes the
%     still energy balance. It is a balance term, not a unit duty.
%   - The air heater's REAL external duty is D.Q_air_heater, computed on
%     the T_coil datum in preheater_heater_awg(), and that is what the
%     thermal-share block reports.
Cp_da_in     = 1005 + 0.05*(P.Tfan_in - 273.15);
Cp_vap_in    = 1860 + 0.12*(P.Tfan_in - 273.15);
Q_fan_thermal = ( P.mdot_da   * Cp_da_in ...
                + mdot_vap_in * Cp_vap_in ) * (P.Tfan_in - P.Ta);   % [W]


% ---- Q_feed_thermal (instantaneous, feed enters continuously) ----
mdot_feed = P.Vfeed * P.rhow_in / 1000 / (P.t_operating * 3600);
% UPPER LIMIT IS THE CASCADE INLET, NOT THE SUPPLY TEMPERATURE.
% This is a BOUNDARY term of the STILL control volume: the feed physically
% crosses that boundary at P.Tfeed_cascade_in (the preheat HX sits
% upstream of the cascade), so its enthalpy relative to ambient must be
% evaluated there. Using the supply temperature instead removes the
% entire preheat enthalpy from Q_in and breaks the closure.
%
% Whether this enthalpy is "external" is a separate, SEC-level question:
% it is internal recycle (evaporator vapour raised from the plant's own
% brine, already charged to the electric heater) and the thermal-share
% block labels it as such. That labelling does not remove it from the
% still's energy balance -- the joules cross the boundary either way.
T_range = linspace(P.Ta, P.Tfeed_cascade_in, 50);
Cp_vals = water_cp(T_range, P.TDSfeed, P);
delta_h = trapz(T_range, Cp_vals);
Q_feed_thermal = mdot_feed * delta_h;                     % [W]
Twallb_iE = results.Twall(end, P.Ns);
Tbase_iE = Tbase_t (end) ;
% ---- Glass boundary losses (instantaneous, at final time) ----
Tsky_end = 0.0552*P.Ta^1.5;
hr_gsky_end = P.eps_g*P.sigma*(S.Tg^2 + Tsky_end^2)*(S.Tg + Tsky_end);

Q_reflect_glass    = P.RF_g * Ig * P.Ag;                       % Mode 3: never absorbed
Q_glass_conv_loss  = P.hwind * (S.Tg - P.Ta) *  P.Ag;           % Mode 1: convective, to ambient
Q_glass_rad_loss   = hr_gsky_end * (S.Tg - Tsky_end) *  P.Ag;   % Mode 2: radiative, to sky

Q_ground_loss_end    = (Tbase_iE - P.Tground) / P.R_floor_ground * P.Ap_floor;
Q_wallstrip_loss_end = (Twallb_iE - P.Ta)      / P.R_wall_out      * P.Awall_bottom_gap;

% ---- Boundary streams ----
% SIGN CONVENTION: positive = energy LEAVING the still, evaluated against
% an ambient datum (T = P.Ta).
%
% These streams are not losses by construction. Whenever the intake air
% is preheated well above ambient and the evaporative load is large, the
% cascade operates BELOW ambient temperature: the films and the exhaust
% are then colder than the surroundings, and the exhaust enthalpy,
% glazing convection and wall conduction all reverse sign and become
% gains. A formulation that partitioned these terms by sign and treated
% only the positive ones as outgoing would move the reversed streams into
% the input side and inflate the apparent efficiency. They are therefore
% summed with their actual signs.
%
% Q_reflect_glass is deliberately EXCLUDED. Q_solar above is the absorbed
% beam power, from which glazing reflection has already been removed;
% listing it again as an outgoing stream would subtract it twice. It is
% reported separately in the solar closure block.
%
% TIME BASIS: every term below is evaluated at the final stored sample.
% Under time-varying irradiance a balance that mixed window-averaged and
% instantaneous quantities would not close, since the two differ by the
% drift of the forcing across the averaging window.
loss_terms = [Q_glass_conv_loss, Q_glass_rad_loss, ...
              Q_wall_end, Q_feed_out_end, Q_vapor_out_end, ...
              Q_ground_loss_end, Q_wallstrip_loss_end];

Q_out_total = sum(loss_terms);             % signed: net energy leaving via the boundary
Q_gain_total = -sum(min(loss_terms,0));    % magnitude of the reversed (incoming) streams,
                                           % reported for transparency only

Q_in_total = Q_solar + Q_fan_thermal + Q_feed_thermal;   % deliberately supplied energy

% ---- Thermal storage ----
% The device is forced by a measured irradiance record and therefore
% never attains a true steady state: the sensible heat content of the
% glazing, films, plates, walls, floor and vapour gaps changes
% continuously over the diurnal cycle. That rate of change is a genuine
% term in the instantaneous power balance and is evaluated here by
% backward differencing the stored state at the final sample. Omitting
% it would leave an unexplained residual of the same order as the
% storage rate itself.
Q_storage_total = instantaneous_storage_rate(results, P);

% ---- Closure ----
% Q_in = Q_out + Q_evap + Q_storage, with the residual reported so that
% the balance can be verified from the printed quantities.
Q_balance_residual = Q_in_total - Q_out_total - Q_evap - Q_storage_total;

% =====================================================================
% DOWNSTREAM CHAIN -- CLOSED AIR LOOP
% ---------------------------------------------------------------------
%   still -> brine -> BRINE EVAPORATOR -> vapour -> FEED PREHEAT HX
%                                                        |
%   still exhaust -> AWG -> ELECTRIC AIR HEATER -> still | (loop closes)
%                     ^                                  |
%                     +------------- surplus vapour <----+
%
% The AIR side of this block is still evaluated once, after the PDE
% solve: T_coil is a fixed setpoint and the air heater is specified by
% outlet temperature, so neither depends on the solution.
%
% The FEED side is NOT feed-forward. The feed preheat HX returns to the
% cascade, and that recycle was already closed INSIDE the integration by
% feed_preheat_T(). What is computed here is the same exchanger
% re-evaluated at the final state for REPORTING; it does not re-close
% anything.
%
% IMPORTANT DISTINCTION, and the reason results.mfw is NOT the product:
% results.mfw is the CASCADE EVAPORATION rate -- water transferred into
% the air stream. The loop is closed, so in steady operation essentially
% all of it is condensed at the coil and collected. It is nonetheless
% not identical to the product, for two reasons that have nothing to do
% with atmospheric rejection:
%   (a) the product ALSO contains the evaporator vapour condensed in the
%       feed preheat HX and the AWG, which never entered the air stream;
%   (b) the loop holds a moisture INVENTORY, so over a finite window the
%       evaporation and the condensation differ by its rate of change.
% The system product is assembled explicitly in results.product below,
% and every SEC figure uses it.
% =====================================================================
Tv_exit  = results.Tv(end,1);            % [K]  air leaving gap 1
wv_exit  = results.wv(end,1);            % [kg/kg]

% ---- Dew point of the cascade exhaust, from its vapour partial pressure ----
Pv_exit  = wv_exit * 101325 / (0.622 + wv_exit);
T_dew    = dew_point_from_Pv(Pv_exit);

% =====================================================================
%  COIL TEMPERATURE -- 'auto' mode
% ---------------------------------------------------------------------
% Set HERE, after the solve, because the exhaust dew point is a solved
% quantity.
%
% CAUTION -- THIS IS NOT A CLOSED LOOP IN T_coil. P.wv_fan_in is
% w_sat(T_coil) (see build_parameters), and it sets P.mdot_da and the
% gap inlet boundary condition, so T_coil DOES feed the solve. Assigning
% a new T_coil here therefore leaves the trajectory that produced the
% dew point inconsistent with the coil reported beside it, unless
% the run is iterated to a fixed point externally. 'auto' mode is a
% one-pass estimate on that basis; 'fixed' mode has no such issue.
%
% A colder coil is not free: the lift (T_cond - T_coil) grows and the
% COP falls, so the compressor draws more. In this plant the heater
% dominates SEC, so the extra condensate is worth far more than the
% extra compressor work -- but the trade is real and is reported.
% =====================================================================
P.T_coil_requested = P.T_coil;          % what the user asked for
P.T_coil_frosted   = false;
if strcmp(P.coil_mode,'auto')
    T_coil_auto = T_dew - P.dT_coil_margin;
    if T_coil_auto < P.T_coil_floor
        P.T_coil_frosted = true;
        warning('ZLDD:coilFrostGuard', ...
            ['Auto coil temperature %.2f K is below the frost guard %.2f K ' ...
             'and has been clamped. The model has NO frost logic, so a ' ...
             'colder coil would be reported but not simulated. The ' ...
             'dew-point margin will be smaller than the %.1f K requested.'], ...
            T_coil_auto, P.T_coil_floor, P.dT_coil_margin);
        T_coil_auto = P.T_coil_floor;
    end
    P.T_coil = T_coil_auto;
    results.P.T_coil = P.T_coil;        % keep the published copy in step
end
results.P.T_coil_requested  = P.T_coil_requested;
results.P.T_coil_frosted    = P.T_coil_frosted;
results.P.coil_mode         = P.coil_mode;
results.P.dT_coil_margin    = P.dT_coil_margin;

wv_coil  = humidity_ratio_from_RH(1.0, P.T_coil, P);   % [kg/kg] saturated at the coil

% ---- Downstream unit balances ----
D = preheater_heater_awg(results.brine_out, Tv_exit, wv_exit, wv_coil, P);

% =====================================================================
% INTEGRATED PLANT PRODUCT  (replaces final-instant x t_operating)
% ---------------------------------------------------------------------
% D above is the final-instant snapshot: it sizes the units, which is
% what procurement needs. It is NOT a daily yield. Multiplying its
% mdot_total by t_operating credits the plant with the closing rate for
% the whole window, and under diurnal forcing the closing rate is not
% the window mean -- the plates are still hot after solar noon, so the
% extrapolation runs high and can push recovery past 100%.
%
% preheater_heater_awg is feed-forward and cheap, so the honest form is
% to evaluate it at every sample in the CLEAN window (post-relaxation,
% the same window mfw_daily uses) and integrate. Production before
% t_relax is treated as zero, exactly as total_distillate_clean does.
%
% The brine struct passed at each sample is built on the SAME basis as
% results.brine_out -- outlet-face rho*u*delta*W with w_salt = C/rho --
% so this is a time series of the identical stream, not a new one.
t_all   = results.t(:);
win_cln = t_all > results.spike_relax_time;
i_cln   = find(win_cln);

mdot_prod_t = zeros(numel(i_cln),1);   % [kg/s] total collected condensate
mdot_vap_t  = zeros(numel(i_cln),1);   % [kg/s] heater vapour
mdot_air_t  = zeros(numel(i_cln),1);   % [kg/s] air-side recovery
W_ext_t     = zeros(numel(i_cln),1);   % [W]    heater + AWG + air heater

% ---- PER-STREAM series, so that every reported stream shares ONE basis --
% The five component streams are integrated over the same clean window as
% the total. On a final-instant basis they would not satisfy the stream
% table's own identity (S12 = S7+S7b+S7c+S8+S10) against a
% window-integrated total: the two sides would differ by the diurnal
% extrapolation error, appearing as a discrepancy between PART A and
% PARTS B/C.
mdot_pre_t   = zeros(numel(i_cln),1);  % [kg/s] S7   feed-HX condensate
mdot_ahx_t   = zeros(numel(i_cln),1);  % [kg/s] S7b  air-preheater condensate
mdot_rec_t   = zeros(numel(i_cln),1);  % [kg/s] S7c  recuperator drain
mdot_avap_t  = zeros(numel(i_cln),1);  % [kg/s] S8   AWG vapour condensate
mdot_aair_t  = zeros(numel(i_cln),1);  % [kg/s] S10  AWG air condensate
mdot_rej_t   = zeros(numel(i_cln),1);  % [kg/s] S11  moisture rejected
mdot_slur_t  = zeros(numel(i_cln),1);  % [kg/s] S6   salt slurry
mdot_salt_t  = zeros(numel(i_cln),1);  % [kg/s] S4s  salt (conserved)
mdot_brn_t   = zeros(numel(i_cln),1);  % [kg/s] S4   cascade brine, total
mdot_brnw_t  = zeros(numel(i_cln),1);  % [kg/s] S4w  cascade brine, water
mdot_evap_t  = zeros(numel(i_cln),1);  % [kg/s] S3   cascade evaporation
W_heat_t     = zeros(numel(i_cln),1);  % [W]    electric heater
W_awg_t      = zeros(numel(i_cln),1);  % [W]    AWG compressor
W_aht_t      = zeros(numel(i_cln),1);  % [W]    electric air heater
f_cond_t     = zeros(numel(i_cln),1);  % [-]    feed-HX condensed fraction

mdot_brine_all = results.rho_out(:,end) .* results.u_out(:,end) .* ...
                 results.delta_out(:,end) .* P.W;             % [kg/s]

for jj = 1:numel(i_cln)
    ii = i_cln(jj);
    b  = struct();
    b.mdot   = max(mdot_brine_all(ii), 0);
    b.C      = results.C_out(ii,end);
    b.rho    = results.rho_out(ii,end);
    b.T      = results.Tw_out(ii,end);
    b.w_salt = b.C / max(b.rho, eps);
    b.basis  = 'final instant (windowed sample)';   % silences the basis guard

    Dj = preheater_heater_awg(b, results.Tv(ii,1), results.wv(ii,1), wv_coil, P);

    % Must carry the SAME four streams as the final-instant assembly
    % below, or kg_per_day and kg_per_day_int will disagree without
        % any check catching it.
    mdot_vap_t(jj)  = Dj.mdot_cond_preheater ...
                    + Dj.mdot_cond_airheat ...
                    + Dj.mdot_cond_awg_vapour;
    % Recuperator drain is air-side moisture, same origin as the coil
    % condensate: the two units split one stream between them.
    mdot_air_t(jj)  = Dj.mdot_cond_air + Dj.mdot_cond_recup;
    mdot_prod_t(jj) = mdot_vap_t(jj) + mdot_air_t(jj);
    W_ext_t(jj)     = Dj.W_heater + Dj.W_awg + Dj.W_air_heater;

    % ---- Per-stream, same sample, same Dj: nothing here is a second
    % evaluation of anything above, so the component integrals sum to the
    % total integral exactly rather than to within a solver tolerance.
    mdot_pre_t(jj)  = Dj.mdot_cond_preheater;
    mdot_ahx_t(jj)  = Dj.mdot_cond_airheat;
    mdot_rec_t(jj)  = Dj.mdot_cond_recup;
    mdot_avap_t(jj) = Dj.mdot_cond_awg_vapour;
    mdot_aair_t(jj) = Dj.mdot_cond_air;
    mdot_slur_t(jj) = Dj.mdot_brine_out;
    mdot_salt_t(jj) = Dj.mdot_salt;
    mdot_brn_t(jj)  = b.mdot;
    mdot_brnw_t(jj) = b.mdot * (1 - b.w_salt);
    W_heat_t(jj)    = Dj.W_heater;
    W_awg_t(jj)     = Dj.W_awg;
    W_aht_t(jj)     = Dj.W_air_heater;
    f_cond_t(jj)    = Dj.f_cond;
    if isfield(Dj,'mdot_moisture_rejected')
        mdot_rej_t(jj) = Dj.mdot_moisture_rejected;
    end
    % Cascade evaporation at this sample: the loop is closed, so what the
    % cascade evaporates is what the recuperator and the coil together
    % condense. Taken from the SAME Dj, so S3 = S7c + S10 holds sample by
    % sample and therefore after integration too.
    mdot_evap_t(jj) = Dj.mdot_cond_recup + Dj.mdot_cond_air;
end

t_cln = t_all(i_cln);
if numel(t_cln) > 1
    results.product.kg_per_day_int = trapz(t_cln, mdot_prod_t);   % [kg] over the clean window
    results.product.vap_kg_int     = trapz(t_cln, mdot_vap_t);
    results.product.air_kg_int     = trapz(t_cln, mdot_air_t);
    results.product.W_ext_mean     = trapz(t_cln, W_ext_t) / (P.t_operating*3600);
    results.product.mdot_total_int = results.product.kg_per_day_int / ...
                                     (P.t_operating*3600);   % [kg/s] day average
    results.product.int_basis      = 'true integral over clean window (diurnal)';
    results.product.int_valid      = true;

    % ---- WINDOW-MEAN RATES for every stream ----------------------------
    % Published as RATES (kg/s), not daily totals, so that the report's
    % g/s and kg/day columns are the same quantity in two units and
    % cannot disagree. dt_win is the integration span actually used.
    % NORMALISED BY THE FULL OPERATING WINDOW, not by the clean-window
    % span. The report multiplies every rate by t_operating to get a
    % kg/day column; dividing by the shorter clean span here would make
    % rate*t_operating exceed the integral it came from, and the two
    % columns of the same row would disagree. Production before
    % t_relax is counted as zero -- the same convention mfw_daily uses.
    dt_win = max(t_cln(end)-t_cln(1), eps);   % [s] diagnostic only
    t_op_s_wm = P.t_operating*3600;
    wm = @(y) trapz(t_cln, y) / t_op_s_wm;
    results.wmean.mdot_preheater  = wm(mdot_pre_t);
    results.wmean.mdot_airheat    = wm(mdot_ahx_t);
    results.wmean.mdot_recup      = wm(mdot_rec_t);
    results.wmean.mdot_awg_vapour = wm(mdot_avap_t);
    results.wmean.mdot_awg_air    = wm(mdot_aair_t);
    results.wmean.mdot_lost_air   = wm(mdot_rej_t);
    results.wmean.mdot_vap        = wm(mdot_vap_t);
    results.wmean.mdot_brine_out  = wm(mdot_slur_t);
    results.wmean.mdot_salt       = wm(mdot_salt_t);
    results.wmean.mdot_brine      = wm(mdot_brn_t);
    results.wmean.mdot_brine_water= wm(mdot_brnw_t);
    results.wmean.mdot_evap       = wm(mdot_evap_t);
    results.wmean.W_heater        = wm(W_heat_t);
    results.wmean.W_awg           = wm(W_awg_t);
    results.wmean.W_air_heater    = wm(W_aht_t);
    % f_cond is an INTENSIVE ratio: its window mean must be duty-weighted,
    % not time-weighted, or a period of near-zero vapour flow carries the
    % same weight as full production.
    results.wmean.f_cond          = trapz(t_cln, f_cond_t .* mdot_vap_t) / ...
                                    max(trapz(t_cln, mdot_vap_t), eps);
    results.wmean.dt_win          = dt_win;
    results.wmean.valid           = true;
else
    results.product.kg_per_day_int = NaN;
    results.product.int_basis      = 'INTEGRATION FAILED: clean window has < 2 samples';
    results.product.int_valid      = false;
    results.wmean.valid            = false;
end
results.product.mdot_prod_t = mdot_prod_t;
results.product.t_prod      = t_cln;

% ---- Aliases, so that the reporting functions below read cleanly ----
mdot_cond = D.mdot_cond_air;                               % [kg/s] air-side recovery
Q_evap_coil  = D.Q_awg_total;                              % [W]
W_compressor = D.W_awg;                                    % [W]
Q_condenser  = D.Q_awg_reject;                             % [W]
T_cond_awg   = D.T_cond_awg;                               % [K]
COP_lift     = D.COP;                                      % [-]

W_compressor_fixCOP = Q_evap_coil / max(P.COP_awg, eps);   % [W] comparator
Q_condenser_fixCOP  = Q_evap_coil + W_compressor_fixCOP;   % [W] comparator

% ---- Cascade evaporation, for the air-side moisture audit ----
mdot_evap_still = steady_check_mdot_evap_total(results);   % [kg/s]

% CLOSED LOOP: the air-side moisture MUST balance. Nothing is rejected
% to atmosphere, so every gram the cascade evaporates is condensed at
% the coil (in steady state) and the residual below is a genuine closure
% check rather than a reported loss.
%
% mdot_moisture_rejected is carried as a field for the sweep drivers and
% the reporting blocks, and it is a STRUCTURAL ZERO: no air leaves the
% loop, so no moisture leaves with it. It is NOT P.mdot_da*wv_coil --
% that expression measures a discharge stream, which this flowsheet does
% not have, and would report a large fictitious loss.
mdot_moisture_rejected = 0;                                         % [kg/s] closed loop: nothing leaves
air_moisture_residual  = P.mdot_da*(wv_exit - P.wv_fan_in) - mdot_evap_still;  % [kg/s] must be ~0

results.awg.Tv_exit                = Tv_exit;
results.awg.wv_exit                = wv_exit;
results.awg.T_dew_return           = T_dew;
results.awg.mdot_condensate        = mdot_cond;
results.awg.Q_evaporator           = Q_evap_coil;
results.awg.W_compressor           = W_compressor;
results.awg.Q_condenser            = Q_condenser;
results.awg.Q_reheat_required      = 0;      % reheat is supplied by the SEPARATE electric air heater,
results.awg.reheat_margin          = NaN;    % not by this unit's condenser -- so no margin screen applies
results.awg.loop_moisture_residual = air_moisture_residual;
results.awg.mdot_moisture_rejected = mdot_moisture_rejected;
results.awg.T_coil                 = P.T_coil;
results.awg.wv_coil                = wv_coil;
results.awg.COP                    = COP_lift;
results.awg.COP_fixed              = P.COP_awg;
results.awg.T_cond                 = T_cond_awg;
results.awg.lift                   = T_cond_awg - P.T_coil;
results.awg.eta_II                 = P.eta_II_awg;
results.awg.W_compressor_fixCOP    = W_compressor_fixCOP;
results.awg.Q_condenser_fixCOP     = Q_condenser_fixCOP;
results.awg.reheat_margin_fixCOP   = NaN;
results.awg.mdot_evap_still        = mdot_evap_still;
results.awg.Q_solar_absorbed       = Q_solar_absorbed;

% ---- Full downstream record ----
results.downstream = D;

% =====================================================================
% SYSTEM PRODUCT
% ---------------------------------------------------------------------
% Three condensate streams are collected and combined (the plasma
% chamber that follows is a treatment step with no mass or energy effect
% and is deliberately not modelled):
%   1. preheater condensate  -- part of the heater vapour, condensed
%                               against the incoming air
%   2. AWG vapour condensate -- the remainder of the heater vapour,
%                               condensed essentially completely
%   3. AWG air condensate    -- the fraction of the cascade's moisture
%                               recovered from the humid air
% Streams 1 and 2 together are the whole heater vapour stream, so the
% product reduces to (heater vapour) + (air-side recovery).
% =====================================================================
% FOUR condensate streams. The air-preheater condensate (2b) is drained
% at that unit and never reaches the AWG, so omitting it here would
% delete real product from every downstream figure.
results.product.mdot_preheater = D.mdot_cond_preheater;    % [kg/s] feed HX
results.product.mdot_airheat   = D.mdot_cond_airheat;      % [kg/s] air preheater
results.product.mdot_recup     = D.mdot_cond_recup;        % [kg/s] recuperator drain
results.product.mdot_awg_vapour= D.mdot_cond_awg_vapour;   % [kg/s]
results.product.mdot_awg_air   = D.mdot_cond_air;          % [kg/s]
results.product.mdot_total     = D.mdot_cond_preheater ...
                               + D.mdot_cond_airheat ...
                               + D.mdot_cond_recup ...
                               + D.mdot_cond_awg_vapour ...
                               + D.mdot_cond_air;          % [kg/s]
results.product.kg_per_day     = results.product.mdot_total * P.t_operating*3600;
results.product.mdot_lost_air  = mdot_moisture_rejected;   % [kg/s] rejected to atmosphere

% =====================================================================
% ONE BASIS FOR EVERY REPORTED NUMBER
% ---------------------------------------------------------------------
% Everything assigned above is the FINAL-INSTANT snapshot: it is what
% sizes the units, and PART C must keep it for exactly that reason. It
% is NOT a yield. Under diurnal forcing the closing rate is above the
% window mean -- the plates are still hot after solar noon -- so
% snapshot x t_operating overstates production.
%
% The snapshot and the window integral must not share a label: if PART A
% reported the integral while PART B, PART C and the SEC denominator read
% the snapshot, one name would carry two numbers with no way for a reader
% to tell them apart.
%
% The snapshot is therefore kept under its own name: it lives in
% results.snapshot.* and PART C's unit specifications read from there.
% Every PERFORMANCE field below is overwritten with the window mean, so
% each stream is one quantity on one basis and the stream-table
% identities close on it.
results.snapshot.mdot_preheater  = results.product.mdot_preheater;
results.snapshot.mdot_airheat    = results.product.mdot_airheat;
results.snapshot.mdot_recup      = results.product.mdot_recup;
results.snapshot.mdot_awg_vapour = results.product.mdot_awg_vapour;
results.snapshot.mdot_awg_air    = results.product.mdot_awg_air;
results.snapshot.mdot_total      = results.product.mdot_total;
results.snapshot.kg_per_day      = results.product.kg_per_day;
results.snapshot.mdot_lost_air   = results.product.mdot_lost_air;
results.snapshot.basis           = 'final instant (unit sizing)';

if isfield(results,'wmean') && isfield(results.wmean,'valid') && results.wmean.valid
    results.product.mdot_preheater  = results.wmean.mdot_preheater;
    results.product.mdot_airheat    = results.wmean.mdot_airheat;
    results.product.mdot_recup      = results.wmean.mdot_recup;
    results.product.mdot_awg_vapour = results.wmean.mdot_awg_vapour;
    results.product.mdot_awg_air    = results.wmean.mdot_awg_air;
    results.product.mdot_lost_air   = results.wmean.mdot_lost_air;
    results.product.mdot_total      = results.product.mdot_total_int;
    results.product.kg_per_day      = results.product.kg_per_day_int;
    results.product.basis           = results.product.int_basis;

    % SELF-CHECK. The five component means are integrals of the same Dj
    % samples that built the total, so this must close to round-off. A
    % failure here means a stream was added to one list and not the
    % other -- the exact inconsistency this check exists to catch.
    resid_basis = results.product.mdot_total - ( ...
        results.product.mdot_preheater + results.product.mdot_airheat + ...
        results.product.mdot_recup + results.product.mdot_awg_vapour + ...
        results.product.mdot_awg_air);
    if abs(resid_basis) > 1e-9*max(results.product.mdot_total, eps)
        warning('ZLDD:ProductBasis', ...
            ['Window-mean product streams do not sum to the window-mean ' ...
             'total (residual %.3e kg/s). A stream is missing from one ' ...
             'of the two lists in the clean-window loop.'], resid_basis);
    end
else
    results.product.basis = 'final-instant rate x t_operating (EXTRAPOLATED)';
end

% ---- Exported so that sweep drivers do not have to recompute them ----
% water_density() and the humidity helpers are LOCAL to this file and are
% not on the path of a calling script, so any quantity a driver needs
% must be published here rather than rebuilt outside.
results.product.feed_water_kg_per_day = ...
    P.Vfeed*P.rhow_in/1000 - P.Vfeed*P.TDSfeed/1000;
% ---- ATMOSPHERIC HARVEST: IDENTICALLY ZERO IN THE CLOSED LOOP ----
% The loop is airtight and admits no make-up, so no atmospheric water
% enters the plant at any point and ALL product is seawater-derived.
%
% This is a structural property of the flowsheet, not a bookkeeping
% choice: with a once-through intake ambient moisture would be a real
% product route, whereas in a closed loop the term vanishes by
% construction. The seawater/harvest split still appears in the headline
% block and the stream table, with the harvest route shown as an explicit
% zero.
%
% The expression is kept in its full form rather than hard-coded to 0 so
% that the reason is legible: w_fan_in IS w_sat(T_coil) in the closed
% loop, so the difference is exactly zero by the definition of the inlet
% state -- not by an assumption applied here.
results.product.mdot_harvest = max(P.mdot_da*(P.wv_fan_in - wv_coil), 0);  % [kg/s] == 0, closed loop
results.product.harvest_kg_per_day = results.product.mdot_harvest * P.t_operating*3600;
results.product.seawater_kg_per_day = ...
    results.product.kg_per_day - results.product.harvest_kg_per_day;
results.product.harvest_is_structural_zero = true;
% Inherits whichever basis the block above selected, so it cannot
% disagree with the recovery printed in PART A.
results.product.recovery_pct   = 100 * results.product.kg_per_day / ...
    max(P.Vfeed*P.rhow_in/1000 - P.Vfeed*P.TDSfeed/1000, eps);


% ---- Energy crossing the boundary of the modelled component ----
% The system boundary of this study encloses the still only. Absorbed
% solar radiation and blower work cross it; the AWG's own electrical
% demand lies outside and is not attributed to this model.
W_blower = results.Qfan_electrical_J / (P.t_operating*3600);   % [W]
results.awg.W_blower       = W_blower;
results.awg.E_still_inputs = Q_solar_absorbed + W_blower;

% =====================================================================
% SYSTEM-LEVEL SPECIFIC ENERGY CONSUMPTION (still + heater + preheater + AWG)
% ---------------------------------------------------------------------
% The control volume encloses the whole plant. FOUR streams of purchased
% energy cross it, the air heater being a separate electrically driven
% unit rather than a recipient of evaporator vapour:
%
%     W_ext = W_blower + W_heater + W_awg + W_air_heater
%
% The intake-air preheat is NOT a purchased stream. It is supplied by the
% electric heater's vapour, whose electrical draw is already counted in
% W_heater, so adding the preheat would double-count the same joules: it
% is an internal transfer between two units that both lie inside the
% boundary.
%
% The denominator is the COLLECTED PRODUCT, not the cascade evaporation.
% These differ because the product also contains evaporator vapour
% condensed in the feed preheat HX and the AWG, which never entered the
% air stream at all, and because the closed loop holds a moisture
% inventory whose rate of change separates the two over a finite window.
% Using mfw_daily here would describe the cascade, not the plant.
% =====================================================================
% WINDOW MEANS, to match the denominator. W_blower is a window mean
% (Qfan_electrical_J / t_operating), so the other three are taken on the
% same basis; a final-instant value for any of them would place an
% instantaneous numerator over an integrated denominator.
if isfield(results,'wmean') && isfield(results.wmean,'valid') && results.wmean.valid
    W_heater_el = results.wmean.W_heater;                      % [W] brine evaporator
    W_awg_el    = results.wmean.W_awg;                         % [W] AWG compressor
    W_airht_el  = results.wmean.W_air_heater;                  % [W] electric air heater
    results.system.W_basis = 'window mean over the clean window';
else
    W_heater_el = results.downstream.W_heater;                 % [W] brine evaporator
    W_awg_el    = results.downstream.W_awg;                    % [W] AWG compressor
    W_airht_el  = results.downstream.W_air_heater;             % [W] electric air heater
    results.system.W_basis = 'final instant (EXTRAPOLATED)';
end
W_external  = W_blower + W_heater_el + W_awg_el + W_airht_el;  % [W]
mdot_prod   = results.product.mdot_total;                      % [kg/s] COLLECTED

results.system.W_blower        = W_blower;
results.system.W_heater        = W_heater_el;
results.system.W_compressor    = W_awg_el;
results.system.W_air_heater    = W_airht_el;
results.system.W_external      = W_external;
results.system.mdot_product    = mdot_prod;
results.system.mdot_evap_still = mdot_evap_still;
results.system.SEC_blower      = W_blower    / max(mdot_prod,eps) / 3600;  % [kWh/m3]
results.system.SEC_heater      = W_heater_el / max(mdot_prod,eps) / 3600;  % [kWh/m3]
results.system.SEC_compressor  = W_awg_el    / max(mdot_prod,eps) / 3600;  % [kWh/m3]
results.system.SEC_air_heater  = W_airht_el  / max(mdot_prod,eps) / 3600;  % [kWh/m3]
results.system.SEC_external    = W_external  / max(mdot_prod,eps) / 3600;  % [kWh/m3]
a_COP_for_band = Q_evap_coil / max(results.downstream.W_awg, eps);   % [-] central COP
results.system.Q_solar_absorbed = Q_solar_absorbed;
results.system.basis            = results.mfw_daily_basis;

% Sensitivity of the AWG term to the eta_II assumption (0.35-0.55).
% The heater and the air heater carry no such band: both are resistance
% loads, ~1.0 efficient by construction, with duties set by the mass and
% energy balances respectively.
%
% W_air_heater MUST BE INCLUDED. Omitting it does not merely shift the
% band -- it puts the band entirely BELOW the central SEC it is supposed
% to bracket, because the air heater is the single largest load in this
% configuration, and a band that does not contain its own central value
% is not an uncertainty statement at all.
% The coil duty must be on the SAME basis as W_awg_el above, or the band
% is computed from a final-instant duty while the central value uses a
% window mean. Recovering it from W_awg_el keeps the two locked together
% whichever branch was taken.
COP_central = a_COP_for_band;
Q_coil_basis = W_awg_el * max(COP_central, eps);
for ii = 1:2
    eta_band = [0.35 0.55];
    COP_b    = eta_band(ii) * P.T_coil / max(T_cond_awg - P.T_coil, 1);
    W_b      = W_blower + W_heater_el + W_airht_el ...
               + Q_coil_basis / max(COP_b, eps);
    results.system.SEC_external_band(ii) = W_b / max(mdot_prod,eps) / 3600;
end

% SELF-CHECK: the band must straddle the central value. If it does not,
% a load is missing from W_b.
if results.system.SEC_external < min(results.system.SEC_external_band) - 1e-9 || ...
   results.system.SEC_external > max(results.system.SEC_external_band) + 1e-9
    warning('ZLDD:SECBand', ...
        ['SEC band [%.2f, %.2f] does not contain the central SEC %.2f ' ...
         'kWh/m3. A load is missing from the band expression.'], ...
         min(results.system.SEC_external_band), ...
         max(results.system.SEC_external_band), results.system.SEC_external);
end

results.steady_check.Q_storage_total    = Q_storage_total;
results.steady_check.Q_balance_residual = Q_balance_residual;
results.steady_check.Q_balance_rel      = Q_balance_residual / max(Q_in_total, eps);
% Q_out_total = -(- Q_reflect_glass - Q_glass_conv_loss - Q_glass_rad_loss- Q_wall_ss - Q_feed_out_ss - Q_vapor_out_ss ...
%              - Q_ground_loss_end - Q_wallstrip_loss_end);

results.steady_check.Q_out_total       =Q_out_total ;
results.steady_check.Q_ground_loss    = Q_ground_loss_end;
results.steady_check.Q_wallstrip_loss = Q_wallstrip_loss_end;
results.steady_check.Q_in_total       = Q_in_total;

results.steady_check.Q_solar          = Q_solar;
results.steady_check.Q_solar_aperture = Q_solar_aperture;
results.steady_check.Q_solar_absorbed = Q_solar_absorbed;
results.steady_check.Q_solar_unused   = Q_solar_unused;

% ---- Thermal input by source ----
% The relative contributions of absorbed solar radiation, intake-air
% preheat and feed preheat determine which energy stream governs
% performance and hence which efficiency definition is appropriate.
Q_src_total = max(Q_solar_absorbed,0) + max(Q_fan_thermal,0) + max(Q_feed_thermal,0);
results.steady_check.solar_share_pct = 100*max(Q_solar_absorbed,0) / max(Q_src_total,eps);
results.steady_check.fan_share_pct   = 100*max(Q_fan_thermal,0)   / max(Q_src_total,eps);
results.steady_check.feed_share_pct  = 100*max(Q_feed_thermal,0)  / max(Q_src_total,eps);

% ---- Gained output ratio on solar-only and total thermal bases ----
% A solar-only denominator is the convention for passive solar stills;
% a total thermal input denominator is the convention for
% humidification-dehumidification systems. Both are reported, since the
% appropriate choice depends on which stream supplies the energy.
results.steady_check.GOR_solar_only  = Q_evap / max(Q_solar_absorbed, eps);
results.steady_check.GOR_total_input = Q_evap / max(Q_src_total, eps);
results.steady_check.Q_fan_thermal = Q_fan_thermal;
results.steady_check.Q_feed_thermal = Q_feed_thermal;

results.steady_check.Q_evap = Q_evap;
results.steady_check.mdot_evap_total = mdot_evap_total;
% ---- Thermal efficiency ----
% The denominator must be the total energy ENTERING the still, which is
% the deliberately supplied energy plus any boundary streams that have
% reversed direction. When the evaporative load drives the cascade below
% ambient temperature the exhaust, glazing convection and wall
% conduction all reverse and carry heat inward; that heat is available
% for evaporation and belongs in the denominator.
Q_entering_total = Q_in_total + Q_gain_total;
results.steady_check.Q_entering_total   = Q_entering_total;
results.steady_check.eta_thermal_steady = Q_evap / max(Q_entering_total, eps) * 100;
results.steady_check.Q_wall_loss  = Q_wall_end;
results.steady_check.Q_feed_out   = Q_feed_out_end;
results.steady_check.Q_vapor_out  = Q_vapor_out_end;
results.steady_check.Q_gain_total = Q_gain_total;
results.steady_check.Q_reflect_glass   = Q_reflect_glass;
results.steady_check.Q_glass_conv_loss = Q_glass_conv_loss;
results.steady_check.Q_glass_rad_loss  = Q_glass_rad_loss;
results.steady_check.Q_in_total        = Q_in_total;




end % extract_results




function D = preheater_heater_awg(brine, Tv_exit, wv_exit, wv_coil, P)
% DOWNSTREAM UNIT BALANCES for the CLOSED-LOOP configuration.
%
%   brine    - results.brine_out struct. The plate-Ns outlet at the
%              FINAL INSTANT: the cascade outlet unaltered, since no
%              unit stands between the last plate and this heater.
%              Do NOT pass the _win (window-mean) fields here -- mdot
%              and w_salt averaged independently do not describe any
%              stream that physically existed, and feeding them under-
%              sizes the crystalliser. See the basis note at brine_out.
%   Tv_exit  - cascade exhaust air temperature [K]
%   wv_exit  - cascade exhaust humidity ratio  [kg/kg]
%   wv_coil  - humidity ratio at the AWG coil, saturated at P.T_coil
%
% Evaluated once, after the solve. Feed-forward: nothing here returns to
% the cascade.

D = struct();
P.Patm   = 101325;
T_sat  = P.bl.T_sat_vap;               % [K] CONDENSING temperature of the vapour (no BPE)
T_boil = P.bl.T_boil;                  % [K] boiling liquor  = T_sat + BPE


% ==================================================================
% 1. ELECTRIC HEATER
% ------------------------------------------------------------------
% Concentrated brine from the floor stage is evaporated until the salt
% mass fraction reaches P.bl.w_salt_target. Salt is conserved, so the
% outlet brine flow follows directly from the salt flow, and the vapour
% raised is the difference.
% ==================================================================
% BASIS GUARD. This inlet must BE the plate-Ns outlet, not a
% reconstruction of it. w_salt is taken as C/rho of that same instant,
% so mdot*w_salt reproduces the audited outlet salt flux; averaging
% mdot and w_salt over a window independently does not.
if isfield(brine,'basis') && ~contains(brine.basis,'final')
    warning('heater:basis', ...
        ['Heater fed a %s brine stream. Expected the final-instant ' ...
         'plate-Ns outlet; sizing will not match PART D.'], brine.basis);
end

mdot_brine_in = max(brine.mdot, 0);                    % [kg/s] TOTAL brine (water + salt)
w_salt_in     = min(max(brine.w_salt, 0), 0.999);      % [-] inlet salt mass fraction
mdot_salt     = mdot_brine_in * w_salt_in;             % [kg/s] conserved
mdot_water_in = mdot_brine_in - mdot_salt;             % [kg/s] water entering the heater
mdot_brine_out= mdot_salt / P.bl.w_salt_target;        % [kg/s] to the salt product
mdot_vap      = max(mdot_brine_in - mdot_brine_out, 0);% [kg/s] vapour raised

% Water closure across the heater: everything that is not salt either
% leaves as vapour or stays as the residual water in the slurry.
D.mdot_water_in  = mdot_water_in;
D.water_residual = mdot_water_in - mdot_vap - (mdot_brine_out - mdot_salt);

% Heater duty: sensible heating of the whole inlet stream from its
% arrival temperature to the vapour temperature, plus latent heat for
% the fraction evaporated. Latent heat is taken at the INLET salinity;
% the salinity correction at the target fraction is outside the
% validity of latent_heat() and is not extrapolated.
%
% THREE ADDITIONAL TERMS relative to a plain latent+sensible duty, all
% consequences of continuous operation at NaCl saturation:
%
%  (a) The brine is raised to T_boil, the LIQUOR temperature, not to the
%      element temperature. Using T_elem here would charge the duty for
%      heating the whole stream to the jacket, which never happens.
%  (b) BPE SUPERHEAT. The steam leaves at T_boil and must be raised the
%      extra BPE degrees above its own saturation temperature. Small
%      (~0.6% of latent) but it is a real parasitic loss: it is spent
%      here and, because the vapour desuperheats in the first
%      centimetres of the feed preheat HX, it returns nothing.
%  (c) CRYSTALLISATION CREDIT. Precipitating NaCl releases enthalpy.
%      This is a genuine (small) reduction in duty and is included for
%      honesty in the energy accounting, since the plant is claimed as
%      zero-liquid-discharge and the salt product is real.
Cp_brine = water_cp(brine.T, brine.C, P);              % [J/kg-K]
hfg_vap  = latent_heat(T_sat, brine.C, P);             % [J/kg]
Cp_v_bpe = 1860 + 0.12*(0.5*(T_boil + T_sat) - 273.15);% [J/kg-K]

Q_sens   = mdot_brine_in * Cp_brine * max(T_boil - brine.T, 0);      % [W] to the LIQUOR temperature
Q_lat    = mdot_vap * hfg_vap;                                       % [W]
Q_bpe    = mdot_vap * Cp_v_bpe * P.bl.BPE;                           % [W] BPE superheat penalty
Q_cryst  = mdot_salt * P.bl.h_cryst;                                 % [W] released -> credit

Q_heater_useful = Q_sens + Q_lat + Q_bpe - Q_cryst;                  % [W]
W_heater = Q_heater_useful * (1 + P.bl.f_heatloss) / P.bl.eta_heater;% [W] electrical

D.Q_heater_bpe   = Q_bpe;
D.Q_heater_cryst = Q_cryst;
D.T_boil         = T_boil;
D.T_sat_vap      = T_sat;
D.BPE            = P.bl.BPE;

D.mdot_brine_in    = mdot_brine_in;
D.w_salt_in        = w_salt_in;
D.mdot_salt        = mdot_salt;
D.mdot_brine_out   = mdot_brine_out;
D.w_salt_out       = P.bl.w_salt_target;
D.mdot_vap         = mdot_vap;
D.Q_heater_sensible= Q_sens;
D.Q_heater_latent  = Q_lat;
D.Q_heater_useful  = Q_heater_useful;
D.W_heater         = W_heater;
D.T_vap_heater     = P.bl.T_vap_heater;

% ==================================================================
% 2. FEED PREHEAT HX  (evaporator vapour -> SEAWATER FEED)
% ------------------------------------------------------------------
% The exchanger sits on the FEED, not on the air. Evaporator vapour
% desuperheats by BPE, then condenses ISOTHERMALLY at T_sat, warming
% incoming seawater on its way to the cascade.
%
% WHICH VARIABLE IS SPECIFIED. Neither the outlet temperature nor a
% condensed fraction: BOTH are outputs. The feed temperature is whatever
% the vapour can support,
%
%     T_feed_out = T_feed_in + Q_absorbed / (mdot_feed*Cp_f)
%
% capped by the pinch at T_sat - dT_pinch. The condensed fraction then
% follows from the duty the feed actually absorbed. The condensed
% fraction is therefore not a free design choice: with the exchanger on
% the feed side, the feed's heat capacity rate decides how much latent
% duty it can take.
%
% Two regimes:
%   pinch NOT binding : the feed absorbs everything the vapour offers,
%                       f_cond -> 1, nothing passes to the AWG.
%   pinch binding     : the feed saturates at T_sat - dT_pinch and the
%                       surplus vapour passes on to the AWG, where it
%                       condenses anyway. Not a loss -- a rerouting.
%
% The feed outlet temperature is TAKEN from feed_preheat_T() via the
% post-solve reconciliation of P.Tfeed -- the same function, the same
% instant, the same value the RHS used. It is deliberately NOT
% recomputed here; see the note at the assignment below.
% ==================================================================
Cp_f_hx      = water_cp(P.Tfeed_in, P.TDSfeed, P);     % [J/kg-K]
Cden_feed    = P.mdot_feed_total * Cp_f_hx;            % [W/K] feed capacity rate

Cp_vap_super = 1860 + 0.12*(0.5*(T_boil + T_sat) - 273.15);
Q_desuper    = mdot_vap * Cp_vap_super * max(T_boil - T_sat, 0);   % [W] BPE superheat only
hfg_sat      = latent_heat(T_sat, 0, P);               % [J/kg] pure-water vapour

Q_cond_capacity = mdot_vap * hfg_sat;                  % [W] available latent duty
Q_vap_total     = Q_desuper + Q_cond_capacity;         % [W] everything the vapour offers

% ---- Feed outlet: TAKEN FROM THE SOLVE, NOT RECOMPUTED ----
% SINGLE POINT OF TRUTH. P.Tfeed was reconciled to feed_preheat_T()'s
% value at the final state before this function was called, so it IS
% the temperature the cascade received. Re-deriving it here would give
% a DIFFERENT number, because this block's mdot_vap comes from the
% reported brine (rho*u*delta at the outlet node) while the RHS used
% the mass-balance brine. Those two estimators differ by the film
% storage rate, so recomputing would put one heater on two flows.
%
% It also keeps the two reports from disagreeing about which ceiling
% binds -- the startup cap or the pinch -- for one exchanger at one
% instant.
T_cap        = T_sat - P.dT_pinch_HX;                  % [K] hard ceiling
dT_unpinched = Q_vap_total / max(Cden_feed, eps);      % [K] if all duty were absorbed
T_feed_out   = P.Tfeed_cascade_in;                     % [K] <- SOLVED cascade inlet (hot side)
pinched      = (P.Tfeed_in + dT_unpinched) > T_cap;

Q_feed_absorbed = Cden_feed * max(T_feed_out - P.Tfeed_in, 0);      % [W]

% Consistency guard: a solved feed temperature above the pinch ceiling is
% thermodynamically inadmissible for this exchanger, and the sizing below
% would then be quoted against an impossible stream.
if T_feed_out > T_cap + 1e-6
    warning('ZLDD:FeedAbovePinch', ...
        ['Solved feed outlet %.2f K exceeds the pinch ceiling %.2f K ' ...
         '(T_sat %.2f K less dT_pinch %.2f K). The exchanger sizing ' ...
         'below is not physical.'], T_feed_out, T_cap, T_sat, P.dT_pinch_HX);
end

% ---- Condensed fraction: an OUTPUT of the duty the feed took ----
% Desuperheating is served first; whatever remains condenses vapour.
if Q_cond_capacity > eps
    f_cond_raw = (Q_feed_absorbed - Q_desuper) / Q_cond_capacity;
else
    f_cond_raw = 0;                                    % no vapour at all
end
f_cond   = min(max(f_cond_raw, 0), 1);
feasible = true;    % cannot be infeasible: the pinch caps the demand

mdot_cond_preheater  = f_cond * mdot_vap;              % [kg/s]
Q_preheater_supplied = min(Q_desuper, Q_feed_absorbed) + mdot_cond_preheater * hfg_sat;
Q_preheater_imbalance = Q_feed_absorbed - Q_preheater_supplied;     % [W] signed
Q_preheater_deficit   = Q_preheater_imbalance;

D.Q_feed_absorbed       = Q_feed_absorbed;
D.Q_air_required        = Q_feed_absorbed;   % alias name for the HX duty
D.Q_vap_total           = Q_vap_total;
D.Q_desuperheat         = Q_desuper;
D.Q_cond_capacity       = Q_cond_capacity;
D.f_cond_raw            = f_cond_raw;
D.f_cond                = f_cond;
D.mdot_cond_preheater   = mdot_cond_preheater;
D.Q_preheater_supplied  = Q_preheater_supplied;
D.Q_preheater_deficit   = Q_preheater_deficit;
D.preheater_feasible    = feasible;
D.Q_preheater_imbalance = Q_preheater_imbalance;   % [W] signed
D.T_feed_in             = P.Tfeed_in;
D.T_feed_out            = T_feed_out;
D.T_feed_cap            = T_cap;
D.feed_pinched          = pinched;
D.dT_feed               = T_feed_out - P.Tfeed_in;
D.Cden_feed             = Cden_feed;

% ==================================================================
% 2b. VAPOUR-FIRED AIR PREHEATER  +  ELECTRIC AIR HEATER
% ------------------------------------------------------------------
% Two units in series on the air line between the AWG and the still.
%
% Surplus evaporator vapour is 373 K saturated steam -- the highest-grade
% stream in the plant. Condensing it against the 295 K AWG coil would be
% a waste of exergy, so it preheats the loop air FIRST and the resistance
% element trims only the shortfall. Its condensate is PRODUCT and must be
% counted; see the product-assembly block in extract_results().
%
% SPECIFICATION, and the two units differ:
%   preheater : DUTY-specified. Takes whatever surplus exists, so T_int
%               is the OUTPUT. It cannot be outlet-T-specified, because
%               the surplus falls to ZERO in the supply-limited regime
%               and the demanded duty would not exist.
%   heater    : OUTLET-T-specified (P.T_air_in). Duty is the OUTPUT. It
%               absorbs the variability, which is what lets the pair hold
%               the setpoint against a swinging steam supply.
%
% SIZING. The element is NOT smaller for having the preheater in front of
% it. The steam supply falls to ZERO whenever the feed HX is supply-
% limited, so the element must be able to carry the whole lift alone.
% Rate it at Q_air_full: it runs at reduced duty, not reduced rating.
% D.f_air_recovered reports what fraction the steam covered on THIS run;
% it is a result and varies with the air flow and the cascade duty.
% ==================================================================
% ==================================================================
% 2a-bis. AIR-TO-AIR RECUPERATOR  (still exhaust <-> AWG outlet)
% ------------------------------------------------------------------
% PASSIVE. No power, no moving parts, no control. The loop circuit is
% still -> AWG -> heater -> still, and the duct passes through one box
% TWICE, separated by a wall: once carrying warm humid air out of the
% still, once carrying cold dry air back from the coil.
% The two never mix; heat crosses the wall from the outgoing stream into
% the returning one.
%
% WHY IT EXISTS. Without it the plant cools the loop air across a span
% and then heats it back across the same span, paying a compressor for
% the first and a resistance element for the second. That is a pure
% irreversibility: the recuperator does not work around it, it removes
% it. Both duties shrink SIMULTANEOUSLY -- the coil sees pre-cooled air,
% the heater sees pre-warmed air.
%
% WHAT IT CANNOT DO. It moves SENSIBLE heat only. The latent duty stays
% entirely with the coil, because only the coil can take the air below
% its dew point. The recoverable ceiling is (Tv_exit - T_coil): the
% exchanger can never drive the hot stream below the cold stream's
% inlet, whatever its area. A narrow spread therefore caps the recovery
% regardless of effectiveness.
%
% CONDENSATION ON THE HOT SIDE. The outgoing stream may cross its dew
% point inside the box, in which case water condenses there and the unit
% needs a drain. That condensate is PRODUCT and is accounted below; the
% latent release is also credited to the hot side, so the cold-side rise
% is not understated. This is checked rather than assumed, because
% whether it happens depends on the operating point.
% ==================================================================
% ---- WHICH SIDE CARRIES Cmin --------------------------------------
% The hot side CONDENSES, so its effective capacity rate is enormous: it
% gives up latent heat at almost constant temperature. The COLD side (dry
% loop air, sensible only) is therefore Cmin, and effectiveness must be
% defined on it:
%
%     eps = (T_cold_out - T_cold_in) / (T_hot_in - T_cold_in)
%
% Defining eps on the hot side instead and then handing the cold side the
% full sensible+latent release is NOT equivalent: it credits the cold
% stream with duty the exchanger's area cannot actually transfer, and the
% resulting outlet can exceed the hot inlet, which no counterflow unit can
% do. A clamp at T_hot_in would hide that rather than prevent it.
Cp_recup_h = moist_air_cp(Tv_exit, wv_exit);                     % [J/kg-K] hot side
Cden_recup = P.mdot_da * Cp_recup_h;                             % [W/K] hot, sensible only

Cp_da_ah = 1005 + 0.05*(0.5*(P.T_air_in + P.T_coil) - 273.15);
Cp_v_ah  = 1860 + 0.12*(0.5*(P.T_air_in + P.T_coil) - 273.15);
Cden_air = P.mdot_da*Cp_da_ah + P.mdot_da*P.wv_fan_in*Cp_v_ah;   % [W/K] cold = Cmin

dT_recup_max    = max(Tv_exit - P.T_coil, 0);                    % [K] ceiling
T_air_recup_out = P.T_coil + P.eps_recup * dT_recup_max;         % [K] cold outlet
dT_air_recup    = T_air_recup_out - P.T_coil;                    % [K]

% Duty is then fixed by the cold side. The exchanger is adiabatic, so the
% hot side must give up EXACTLY this much -- no more, no less.
Q_recup_total = Cden_air * dT_air_recup;                         % [W]

% ---- HOT-SIDE OUTLET: SOLVED, NOT ASSUMED -------------------------
% The hot stream sheds Q_recup_total as sensible cooling PLUS whatever it
% condenses on the way. Condensation depends on the outlet temperature it
% is trying to find, so this is implicit:
%
%   Q = mdot_da*Cp*(Tv_exit - T_out) + mdot_da*(wv_exit - wsat(T_out))*hfg
%
% Monotonic in T_out, so bisection is robust and needs no derivative.
% Bracketed by the cold inlet (the hot stream cannot be driven below it)
% and its own inlet.
lo = P.T_coil;  hi = Tv_exit;
for it_recup = 1:60
    T_try   = 0.5*(lo + hi);
    w_sat_t = humidity_ratio_from_RH(1.0, T_try, P);
    m_c_try = max(P.mdot_da*(wv_exit - w_sat_t), 0);             % condensed by T_try
    Q_try   = Cden_recup*max(Tv_exit - T_try, 0) ...
              + m_c_try*latent_heat(T_try, 0, P);
    if Q_try > Q_recup_total, lo = T_try; else, hi = T_try; end
end
T_recup_hot_out = 0.5*(lo + hi);                                 % [K] -> AWG coil

wv_recup_sat    = humidity_ratio_from_RH(1.0, T_recup_hot_out, P);
mdot_cond_recup = max(P.mdot_da*(wv_exit - wv_recup_sat), 0);    % [kg/s] drained here
hfg_recup       = latent_heat(T_recup_hot_out, 0, P);            % [J/kg]
Q_recup_latent  = mdot_cond_recup * hfg_recup;                   % [W]
Q_recup_sensible= Cden_recup * max(Tv_exit - T_recup_hot_out, 0);% [W]
dT_recup        = Tv_exit - T_recup_hot_out;                     % [K]

% Adiabatic closure: what the hot side gave up IS what the cold side took.
resid_recup = (Q_recup_sensible + Q_recup_latent) - Q_recup_total;
if abs(resid_recup) > 1e-6*max(Q_recup_total, eps)
    warning('ZLDD:RecupEnergy', ...
        'Recuperator energy does not close: %.3e W of %.3e W.', ...
        resid_recup, Q_recup_total);
end

% Humidity entering the coil is whatever survived the recuperator.
wv_to_coil = wv_exit - mdot_cond_recup / max(P.mdot_da, eps);    % [kg/kg]

% The lift the downstream units must still supply, AFTER recovery.
Q_air_full  = Cden_air * max(P.T_air_in - T_air_recup_out, 0);   % [W] remaining
Q_air_gross = Cden_air * max(P.T_air_in - P.T_coil, 0);          % [W] with no recuperator

D.eps_recup            = P.eps_recup;
D.T_recup_hot_out      = T_recup_hot_out;
D.T_air_recup_out      = T_air_recup_out;
D.dT_recup             = dT_recup;
D.dT_recup_max         = dT_recup_max;
D.dT_air_recup         = dT_air_recup;
D.Q_recup_sensible     = Q_recup_sensible;
D.Q_recup_latent       = Q_recup_latent;
D.Q_recup_total        = Q_recup_total;
D.mdot_cond_recup      = mdot_cond_recup;
D.wv_to_coil           = wv_to_coil;
D.Q_air_gross          = Q_air_gross;

% Condensate subcooling coefficient. It is defined here because this unit
% needs it first; section 3 aliases this value rather than recomputing
% it, so the two cannot drift apart.
Cp_cond_est = water_cp(0.5*(T_sat + P.T_coil), 0, P);            % [J/kg-K]

% ---- SUBCOOLING DATUM: THE COLD STREAM THIS UNIT ACTUALLY SEES -------
% The coil is not part of this exchanger. The cold stream here is loop
% air leaving the RECUPERATOR at T_air_recup_out, and in counterflow the
% condensate cannot leave colder than the cold-stream INLET. Releasing
% hfg + Cp*(T_sat - T_coil) would credit the unit with
% (T_air_recup_out - T_coil) degrees of subcooling that no surface in it
% can perform, overstating the steam-fired duty and understating the
% electric trim by the same amount -- the trim is formed as a residual
% (Q_air_full - Q_air_cond) below, so such a bias would propagate
% straight into the reported SEC.
%
% The AWG in section 3 uses the T_coil datum: that unit does reject to
% the coil, so its subcooling term belongs on that datum.
T_sub_datum = T_air_recup_out;                                   % [K] cold-side inlet
h_release   = hfg_sat + Cp_cond_est*max(T_sat - T_sub_datum, 0); % [J/kg]

% Surplus steam reaching this unit -- the vapour the FEED did not take.
mdot_vap_surplus = max(mdot_vap - mdot_cond_preheater, 0);       % [kg/s]
Q_surplus_avail  = mdot_vap_surplus * h_release;                 % [W]

D.h_release_airheat = h_release;
D.T_sub_datum       = T_sub_datum;
D.mdot_vap_surplus  = mdot_vap_surplus;
D.Q_surplus_avail   = Q_surplus_avail;

% Approach limit: air cannot be driven closer to the steam than this.
% Inactive at the current Qfan, present so the unit cannot misbehave if
% the air flow is reduced.
T_int_cap  = T_sat - P.dT_pinch_cond;                            % [K]

% ---- DATUM: THE PREHEATER SITS DOWNSTREAM OF THE RECUPERATOR --------
% This unit does not take air straight off the coil: the recuperator is
% upstream and has already lifted it to T_air_recup_out, so both the
% approach cap and T_int are referred to that state rather than to
% P.T_coil. Referring them to T_coil would double-count the recuperated
% span in the cap and would place T_int below the temperature the air has
% already reached.
Q_cond_cap = Cden_air * max(min(P.T_air_in, T_int_cap) - T_air_recup_out, 0);

Q_air_cond = min([Q_surplus_avail, Q_air_full, Q_cond_cap]);     % [W] recovered
T_int      = T_air_recup_out + Q_air_cond/max(Cden_air, eps);    % [K] OUTPUT

Q_air_trim   = max(Q_air_full - Q_air_cond, 0);                  % [W]
W_air_heater = Q_air_trim / P.eta_air_heater;                    % [W] ELECTRICAL

% Vapour actually consumed here. The remainder still goes to the AWG.
mdot_cond_airheat = Q_air_cond / max(h_release, eps);            % [kg/s]

D.Q_air_full         = Q_air_full;
D.Q_air_cond         = Q_air_cond;
D.Q_air_trim         = Q_air_trim;
D.Q_air_heater       = Q_air_full;      % alias name for the TOTAL lift
D.T_int              = T_int;
D.W_air_heater       = W_air_heater;
D.mdot_cond_airheat  = mdot_cond_airheat;
D.f_air_recovered    = Q_air_cond / max(Q_air_full, eps);
D.T_air_in           = P.T_air_in;
% Spans on the datum each unit actually works across. dT_air_cond is the
% rise the STEAM delivers (from the recuperator outlet, not the coil);
% dT_air_gross is the coil-to-setpoint span the loop needs in total, of
% which the passive recuperator supplies dT_air_recup for free.
D.dT_air_cond        = T_int - T_air_recup_out;
D.dT_air_trim        = P.T_air_in - T_int;
D.dT_air_heater      = P.T_air_in - T_air_recup_out;   % lift left after the recuperator
D.dT_air_gross       = P.T_air_in - P.T_coil;          % coil -> still inlet, all sources
% ==================================================================
% 3. AWG  (two inlet streams)
% ------------------------------------------------------------------
% Stream A -- the vapour NOT condensed in the preheater. Essentially
%             pure water vapour, so it condenses completely; the
%             condensate is then subcooled to the coil temperature.
% Stream B -- the humid air leaving the cascade. Cooled to T_coil and
%             leaves SATURATED at that temperature, so recovery is set
%             by the humidity-ratio difference. Guarded so that an
%             exhaust drier than the coil saturation state gives zero
%             recovery rather than negative condensate.
% ==================================================================
% THREE consumers of the raised vapour: the feed HX, the vapour-fired air
% preheater (section 2b), and finally the AWG. Omitting the air preheater
% here would double-count its condensate.
mdot_vap_to_awg = max(mdot_vap - mdot_cond_preheater - mdot_cond_airheat, 0);  % [kg/s] stream A
hfg_coil        = latent_heat(P.T_coil, 0, P);         % [J/kg]
Cp_cond         = Cp_cond_est;                         % hoisted in 2b; do NOT redefine
Q_awg_vapour    = mdot_vap_to_awg * ( hfg_coil + Cp_cond*max(T_sat - P.T_coil,0) );  % [W]

% The vapour split must close: feed HX + air preheater + AWG == raised.
% Asserted rather than trusted, because the three terms are computed in
% three different blocks and a silent leak here shows up as missing
% product rather than as an error.
resid_split = mdot_vap - (mdot_cond_preheater + mdot_cond_airheat + mdot_vap_to_awg);
if abs(resid_split) > 1e-9*max(mdot_vap, eps)
    warning('ZLDD:VapourSplit', ...
        'Vapour split does not close: residual %.3e kg/s of %.3e kg/s raised.', ...
        resid_split, mdot_vap);
end

% Stream B: dew-point guard. No recovery unless the coil is genuinely
% below the incoming dew point.
%
% THE COIL SEES THE RECUPERATOR OUTLET, NOT THE CASCADE EXHAUST. The air
% arrives pre-cooled and, if it crossed its dew point in the recuperator,
% already partly dried. Charging this block against the raw exhaust state
% would bill the compressor for sensible cooling the passive exchanger
% already did, and would double-count any moisture the recuperator drained.
% Total condensate over the two units is unchanged; only its SPLIT moves.
mdot_cond_air = max( P.mdot_da * (wv_to_coil - wv_coil), 0 );       % [kg/s]
Cp_ret        = moist_air_cp(T_recup_hot_out, wv_to_coil);
Q_awg_air     = P.mdot_da*(1 + wv_to_coil) * Cp_ret ...
                * max(T_recup_hot_out - P.T_coil, 0) ...
                + mdot_cond_air * hfg_coil;                          % [W]

% Air-side moisture must close across BOTH units: what the cascade
% evaporated leaves either in the recuperator drain or on the coil.
resid_air = (mdot_cond_recup + mdot_cond_air) - P.mdot_da*(wv_exit - wv_coil);
if abs(resid_air) > 1e-9*max(P.mdot_da*wv_exit, eps)
    warning('ZLDD:RecupMoisture', ...
        'Recuperator/coil moisture split does not close: %.3e kg/s.', resid_air);
end

Q_awg_total = Q_awg_vapour + Q_awg_air;                              % [W]

% Heat-pump lift. The reheat that returns loop air to the still inlet is
% supplied by the SEPARATE electric air heater, not by this unit's
% condenser, so the AWG still rejects to AMBIENT and the lift remains
% (Ta + approach) - T_coil.
%
% The DUTY is set by the loop, not by a single-pass throughput: the coil
% processes recirculated air on every pass, so Q_awg_air is charged
% continuously against the full loop flow.
T_cond_awg = P.Ta + P.dT_cond_approach;                              % [K]
COP        = P.bl.eta_II_awg * P.T_coil / max(T_cond_awg - P.T_coil, 1);
W_awg      = Q_awg_total / max(COP, eps);                            % [W]

D.mdot_cond_awg_vapour = mdot_vap_to_awg;
D.mdot_cond_air        = mdot_cond_air;
D.T_coil_inlet         = T_recup_hot_out;   % what the coil actually sees
D.Q_awg_vapour         = Q_awg_vapour;
D.Q_awg_air            = Q_awg_air;
D.Q_awg_total          = Q_awg_total;
D.T_cond_awg           = T_cond_awg;
D.COP                  = COP;
D.W_awg                = W_awg;
D.Q_awg_reject         = Q_awg_total + W_awg;
D.wv_coil              = wv_coil;
D.wv_exit              = wv_exit;

% ---- Combined electrical demand of the downstream chain ----
% THREE electrical loads: the brine evaporator, the AWG compressor, and
% the electric air heater, which is a separate unit in this flowsheet.
D.W_downstream = W_heater + W_awg + W_air_heater;                    % [W]

% ---- FEED PREHEAT CEILING --------------------------------------------
% The largest FEED temperature rise the evaporator vapour can support,
% and the pinch ceiling that sits above it. Two distinct limits:
%   dT_preheat_max : supply-limited -- all the vapour duty absorbed
%   dT_pinch_max   : approach-limited -- the feed cannot get closer to
%                    the condensing steam than dT_pinch_HX
% Whichever is smaller binds. When the PINCH binds, surplus vapour is
% rerouted to the AWG rather than lost, so this is not an infeasibility.
D.dT_preheat_max    = Q_vap_total / max(Cden_feed, eps);             % [K] supply limit
D.dT_pinch_max      = T_cap - P.Tfeed_in;                            % [K] approach limit
D.dT_preheat_actual = T_feed_out - P.Tfeed_in;                       % [K]
D.dT_preheat_slack  = min(D.dT_preheat_max, D.dT_pinch_max) - D.dT_preheat_actual;
% Three regimes, not two. The DESIGN target can bind before either
% physical limit does, and the label must be able to say so.
if isfield(P,'Tfeed_regime')
    D.binding_limit = P.Tfeed_regime;        % 'supply' | 'target' | 'pinch'
else
    D.binding_limit = ternary(D.dT_pinch_max < D.dT_preheat_max, ...
                              'pinch', 'vapour supply');
end
end % preheater_heater_awg

function [Q_fan_electrical, Q_fan_thermal] = fan_energy(P,~)
% Q_fan_electrical: blower work to push air through the duct [J]
% Q_fan_thermal:    net EXTERNAL sensible heat carried in by the intake
%                    air relative to ambient (zero if Tfan_in == Ta;
%                    the Tfan_in -> Tv(1) rise is NOT counted here since
%                    that heat is internally sourced from solar, already
%                    tracked via q_conv_wv/q_conv_vp in the plate balance)

% ---- Electrical: duct friction + U-turn (bend) losses -> blower power ----
rho_a  = air_density(P.Tfan_in);
mu_a   = air_viscosity(P.Tfan_in);
u_duct = P.Qfan ./ P.Acomp;
Re     = rho_a .* u_duct .* P.Dh ./ mu_a;

f = zeros(length(Re),1);
for i=1:length(Re)
    if Re(i) < 2300
        f(i) = 64./Re(i);
    else              
        f(i) = 0.184 * Re(i).^(-0.2); 
    end
end

dP_friction = f .* (P.L_stage ./ P.Dh) .* 0.5 .* rho_a .* u_duct.^2;

% -- 180-deg U-turn loss at each of the (Ns-1) gap-to-gap transitions --
K_bend      = 1.8;   % sharp-edged close-return bend, rectangular duct -- tune to your geometry
dP_bends    = K_bend * 0.5 * rho_a * u_duct(1:end-1).^2;   % Ns-1 turns

% -- optional: entrance (ambient -> gap Ns) and exit (gap 1 -> ambient) losses --
K_entrance  = 0.5;
K_exit      = 1.0;
dP_ends     = (K_entrance + K_exit) * 0.5 * rho_a * mean(u_duct.^2);

dP_total = sum(dP_friction) + sum(dP_bends) + dP_ends;

eta_fan = 0.6;
Q_fan_electrical = (dP_total * P.Qfan / eta_fan) * 3600 * P.t_operating;

% ---- Thermal: the ELECTRIC AIR HEATER's duty --------------------------
% DATUM: T_coil. The loop is closed and no ambient air enters, so the air
% arriving at the still is not drawn from outside -- it comes from the
% AWG coil. Referencing Ta here would measure a stream that does not
% exist and would charge (or credit) the still with joules no unit
% supplies. The quantity formed on the T_coil datum is the electric air
% heater's duty, which is a genuine EXTERNAL input to the loop.
mdot_da   = P.mdot_da;
% Evaluated per component (dry air and water vapour separately) so that
% this integrated quantity is consistent with the instantaneous stream
% used in the steady power balance. Using the dry-air specific heat
% alone would omit the sensible enthalpy carried by the vapour fraction,
% and the two reports would disagree on the same stream.
T_mean_ah  = 0.5*(P.T_air_in + P.T_coil);
Cp_da_fan  = 1005 + 0.05*(T_mean_ah - 273.15);
Cp_vap_fan = 1860 + 0.12*(T_mean_ah - 273.15);
mdot_vap_fan = P.mdot_da * P.wv_fan_in;

Q_fan_thermal = ( mdot_da      * Cp_da_fan ...
                + mdot_vap_fan * Cp_vap_fan ) ...
                * (P.T_air_in - P.T_coil) * 3600 * P.t_operating;   % [J], signed
end

function Q_feed = feed_thermal_energy(P)
% Rigorous sensible-heat input from feed stream, integrating Cp(T) at
% FIXED feed salinity (concentration hasn't started yet at the inlet
% boundary -- that happens progressively down the cascade, handled
% separately by the local Cp_w term inside d_Tw).

% DATUM: the SUPPLIED seawater temperature, P.Tfeed_in, not ambient.
% The feed preheat HX raises the seawater from Tfeed_in to P.Tfeed
% (solved), and that rise is supplied INTERNALLY by evaporator vapour.
% Integrating from Ta instead would mix an ambient datum into a term
% whose two ends are both plant streams.
mdot_feed = P.Vfeed * P.rhow_in / 1000 / (P.t_operating * 3600);   % [kg/s]

% Numerically integrate Cp(T, Cfeed) from the supplied seawater
% temperature to the preheated cascade inlet temperature.
T_range = linspace(P.Tfeed_in, P.Tfeed_cascade_in, 50);   % cold side -> hot side
Cp_vals = water_cp(T_range, P.TDSfeed, P);   % Cp at each T, fixed feed salinity

% Trapezoidal integration: this is int_{Ta}^{Tfeed} Cp dT
delta_h = trapz(T_range, Cp_vals);   % [J/kg] -- specific enthalpy change

Q_feed = mdot_feed * delta_h * (P.t_operating * 3600);   % [J], SIGNED (negative if Tfeed < Ta)
end

function MB = compute_mass_balance(results)
% Final-time water and salt mass-balance residuals, per plate and
% cascade-wide: mdot_in - mdot_out - evaporation - storage should be
% ~0 for water, and mdot_in - mdot_out should be ~0 for salt.
P  = results.P;
iE = numel(results.t);
Np = P.Ns;   % NOTE: local named "Np" for brevity below; it spans all Ns stages, incl. the floor

Tw_f    = results.Tw(:,:,iE);
delta_f = results.delta(:,:,iE);
C_f     = results.C(:,:,iE);
u_f     = results.u(:,:,iE);

rho_f = water_density(Tw_f, C_f);
mdot_w_field = rho_f .* u_f .* delta_f .* P.W;
mdot_s_field = C_f   .* u_f .* delta_f .* P.W;


mdot_w_in_bc = zeros(Np,1);
mdot_s_in_bc = zeros(Np,1);
for k = 1:Np
    if k == 1
        % Same self-consistent inlet the RHS used, at the SOLVED cascade
        % inlet temperature (the hot side).
        Tw_in_k = P.Tfeed_cascade_in;
        [M_in_k, Ms_in_k] = feed_inlet_state(P, Tw_in_k);
    else
       
        M_in_k  = rho_f(k-1,end) * delta_f(k-1,end);
        Ms_in_k = C_f(k-1,end)   * delta_f(k-1,end);
        Tw_in_k = Tw_f(k-1,end);
    end
    [delta_in_k, C_in_k, rho_in_k] = invert_holdup(M_in_k, Ms_in_k, Tw_in_k, P);
    mu_in_k = water_viscosity(Tw_in_k, C_in_k);
    % Mirrors the inlet BC imposed by rhs(): Nusselt closure at the feed
    % for stage 1, upstream SOLVED velocity for every downstream stage.
    % Re-deriving the Nusselt value in this diagnostic would report an
    % inlet flux the solver never used and would make mdot_w_in(k) differ
    % from mdot_w_out(k-1) by the closure error -- a cascade-total
    % residual invisible in the per-stage rows. With the inlet taken from
    % the solved upstream outlet,
    %       mdot_w_in(k) == mdot_w_out(k-1)  to machine precision,
    % so the cascade total tests the PDE solution itself rather than a
    % reporting inconsistency.
    if k == 1
        u_in_k = rho_in_k * P.g * sin(P.theta) * delta_in_k^2 / (3*mu_in_k);
    else
        u_in_k = max(u_f(k-1,end), 0);
    end
    mdot_w_in_bc(k) = M_in_k  * u_in_k * P.W;
    mdot_s_in_bc(k) = Ms_in_k * u_in_k * P.W;
end

MB.time         = results.t(iE);
% LABELLING. mdot_w_in_bc = M*u*W with M = rho*delta is the TOTAL brine
% holdup flux -- water AND salt -- so it is reported as brine, not water.
% Differencing it as though it were water still closes, because salt is
% conserved to ~1e-6 and the salt terms cancel between inlet and outlet,
% but the quantity balanced is total brine. PART B carries the water-only
% balance (S1w); the two blocks apply the same check to two different
% quantities.
MB.mdot_br_in   = mdot_w_in_bc;                     % [kg/s] TOTAL brine
MB.mdot_br_out  = mdot_w_field(:,end);              % [kg/s] TOTAL brine
MB.mdot_s_in    = mdot_s_in_bc;
MB.mdot_s_out   = mdot_s_field(:,end);
% Genuine water flows: total less the conserved salt.
MB.mdot_w_in    = mdot_w_in_bc      - mdot_s_in_bc;
MB.mdot_w_out   = mdot_w_field(:,end) - mdot_s_field(:,end);
MB.mevap_plate  = results.mevap_plate(iE,:).';

nT = numel(results.t);
if nT >= 3
    i0 = iE - 1; i1 = iE;
    dt = results.t(i1) - results.t(i0);
else
    i0 = 1; i1 = iE;
    dt = max(results.t(i1) - results.t(i0), eps);
end

rho_i0 = water_density(results.Tw(:,:,i0), results.C(:,:,i0));
rho_i1 = water_density(results.Tw(:,:,i1), results.C(:,:,i1));

mass_i0 = sum(rho_i0 .* results.delta(:,:,i0), 2) .* P.W .* P.dx_stage;
mass_i1 = sum(rho_i1 .* results.delta(:,:,i1), 2) .* P.W .* P.dx_stage;

MB.storage_rate_total = (mass_i1 - mass_i0) / max(dt, eps);   % [kg/s] brine holdup
% storage_rate_salt is computed below; the water-only storage is formed
% there, once both parts exist, so the residual differences like with
% like: water in, water out, evaporation, WATER storage.

% Read the CONSERVED SALT STATE directly. Rebuilding it as
%     Ms = results.C .* results.delta
% would form the salt holdup from the CLAMPED concentration, which caps
% at P.C_saturation and therefore under-reports the salt inventory
% wherever the clamp is active -- an apparent salt-balance residual that
% the governing equations do not contain (d_Ms has no sink term, so salt
% is conserved exactly by construction).
if isfield(results,'Ms')
    Ms_i0 = results.Ms(:,:,i0);
    Ms_i1 = results.Ms(:,:,i1);
else   % results struct without the stored salt state
    Ms_i0 = results.C_raw(:,:,i0) .* results.delta(:,:,i0);
    Ms_i1 = results.C_raw(:,:,i1) .* results.delta(:,:,i1);
end
MB.storage_rate_salt = sum(Ms_i1 - Ms_i0, 2) .* P.W .* P.dx_stage / max(dt, eps);

% WATER storage = brine storage less salt storage. This is the water term
% the water residual requires; the brine and salt parts are kept
% separately above.
MB.storage_rate = MB.storage_rate_total - MB.storage_rate_salt;

MB.res_water = MB.mdot_w_in - MB.mdot_w_out - MB.mevap_plate - MB.storage_rate;

MB.res_salt = MB.mdot_s_in - MB.mdot_s_out - MB.storage_rate_salt;



MB.rel_water = MB.res_water ./ max(MB.mdot_w_in, eps);
MB.rel_salt  = MB.res_salt  ./ max(MB.mdot_s_in, eps);

MB.res_water_no_storage = MB.mdot_w_in - MB.mdot_w_out - MB.mevap_plate;
MB.rel_water_no_storage = MB.res_water_no_storage ./ max(MB.mdot_w_in, eps);

MB.total.feed_water  = MB.mdot_w_in(1);
MB.total.brine_water = MB.mdot_w_out(end);
% TOTAL-brine terminals, carried alongside the water ones. brine_water is
% water only, so any consumer comparing against results.brine_out.mdot --
% which is TOTAL -- must use brine_total, or it will read the salt
% fraction as a mismatch.
MB.total.feed_total  = MB.mdot_br_in(1);
MB.total.brine_total = MB.mdot_br_out(end);
MB.total.evap_water  = sum(MB.mevap_plate);
MB.total.storage     = sum(MB.storage_rate);
MB.total.res_water   = MB.total.feed_water - MB.total.brine_water - MB.total.evap_water - MB.total.storage;
MB.total.rel_water   = MB.total.res_water / max(MB.total.feed_water, eps);



MB.total.feed_salt   = MB.mdot_s_in(1);
MB.total.brine_salt  = MB.mdot_s_out(end);
MB.total.storage_salt = sum(MB.storage_rate_salt);
MB.total.res_salt = MB.total.feed_salt - MB.total.brine_salt - MB.total.storage_salt;
MB.total.rel_salt = MB.total.res_salt / max(MB.total.feed_salt, eps);

end % compute_mass_balance

function plot_performance(results)
% Figure: cumulative evaporation captured by the counter-current air
% stream vs. plate number, and total freshwater production rate vs. time.
P = results.P;

cum_evap = fliplr(cumsum(fliplr(results.mevap_plate(end,:))));

figure

plot(1:P.Ns, cum_evap*1000, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Plate number (k)');
ylabel('Cumulative evaporation in air stream [g/s]');
title('Cumulative Evaporation vs. Plate Number (air flow: N_p \rightarrow 1)');
grid on;


figure
plot(results.t, results.mfw*1000, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Freshwater production [g/s]');
title('Freshwater Production vs. Time'); grid on;


end % plot_performance

function print_headline(results)
% PART A -- performance summary on the whole-plant control volume.
%
% Every printed number carries a unit and, where the value depends on it,
% a basis: [win] window integral or window mean over the clean operating
% window, [end] the final-instant state. Cascade-only figures sit in their
% own block at the foot of this part and use a different control volume.
P  = results.P;
pr = results.product;
a  = results.awg;
D  = results.downstream;
sysE = [];
if isfield(results,'system'), sysE = results.system; end
t_op_s = P.t_operating*3600;
kgday  = @(m) m*t_op_s;

m_feed_water = P.Vfeed*(P.rhow_in - P.TDSfeed)/1000;              % [kg/day]
m_loop_vapour = kgday(P.mdot_da*P.wv_fan_in);                     % [kg/day]

if isfield(pr,'int_valid') && pr.int_valid
    m_prod_d = pr.kg_per_day_int;  prod_tag = '[win]';
else
    m_prod_d = pr.kg_per_day;      prod_tag = '[end x t_op]';
end

m_salt_d     = P.Vfeed*P.TDSfeed/1000;                            % [kg/day]
m_residual_d = (m_salt_d/P.bl.w_salt_target)*(1 - P.bl.w_salt_target);
rec_max      = 100*(m_feed_water - m_residual_d)/max(m_feed_water,eps);
rec_pct      = 100*m_prod_d/max(m_feed_water,eps);

hdr('PART A -- PERFORMANCE SUMMARY', ...
    'Closed air loop, zero make-up. Control volume: whole plant.');

% ---------------------------------------------------------------- headline
sec('HEADLINE');
prow('Freshwater product',        '%.2f', m_prod_d,               'kg/day', prod_tag);
prow('Seawater recovery',         '%.2f', rec_pct,                '%',      sprintf('%s  ceiling %.2f %%', prod_tag, rec_max));
if ~isempty(sysE)
prow('Specific energy, plant',    '%.2f', sysE.SEC_external,         'kWh/m3', '[win]');
prow('Solar share of energy in',  '%.2f', 100*sysE.Q_solar_absorbed/max(sysE.Q_solar_absorbed+sysE.W_external,eps), ...
                                                                  '%',      '[win]');
end
prow('Brine concentration factor','%.2f', results.brine_out.CF,   '-',      '[end]  cascade outlet');
prow('Salt recovered, anhydrous', '%.2f', kgday(pick(results,'wmean','mdot_salt', D.mdot_salt)), ...
                                                                  'kg/day', '[win]');
if rec_pct > rec_max + 1e-6
    note('Recovery exceeds the slurry ceiling. Seawater is the only water');
    note('inlet in a closed loop, so this is an accounting fault, not yield.');
end

% ---------------------------------------------------------------- geometry
sec('CONFIGURATION');
prow('Cascade stages, plates + floor','%d',   P.Ns,               '-',    '');
prow('Plate length x width',      '%.3f',     P.L,                'm',    sprintf('x %.3f m', P.W));
prow('Wetted plate area, total',  '%.2f',     sum(P.Ap_stage),    'm2',   '');
prow('Glazing aperture area',     '%.2f',     P.Ag,               'm2',   '');
prow('Stack height, front',       '%.3f',     P.stack_front,      'm',    sprintf('rear %.3f m', P.stack_back));
prow('Plate inclination',         '%.2f',     P.theta*180/pi,     'deg',  sprintf('glazing %.2f deg', P.beta*180/pi));
prow('Seawater feed',             '%.2f',     P.Vfeed,            'L/day',sprintf('at %.0f kg/m3', P.TDSfeed));
prow('Air volumetric flow',       '%.3f',     P.Qfan,             'm3/s', '');
prow('Operating window',          '%.2f',     P.t_operating,      'h',    '');

% ------------------------------------------------------------ temperatures
sec('TEMPERATURES');
sub('ambient and feed');
prow('Ambient air',               '%.2f',     P.Ta,               'K',    'input');
prow('Ambient relative humidity', '%.1f',     100*P.RH_amb,       '%',    'input, radiative terms only');
prow('Ground',                    '%.2f',     P.Tground,          'K',    'input');
prow('Seawater supplied',         '%.2f',     P.Tfeed_in,         'K',    'input');

sub('air path, closed loop');
prow('AWG coil outlet',           '%.2f',     P.T_coil,           'K',    'input, saturated');
prow('Air heater outlet',         '%.2f',     P.T_air_in,         'K',    'setpoint, still inlet');
prow('Rise, coil to still inlet', '%.2f',     P.T_air_in-P.T_coil,'K',    'recuperator + steam + electric');
prow('Still inlet humidity ratio','%.5f',     P.wv_fan_in,        'kg/kg','= w_sat(T_coil)');
prow('Cascade exhaust',           '%.2f',     a.Tv_exit,          'K',    'solved');
prow('Exhaust dew point',         '%.2f',     a.T_dew_return,     'K',    'solved');
prow('Dew-point margin at coil',  '%.2f',     a.T_dew_return-a.T_coil, 'K', ...
     ternary(a.T_dew_return-a.T_coil > 0, '[end] condensing', '[end] NO RECOVERY'));
if isfield(P,'T_coil_frosted') && P.T_coil_frosted
    note(sprintf('Coil clamped at the %.1f K frost guard; no frost logic is modelled.', P.T_coil_floor));
end

sub('feed path');
if isfield(results,'feedpre') && results.feedpre.dynamic
    fp = results.feedpre;
    prow('Cascade feed inlet',    '%.2f',     fp.Tfeed_out,       'K',    'solved in RHS');
    prow('Rise across the feed HX','%.2f',    fp.dT,              'K',    ...
         sprintf('binding: %s', ternary(fp.capped,'STARTUP CAP',fp.regime)));
    prow('Limit, supply',         '%.2f',     fp.dT_raw,          'K',    '');
    prow('Limit, design target',  '%.2f',     fp.T_target-fp.Tfeed_in, 'K', '');
    prow('Limit, pinch',          '%.2f',     fp.dT_pinch,        'K',    '');
    prow('Feed capacity rate',    '%.3f',     fp.Cden,            'W/K',  '');
    prow('Condensing steam',      '%.2f',     fp.T_cond,          'K',    'BPE stays in the crystalliser');
    if fp.capped
        note('The startup guard is still clamping at t_end. It is a solver aid,');
        note('so this feed temperature is (Tfeed_in + cap), not a solved result.');
        note('Shorten FC.t_preheat_guard or raise FC.dT_preheat_cap.');
    end
else
    prow('Cascade feed inlet',    '%.2f',     P.Tfeed,            'K',    'input');
end

sub('downstream units');
prow('Cascade brine outlet',      '%.2f',     results.brine_out.T,'K',    'solved');
prow('Crystalliser vapour',       '%.2f',     D.T_vap_heater,     'K',    'element side');
prow('Vapour saturation',         '%.2f',     D.T_sat_vap,        'K',    'at the evaporator pressure');
prow('AWG condenser reject',      '%.2f',     a.T_cond,           'K',    sprintf('Ta + %.1f K approach', P.dT_cond_approach));
prow('Heat-pump lift',            '%.2f',     a.lift,             'K',    sprintf('COP %.2f', a.COP));

sub('cascade film, final state');
prow('Film inlet, stage 1',       '%.2f',     results.Tw(1,1,end),'K',    '');
prow('Film outlet, last stage',   '%.2f',     results.Tw(P.Ns,end,end), 'K', '');
prow('Film minimum, all stages',  '%.2f',     min(min(results.Tw(:,:,end))), 'K', ...
     sprintf('maximum %.2f K', max(max(results.Tw(:,:,end)))));
prow('Glazing',                   '%.2f',     results.Tg(end),    'K',    '');

% ------------------------------------------------------------------ water
sec('WATER');
prow('Seawater fed, water only',  '%.2f',     m_feed_water,       'kg/day','[win]');
prow('Freshwater collected',      '%.2f',     m_prod_d,           'kg/day',sprintf('%s  %.2f L/h', prod_tag, pr.mdot_total*3600));
prow('Seawater recovery',         '%.2f',     rec_pct,            '%',     prod_tag);
prow('Recovery ceiling',          '%.2f',     rec_max,            '%',     sprintf('water held in %.0f wt%% slurry', 100*P.bl.w_salt_target));
prow('Moisture circulating',      '%.2f',     m_loop_vapour,      'kg/day','loop inventory, not an inlet');
note('The loop admits no make-up, so seawater is the only route by which');
note('water enters the plant. There is no atmospheric harvest term.');

% ---------------------------------------------------------- brine and salt
sec('BRINE AND SALT');
prow('Brine leaving the cascade', '%.2f',     results.brine_out.C,'kg/m3', '[end]');
prow('Concentration factor',      '%.2f',     results.brine_out.CF,'-',    '[end] crystalliser inlet');
prow('Salt recovered, anhydrous', '%.2f',     kgday(pick(results,'wmean','mdot_salt', D.mdot_salt)), 'kg/day','[win]');
prow('Salt slurry',               '%.2f',     kgday(pick(results,'wmean','mdot_brine_out', D.mdot_brine_out)), ...
                                                                  'kg/day',sprintf('[win] at %.0f wt%%', 100*D.w_salt_out));
note('ZLD closure occurs in the crystalliser. The cascade alone does not');
note('reach salt saturation and must not be described as achieving ZLD.');

% ----------------------------------------------------------------- energy
if ~isempty(sysE)
    sec('ENERGY');
    sub('purchased work');
    prow('Blower',                '%.2f',     sysE.W_blower,         'W',    '[win]');
    prow('Crystalliser heater',   '%.2f',     sysE.W_heater,         'W',    '[win]');
    prow('AWG compressor',        '%.2f',     sysE.W_compressor,     'W',    sprintf('[win] COP %.2f, lift %.1f K', a.COP, a.lift));
    prow('Electric air heater',   '%.2f',     sysE.W_air_heater,     'W',    sprintf('[win] %.1f K trim', D.dT_air_trim));
    prule();
    prow('Total purchased work',  '%.2f',     sysE.W_external,       'W',    '[win]');
    prow('Absorbed solar, free',  '%.2f',     sysE.Q_solar_absorbed, 'W',    '[win]');

    sub('specific energy on the product');
    prow('Blower',                '%.2f',     sysE.SEC_blower,       'kWh/m3','');
    prow('Crystalliser heater',   '%.2f',     sysE.SEC_heater,       'kWh/m3','');
    prow('AWG compressor',        '%.2f',     sysE.SEC_compressor,   'kWh/m3','');
    prow('Electric air heater',   '%.2f',     sysE.SEC_air_heater,   'kWh/m3','');
    prule();
    prow('SEC, PLANT',            '%.2f',     sysE.SEC_external,     'kWh/m3','[win] quote this figure');
    prow('  band, eta_II 0.55',   '%.2f',     min(sysE.SEC_external_band), 'kWh/m3', ...
         sprintf('to %.2f kWh/m3 at eta_II 0.35', max(sysE.SEC_external_band)));
    note('Reference: RO 3-4 kWh/m3, MED 60-70 kWh/m3 thermal-equivalent.');
    note('The case for this plant rests on brine concentration and off-grid');
    note('operation, not on specific energy.');
end

% ------------------------------------------------------------ cascade only
sec('CASCADE ONLY');
note('Control volume: the cascade alone, excluding crystalliser, feed HX');
note('and AWG. These figures are not the plant figures above.');
prow('Cascade evaporation',       '%.2f',     results.mfw_daily,  'kg/day','[win]');
prow('Cascade evaporation',       '%.4f',     sum(results.mevap_plate(end,:))*1000, 'g/s', '[end]');
prow('Evaporation / water fed',   '%.2f',     results.water_recovery_pct, '%', '[win]');
prow('Brine TDS, last stage',     '%.2f',     squeeze(results.C(P.Ns,end,end)), 'kg/m3','[end]');
prow('Cascade residence time',    '%.2f',     results.retention_time_total/60, 'min','[end]');
prow('Thermal efficiency',        '%.2f',     results.eta_overall,'%',    '[win] E_evap / E_in, ambient datum');
prow('GOR, total boundary supply','%.3f',     results.steady_check.GOR_total_input, '-','[end] quote this figure');
prow('SEC, cascade only',         '%.2f',     results.SEC_kWh_per_m3, 'kWh/m3','[win] different control volume');
note('The cascade SEC counts the intake-air preheat as an input. The work');
note('producing that preheat is already in the plant SEC above; the two');
note('must not be added.');
foot();
end % print_headline


function print_stream_table(results)
% PART B -- stream table and the identities that reconcile it.
%
% Every mass figure elsewhere in the report is one of these rows. All rows
% share one basis, stated in the header, so the g/s and kg/day columns are
% the same quantity in two units. The cascade-node identity is the single
% exception and is labelled where it is printed.
P  = results.P;
a  = results.awg;
D  = results.downstream;
pr = results.product;
t_op_s = P.t_operating*3600;
kgday  = @(m) m*t_op_s;

m_feed   = P.Vfeed*P.rhow_in/1000/t_op_s;                          % [kg/s]
m_feed_w = m_feed - P.Vfeed*P.TDSfeed/1000/t_op_s;                 % [kg/s]
m_harv   = P.mdot_da*(P.wv_fan_in - a.wv_coil);                    % [kg/s]
m_loop   = P.mdot_da*P.wv_fan_in;                                  % [kg/s]

if isfield(results,'wmean') && isfield(results.wmean,'valid') && results.wmean.valid
    wmn        = results.wmean;
    m_evap     = results.mfw_daily/t_op_s;
    m_vap      = wmn.mdot_vap;
    m_slurry   = wmn.mdot_brine_out;
    m_salt_str = wmn.mdot_salt;
    m_brine    = wmn.mdot_brine;
    m_brine_w  = wmn.mdot_brine_water;
    m_recup    = pr.mdot_recup;
    m_coilcond = pr.mdot_awg_air;
    m_reject   = pr.mdot_lost_air;
    tab_basis  = 'window mean over the clean window';
else
    m_evap     = a.mdot_evap_still;
    m_vap      = D.mdot_vap;
    m_slurry   = D.mdot_brine_out;
    m_salt_str = D.mdot_salt;
    m_brine    = results.brine_out.mdot;
    m_brine_w  = results.brine_out.mdot_water;
    m_recup    = D.mdot_cond_recup;
    m_coilcond = a.mdot_condensate;
    m_reject   = a.mdot_moisture_rejected;
    tab_basis  = 'final instant, window integration unavailable';
end

hdr('PART B -- STREAM TABLE', sprintf('Basis: %s.', tab_basis));

srow = @(id,name,from,to,m) fprintf('  %-4s %-24s %-13s %-13s %8.4f %9.2f\n', ...
                                    id, name, from, to, m*1000, kgday(m));
% Stream rows carry their units in the column header; every other printed
% quantity carries its unit on the line.

fprintf('\n  WATER\n');
fprintf('  %-4s %-24s %-13s %-13s %8s %9s\n','ID','stream','from','to','g/s','kg/day');
prule();
srow('S1',  'seawater feed, total',    'boundary',   'cascade',      m_feed);
srow('S1w', '  of which water',        'boundary',   'cascade',      m_feed_w);
srow('S2',  'moisture in loop air',    'loop',       'cascade',      m_loop);
prule();
srow('S3',  'cascade evaporation',     'cascade',    'awg',          m_evap);
srow('S4',  'cascade brine, total',    'cascade',    'crystalliser', m_brine);
srow('S4w', '  of which water',        'cascade',    'crystalliser', m_brine_w);
srow('S4s', '  of which salt',         'cascade',    'crystalliser', m_salt_str);
srow('S5',  'crystalliser vapour',     'crystalliser','feed HX',     m_vap);
srow('S6',  sprintf('salt slurry, %.0f wt%%',100*D.w_salt_out), ...
                                       'crystalliser','boundary',    m_slurry);
srow('S7',  'feed-HX condensate',      'feed HX',    'plasma',       pr.mdot_preheater);
srow('S7b', 'air-preheat condensate',  'air preheat','plasma',       pr.mdot_airheat);
srow('S7c', 'recuperator drain',       'recuperator','plasma',       m_recup);
srow('S8',  'vapour to AWG',           'feed HX',    'awg',          pr.mdot_awg_vapour);
srow('S10', 'AWG coil condensate',     'awg',        'plasma',       m_coilcond);
prule();
srow('S12', 'PRODUCT, all condensate', 'plasma',      'boundary',     pr.mdot_total);
prule();
if abs(m_harv) <= eps && abs(m_reject) <= eps
    note('S9 ambient moisture captured and S11 moisture rejected to');
    note('atmosphere are structurally zero in a closed loop and are omitted.');
else
    srow('S9',  'ambient moisture captured','boundary','awg',        m_harv);
    srow('S11', 'moisture to atmosphere', 'awg',      'boundary',    m_reject);
end

% ------------------------------------------------------------- identities
b4 = results.brine_out;
if isfield(results,'mass_balance') && isfield(results.mass_balance,'total')
    stor_c = results.mass_balance.total.storage;                   % [kg/s] water
else
    stor_c = 0;
end
m_evap_snap = a.mdot_evap_still;
res_casc = m_feed_w - m_evap_snap - b4.mdot_water - stor_c;

fprintf('\n  IDENTITIES\n');
prule();
idn('S1w = S3 + water(S4) + storage', 'cascade node, water only');
ideq(m_feed_w, [m_evap_snap, b4.mdot_water, stor_c], res_casc);
note(sprintf('Final-instant basis: the storage term is a rate from the'));
note(sprintf('discretisation and has no window mean. Salt does not evaporate,'));
note(sprintf('so total brine %.4f g/s enters as water %.4f g/s at w_salt %.4f.', ...
     b4.mdot*1000, b4.mdot_water*1000, b4.w_salt));
note('Storage is holdup still filling the films at t_end, not a closure error.');

fprintf('\n');
idn('S3 = S7c + S10', 'closed loop: what evaporates, condenses');
ideq(m_evap, [m_recup, m_coilcond], m_evap - m_recup - m_coilcond);
note('Two units dry the loop air: the recuperator drains its hot side and');
note('the coil takes the rest. S2 circulates and cancels from both sides.');
note('A residual here is loop moisture inventory changing over the window.');

fprintf('\n');
idn('S5 = S7 + S7b + S8', 'crystalliser vapour splits three ways');
ideq(m_vap, [pr.mdot_preheater, pr.mdot_airheat, pr.mdot_awg_vapour], ...
     m_vap - pr.mdot_preheater - pr.mdot_airheat - pr.mdot_awg_vapour);
note(sprintf('f_cond = %.4f at the feed HX. All three streams are product.', ...
     pick(results,'wmean','f_cond', D.f_cond)));

fprintf('\n');
idn('S4 -> crystalliser', 'one stream, two bases');
note(sprintf('S4 is the last-stage outlet and the crystalliser inlet unaltered.'));
note(sprintf('Yield basis %.4f g/s here; rating basis %.4f g/s in PART C.', ...
     m_brine*1000, results.brine_out.mdot*1000));

% -------------------------------------------------------------------- salt
fprintf('\n  SALT\n');
prule();
prow('In with the feed',    '%.4f', P.Vfeed*P.TDSfeed/1000/t_op_s*1000, 'g/s', ...
     sprintf('%.2f kg/day', P.Vfeed*P.TDSfeed/1000));
prow('Out in the slurry',   '%.4f', m_salt_str*1000, 'g/s', ...
     sprintf('%.2f kg/day, anhydrous', kgday(m_salt_str)));
prow('Film storage rate',   '%.3e', results.mass_balance.total.storage_salt*1000, 'g/s', ...
     ternary(results.mass_balance.total.storage_salt > 0, ...
             '[end] films accumulating', '[end] films discharging'));
note('Salt in and salt out differ because the film salt inventory is');
note('changing, not because salt is lost. The storage-corrected closure is');
note('in the mass-balance block and is of order 1e-6 relative.');
foot();
end % print_stream_table


function print_verification(results)
% =====================================================================
%  PART D -- MODEL VERIFICATION
% ---------------------------------------------------------------------
% Solver audit, not a result. Everything here checks that the model did
% what it was asked to do: closure residuals, interface reciprocity,
% flow regimes, per-plate states, stationarity. None of it belongs in
% the paper's results section; it belongs in an appendix or in a reply
% to a reviewer.
% =====================================================================
P = results.P;

hdr('PART D -- MODEL VERIFICATION', ...
    'Solver audit, not results. Closure, reciprocity, regimes, stationarity.');

sec('RESIDENCE TIME AND TERMINAL SALINITY   (final time)');
fprintf('Total cascade retention time (feed -> Plate %d outlet): %.3f s (%.2f min)\n', ...
    P.Ns, results.retention_time_total, results.retention_time_total/60);
for k = 1:P.Ns
    fprintf('   Plate %2d : %7.4f s   (cumulative: %8.4f s)\n', ...
        k, results.retention_time_plate(k), results.retention_time_cumulative(k));
end

TDS_out_per_plate = squeeze(results.C(:,end,end));
fprintf('\nTDS at outlet of each plate (final time):\n');
for k = 1:P.Ns
    fprintf('   Plate %2d : %.3f kg/m^3\n', k, TDS_out_per_plate(k));
end
foot();

print_plate_table(results);
print_gap_diagnostics(results);
print_reynolds_diagnostics(results);
print_mass_balance(results);
print_energy_balance(results);
Print_AWG_Duty(results);
steady_state_power_balance(results);
print_vapor_pressure_diagnostics(results);
print_system_energy(results);
run_diagnostics(results);
end % print_verification


function print_mass_balance(results)
% Console report of compute_mass_balance()'s results: per-plate and
% cascade-total water/salt residuals, plus the storage-rate contribution.
P  = results.P;
MB = results.mass_balance;

hdr('MASS BALANCE   (final time)', ...
    'Convention: in - out - evap - storage ~ 0 for water; in - out ~ 0 for salt.');
fprintf('  m_in is the TRUE inlet boundary flux, Variable_in(k) =\n');
fprintf('  Variable_out(k-1), matching the PDE boundary condition rather\n');
fprintf('  than the node-1 cell average.\n');

fprintf('\n--- Per-plate residuals ---\n');
fprintf('  (m_w_* are WATER flows: total brine less the conserved salt.\n');
fprintf('   Storage is the water part of the holdup rate. Same basis as\n');
fprintf('   the S1w identity in PART B.)\n');
fprintf('  k |  m_w_in     m_w_out    m_evap    res_water    rel_water\n');
fprintf('    |  [g/s]       [g/s]      [g/s]     [g/s]         [-]\n');
for k = 1:P.Ns
    fprintf('  %2d | %9.4f  %9.4f  %8.5f  %+10.4e  %+9.2e\n', ...
        k, MB.mdot_w_in(k)*1000, MB.mdot_w_out(k)*1000, MB.mevap_plate(k)*1000, ...
           MB.res_water(k)*1000, MB.rel_water(k));
end

fprintf('\n--- Storage-corrected per-plate water residual ---\n');
fprintf('  k |  storage_rate   res_water(no storage)   res_water(WITH storage)\n');
fprintf('    |    [g/s]              [g/s]                     [g/s]\n');
for k = 1:P.Ns
    fprintf('  %2d | %+11.4e     %+14.4e         %+14.4e\n', ...
        k, MB.storage_rate(k)*1000, MB.res_water_no_storage(k)*1000, MB.res_water(k)*1000);
end

fprintf('\n  k |  m_s_in     m_s_out               res_salt    rel_salt\n');
fprintf('    |  [g/s]       [g/s]                 [g/s]         [-]\n');
for k = 1:P.Ns
    fprintf('  %2d | %9.6f  %9.6f              %+10.4e  %+9.2e\n', ...
        k, MB.mdot_s_in(k)*1000, MB.mdot_s_out(k)*1000, ...
           MB.res_salt(k)*1000,  MB.rel_salt(k));
end

fprintf('\n--- Cascade totals ---\n');
fprintf('Water in (feed, plate 1 inlet)          : %.5f g/s\n', MB.total.feed_water*1000);
fprintf('Water out (brine, plate %d outlet)       : %.5f g/s\n', P.Ns, MB.total.brine_water*1000);
fprintf('Water out (distillate, evaporation sum) : %.5f g/s\n', MB.total.evap_water*1000);
fprintf('  Storage rate (water, all stages)        : %+.4e g/s\n', MB.total.storage*1000);
fprintf('  Water residual (in - out - evap - storage): %+.4e g/s   (relative: %+.2e)\n', ...
    MB.total.res_water*1000, MB.total.rel_water);

fprintf('Salt in  (feed)                         : %.6f g/s\n', MB.total.feed_salt*1000);
fprintf('Salt out (brine, plate %d outlet)        : %.6f g/s\n', P.Ns, MB.total.brine_salt*1000);
% The cascade-level salt residual includes the rate of change of salt
% held within the films, which is non-zero while the concentration
% profile is still developing. The residual is therefore reported both
% with and without the storage term; only the storage-corrected value
% is expected to approach machine precision.
fprintf('  Storage rate (salt, all stages)         : %+.4e g/s\n', MB.total.storage_salt*1000);
fprintf('  Salt residual (in - out, NO storage)    : %+.4e g/s\n', ...
    (MB.total.feed_salt - MB.total.brine_salt)*1000);
fprintf('  Salt residual (in - out - storage)      : %+.4e g/s   (relative: %+.2e)\n', ...
    MB.total.res_salt*1000, MB.total.rel_salt);

foot();

end % print_mass_balance

function print_gap_diagnostics(results)
% Console report of final-time vapor-gap state and duct heat/mass
% transfer coefficients (hc, hm) per gap, alongside the ambient-property
% reference values for comparison.
P  = results.P;
iE = numel(results.t);

S.Tg    = results.Tg(iE);
S.Tv    = results.Tv(iE,:).';
S.wv    = results.wv(iE,:).';
S.Tw    = results.Tw(:,:,iE);
S.delta = results.delta(:,:,iE);
S.C     = results.C(:,:,iE);
S.Tp    = results.Tp(:,:,iE);

V_gap = P.V_gap_vec;   % Np x 1, per-compartment (gap 1 differs -- see build_parameters)

% ---- Exact local mass-flow vector at final time ----
mdot_local_vec_iE = P.mdot_da + P.mdot_vapor_in + ...
    flipud(cumsum(flipud(results.mevap_plate(iE,:).')));

hdr('VAPOUR-GAP PROPERTIES   (final time)', ...
    'Per-gap state and duct transfer coefficients against ambient reference.');
fprintf('Reference (ambient / constant-property model):\n');
fprintf('  rho_air(Ta) = %.4f kg/m^3   Cp_air(Ta) = %.2f J/kg-K\n', P.rho_air, P.Cp_air);
fprintf('  mu_air(Ta)  = %.4e Pa*s     k_air(Ta)  = %.4f W/m-K\n\n', P.mu_air, P.k_air);

fprintf('  k |   Tv     wv     |   hc        hm       |  Cp_v     rho_v    Mv\n');
fprintf('    |  [K]    [-]     | [W/m2K]   [m/s]     | [J/kgK] [kg/m3]  [kg]\n');

for k = 1:P.Ns
    [hc_wv, hm_wv] = film_gap_coeffs(k, S, P, mdot_local_vec_iE(k));
    hc_k = mean(hc_wv);
    hm_k = mean(hm_wv);

    Tv_k  = S.Tv(k);
    wv_k  = S.wv(k);
    Cp_k  = moist_air_cp(Tv_k, wv_k);
    rho_k = moist_air_density(Tv_k, wv_k);
    Mv_k  = rho_k * V_gap(k);

    fprintf('  %2d | %6.2f  %.5f | %7.3f  %.4e | %7.2f  %6.4f  %.4e\n', ...
        k, Tv_k, wv_k, hc_k, hm_k, Cp_k, rho_k, Mv_k);
end

foot();

end % print_gap_diagnostics

function print_plate_table(results)
% Console report of final-time per-plate inlet/outlet conditions (film
% temperature, thickness, velocity, TDS, density; plate & vapor gap
% temperatures; solar input and evaporation per plate).
P  = results.P;
iE = numel(results.t);

hdr('PER-STAGE INLET / OUTLET   (final time)', ...
    'Film inlet at x = 0 (top of plate); outlet at x = L.');

% CLOSED LOOP: the bottom gap is fed by the air heater at T_air_in, not by
% ambient, so the table below reports Tv_in(Ns) = T_air_in.
fprintf('Air flows COUNTER-CURRENT to the film: LOOP air enters the bottom gap (k = %d)\n', P.Ns);
fprintf('at T_air_in = %.2f K and is drawn UPWARD by the suction fan, exiting at gap 1.\n', P.T_air_in);
fprintf('Tv_in(k) is therefore Tv_out(k+1) for k < Np, and T_air_in for k = Np.\n\n');

fprintf('--- Water film ---\n');
fprintf('  k  |  Tw_in  Tw_out |  delta_in  delta_out |   u_in     u_out   |  TDS_in  TDS_out |  rho_in  rho_out\n');
fprintf('     |   [K]    [K]   |   [mm]      [mm]     | [mm/s]   [mm/s]    | [kg/m3] [kg/m3]  | [kg/m3] [kg/m3]\n');
for k = 1:P.Ns
    fprintf('  %2d | %6.2f  %6.2f | %8.4f  %8.4f   | %7.4f  %7.4f   | %6.2f  %6.2f   | %7.2f  %7.2f\n', ...
        k, results.Tw_in_true(iE,k), results.Tw_out(iE,k), ...
        results.delta_in_true(iE,k)*1000, results.delta_out(iE,k)*1000, ...
        results.u_in_true(iE,k)*1000,     results.u_out(iE,k)*1000, ...
        results.C_in_true(iE,k),          results.C_out(iE,k), ...
        results.rho_in_true(iE,k),        results.rho_out(iE,k));
end

fprintf('\n--- Absorber plate & vapor gap ---\n');
fprintf('   k |  Tp_in   Tp_out |  Tv_in   Tv_out\n');
fprintf('     |   [K]     [K]   |   [K]     [K]\n');
for k = 1:P.Ns
    fprintf('  %2d | %6.2f  %6.2f  | %6.2f  %6.2f\n', ...
        k, results.Tp_in(iE,k), results.Tp_out(iE,k), ...
        results.Tv_in(iE,k),   results.Tv_out(iE,k));
end

fprintf('\n--- Solar irradiance reaching each water layer & evaporation per plate ---\n');
fprintf('   k |  I_k        Q_solar_w    Q_solar_p    m_evap\n');
fprintf('     | [W/m^2]       [W]           [W]        [g/s]\n');
for k = 1:P.Ns
    fprintf('  %2d | %7.2f    %9.4f    %9.4f    %8.5f\n', ...
        k, results.I_layer(iE,k), ...
           results.q_solar_w_plate(iE,k), ...
           results.q_solar_p_plate(iE,k), ...
           results.mevap_plate(iE,k)*1000);
end



foot();

end % print_plate_table

function print_reynolds_diagnostics(results)
% Console report of final-time film and duct Reynolds numbers per plate,
% with a laminar/wavy-laminar/turbulent flow-regime label for each.
P  = results.P;
iE = numel(results.t);

hdr('FLOW REGIMES   (final time)', ...
    'Film and duct Reynolds numbers with the correlation branch taken.');
fprintf('  k |  Re_film    regime_film   |  Re_gap     regime_gap\n');
fprintf('    |   [-]                     |   [-]\n');
for k = 1:P.Ns
    Ref = results.Re_film(iE,k);
    Reg = results.Re_gap(iE,k);

    if Ref < 30
        regF = 'laminar (Nusselt)';
    elseif Ref < 1600
        regF = 'wavy-laminar';
    else
        regF = 'turbulent';
    end

    if Reg < 2300
        regG = 'laminar';
    else
        regG = 'turbulent';
    end

    fprintf('%s\n', deblank(sprintf('  %2d | %9.4f  %-18s | %9.2f  %-10s', ...
            k, Ref, regF, Reg, regG)));
end
foot();

end % print_reynolds_diagnostics

function print_energy_balance(results)
P    = results.P;
iE   = numel(results.t);
Np   = P.Ns;   % NOTE: local named "Np" for brevity below; it spans all Ns stages, incl. the floor (stage Np==Ns)
t=results.t;
hdr('ENERGY BALANCE BY CONTROL VOLUME   (final time)', ...
    'Each block sums to zero to within the solver tolerance.');

% ---- Exact instantaneous state derivative at the final time, used below
% for every storage term (glass, wall nodes, floor) -- computed ONCE here
% (single point of truth) rather than re-calling rhs() separately in each
% section.
dY_end = rhs(t(iE), results.Y_final, P);

Tv1   = results.Tv(iE,1);   Tg = results.Tg(iE);

% ---- hc_vg recomputed here EXACTLY as glass_rhs() computes it internally
% (flat P.Acomp/P.Dh reference), NOT taken from results.hc_vg -- that
% diagnostic field uses a different, more detailed node-resolved taper
% formula (see extract_results), so it does not exactly equal the term
% the ODE solver actually used for d_Tg. Recomputing here, matching
% glass_rhs precisely, is what lets the storage-corrected SUM below close
% to ~0, following the same convention as the floor and wall sections:
% each diagnostic reproduces the coefficient its own ODE used.
% The glazing receives its convective exchange with gap 1 as an absolute
% power evaluated by vapor_gap_rhs from film_gap_coeffs(1,...). This
% diagnostic reproduces that evaluation so the balance below closes.
mdot_local_gap1_iE = P.mdot_da + P.mdot_vapor_in + sum(results.mevap_plate(iE,:));
S_g.Tv = results.Tv(iE,:).';
S_g.Tw = results.Tw(:,:,iE);
hc_vg  = mean(film_gap_coeffs(1, S_g, P, mdot_local_gap1_iE));

q_solar_g   =  P.alpha_g * results.I(iE) * P.Ag;
q_rad_wg    =  results.q_glass_rad_wg(iE);
q_conv_vg   =  hc_vg * (Tv1 - Tg) * P.Ag;
q_conv_amb  = -results.q_glass_conv(iE);
q_rad_sky   = -results.q_glass_rad_sky(iE);

% ---- Storage term: the glass has real thermal capacitance (P.mCp_g) and
% Tg is a genuine dynamic ODE state, NOT assumed quasi-steady -- exactly
% like the floor below. Early in a transient run the glass is still
% actively heating/cooling, so the raw flux imbalance is NOT expected to
% vanish on its own; this term closes that gap explicitly, the same way
% the floor's storage-corrected closure does further down.
d_Tg_true   = dY_end(P.idx.Tg);         % [K/s], areal formula (see glass_rhs)
q_storage_g = P.mCp_g * P.Ag * d_Tg_true;   % [W], + = glass heating up

sum_g_quasi_steady = q_solar_g + q_rad_wg + q_conv_vg + q_conv_amb + q_rad_sky;
sum_g              = sum_g_quasi_steady - q_storage_g;

fprintf('\n--- Glass cover ---\n');
fprintf('  q_solar_g   (absorbed solar)   : %+10.3f W\n', q_solar_g);
fprintf('  q_rad_wg    (from film 1)      : %+10.3f W\n', q_rad_wg);
fprintf('  q_conv_vg   (from gap 1)       : %+10.3f W\n', q_conv_vg);
fprintf('  q_conv_gamb (to ambient)       : %+10.3f W\n', q_conv_amb);
fprintf('  q_rad_gsky  (to sky)           : %+10.3f W\n', q_rad_sky);
fprintf('  q_storage   (glass heating up) : %+10.3f W\n', q_storage_g);
fprintf('  SUM (should be ~0)             : %+10.3f W\n', sum_g);


fprintf('\n--- Water films (per stage; stage %d is the floor''s own film) ---\n', P.Ns);
fprintf('  k |    q_sol     q_p2w     q_wv       q_evap      q_rad      adv_net       SUM\n');
fprintf('    |     [W]       [W]      [W]         [W]         [W]         [W]         [W]\n');


for k = 1:Np
    qs = +results.q_solar_w_plate(iE,k);
    qp = +results.q_pw_plate(iE,k);
    qc = -results.q_conv_wv_plate(iE,k);
    qe = -results.q_evap_plate(iE,k);
    qr = -results.q_rad_wp_plate(iE,k);

    % ---- Node-resolved advective term (replaces 2-point lumped adv) ----
    Tw_k    = results.Tw(k,:,iE);
    delta_k = max(results.delta(k,:,iE), 1e-8);
    u_k     = results.u(k,:,iE);
    rho_w_k = results.rho_w(k,:,iE);
    Cp_w_k  = results.Cp_w(k,:,iE);

    if k == 1
        % HOT SIDE. The film enters stage 1 at the preheat HX outlet, not
        % at the supply temperature. Using P.Tfeed here understates the
        % advective inflow by mdot*Cp*(Tfeed_cascade_in - Tfeed) and
        % leaves the stage-1 row failing to close by that amount while
        % every other stage closes. The RHS always uses the solved
        % value; only this diagnostic can disagree.
        Tw_in_k = P.Tfeed_cascade_in;
    else
        Tw_in_k = results.Tw_out(iE,k-1);
    end

    dTw_dx = upwind_deriv_1field(Tw_k, P.dx_stage(k), Tw_in_k);   % same stencil as rhs()

    % local areal advective rate [W/m^2], sign matches "-u*dTw_dx" in rhs()
    adv_areal = -rho_w_k .* Cp_w_k .* delta_k .* u_k .* dTw_dx;

    adv = sum(adv_areal) * P.dx_stage(k) * P.W;                    % [W]

    row = qs + qp + qc + qe + qr + adv;
    fprintf(' %3d | %+7.3f   %+7.3f   %+7.3f   %+9.3f   %+8.3f   %+8.3f   %+9.4f\n', ...
        k, qs, qp, qc, qe, qr, adv, row);
end
fprintf('\n--- Absorber plates (real plates only, k=1..%d) ---\n', P.Np);
fprintf('  k |  q_sol_p   q_p2w    q_conv_vp   q_rad_gain    SUM\n');
fprintf('    |    [W]      [W]        [W]         [W]        [W]\n');
for k = 1:P.Np
    qs = +results.q_solar_p_plate(iE,k);
    qp = -results.q_pw_plate(iE,k);
    qv = +results.q_conv_vp_plate(iE,k);
    qg = +results.q_rad_gain_plate(iE,k);
    row = qs + qp + qv + qg;
    fprintf(' %3d | %+8.3f %+8.3f %+10.3f %+10.3f %+9.4f\n', ...
        k, qs, qp, qv, qg, row);
end
fprintf('(Stage %d, the floor, is NOT shown here -- its own conductive back-\n', P.Ns);
fprintf(' side losses (ground + wall-strip) are a different physical\n');
fprintf(' mechanism than this table''s "q_rad_gain" column measures for a\n');
fprintf(' real plate, so it gets its own breakdown further below instead\n');
fprintf(' of a misleading row here.)\n');

fprintf('\n--- Vapor gaps (per gap) ---\n');
fprintf('  k |  H_adv_in    H_adv_out    q_conv_film   q_conv_up    q_latent     q_wall      SUM\n');
fprintf('    |    [W]         [W]           [W]           [W]         [W]         [W]        [W]\n');

evap_per_plate = results.mevap_plate(iE,:).';

% Inlet state of the fan stream. The stored humidity ratio is used
% directly rather than reconstructed from P.phi_fan_in, because the
% relative humidity at the still inlet is itself a DERIVED quantity
% (P.phi_fan_in = w_sat(T_coil)/w_sat(Tfan_in)) and the mapping between
% w and phi is nonlinear:
%
%     w = 0.622*phi*Psat / (Patm - phi*Psat)   is not   phi * w_sat
%
% so the round trip w -> phi -> w does not return the original value.
% The discrepancy is only a few per cent in w, but it enters the
% terminal gap's inlet enthalpy Min*Cp_in*(Tin - Tref), which is of
% order 24 kW, and therefore appears as a residual of order 10 W in
% that gap's balance alone -- the terminal gap being the only one whose
% inlet is the fan stream rather than the gap below it.
w_fan_in  = P.wv_fan_in;
Cp_fan_in = moist_air_cp(P.Tfan_in, w_fan_in);

Awall_vec = wall_area_per_gap(P);

for j = 1:Np

    Tv_j = results.Tv(iE,j);
    wv_j = results.wv(iE,j);
    Cp_j = moist_air_cp(Tv_j, wv_j);

    if j < Np
        Tin = results.Tv(iE,j+1);
        wIn = results.wv(iE,j+1);
        Cp_in = moist_air_cp(Tin, wIn);
        Min  = P.mdot_da + P.mdot_vapor_in + sum(evap_per_plate(j+1:end));
    else
       Tin = P.Tfan_in; Cp_in = Cp_fan_in; Min = P.mdot_da + P.mdot_vapor_in;   % inlet mass flow includes the vapour carried by the inlet air
    end
    Mout = P.mdot_da + P.mdot_vapor_in + sum(evap_per_plate(j:end));

    % Enthalpy flux terms are referenced to P.Tref, exactly as in the
    % solved ODE (vapor_gap_rhs, d_Tv), rather than to absolute T.
    % ASSUMPTION: since Min != Mout (Mout picks up evap_per_plate(j) more
    % than Min), referencing to absolute T would introduce a spurious
    % residual of order Tref*Cp*(Min-Mout) -- this diagnostic must mirror
    % the same convention as the solved equation to avoid that.
    H_in  = +Min  * Cp_in * (Tin - P.Tref);
    H_out = -Mout * Cp_j  * (Tv_j - P.Tref);
    % The gap receives the same nodal integrals used by the film and
    % plate balances. The lumped alternative
    % mean(hc)*(mean(Tw)-Tv)*W*L_stage differs from these by the
    % covariance mean(hc'*Tw'), and for the terminal gap also by the
    % interfacial area, since the partner plate spans L rather than
    % floor_L.
    if P.opt.interface_consistent
        q_from_film = +results.q_conv_wv_plate(iE,j);      % film j -> gap j
        if j == 1
            % gap 1's top is the glass: same expression glass_rhs receives
            q_to_above = -hc_vg * (Tv_j - results.Tg(iE)) * P.Ag;
        else
            % gap j's top is plate (j-1): the value plate (j-1) gained
            q_to_above = -results.q_conv_vp_plate(iE,j-1);
        end
    else
        q_from_film = +results.hc_wv(iE,j) * ...
                       (mean(results.Tw(j,:,iE)) - Tv_j) * P.W * P.L_stage(j);
        if j == 1
            Tabove = results.Tg(iE);
        else
            Tabove = mean(results.Tp(j-1,:,iE));
        end
        q_to_above = -results.hc_wv(iE,j) * (Tv_j - Tabove) * P.W * P.L_stage(j);
    end
    Cp_vapor_j = 1860 + 0.12*(mean(results.Tw(j,:,iE)) - 273.15); 
    q_lat = + evap_per_plate(j)*Cp_vapor_j*(mean(results.Tw(j,:,iE))-P.Tref) ;

    % Reuses the already-stored forced-convection coefficient
    % (results.hc_wv, computed by film_gap_coeffs inside extract_results)
    % rather than calling film_gap_coeffs again, since no S struct is in
    % scope here and hc_wv(iE,j) is that same coefficient.
    Twall_j = results.Twall(iE,j);
    q_wall_j = -results.hc_wv(iE,j) * (Tv_j - Twall_j) * Awall_vec(j);   % vapor -> wall, negative = loss

    row = H_in + H_out + q_from_film + q_to_above + q_lat + q_wall_j;

    fprintf(' %3d | %+10.3f %+11.3f %+13.3f %+12.3f %+10.3f %+10.3f %+9.4f\n', ...
        j, H_in, H_out, q_from_film, q_to_above, q_lat, q_wall_j, row);
end

fprintf('\n--- Wall / insulation nodes (per gap) ---\n');
fprintf('(mCp_wall_areal*Awall,j*dTwall/dt), since the wall node has real capacitance.\n');
fprintf('The q_rad column is the wall''s net radiative loss to the\n');
fprintf('film+top-surface enclosure (radiosity_all_gaps), in addition to\n');
fprintf('its convective and conductive terms.\n\n');
fprintf('  k |    Tv       Twall       Tamb   |   q_in(v->wall)   q_out(wall->amb)   q_rad(wall->encl)   q_cond_net   q_storage      SUM\n');
fprintf('    |   [K]        [K]        [K]    |        [W]              [W]                [W]             [W]          [W]         [W]\n');

Awall_vec = wall_area_per_gap(P);
S_iE.Tv = results.Tv(iE,:).';
S_iE.Tw = results.Tw(:,:,iE);

% Axial conduction along the continuous wall sheet, mirroring wall_rhs()
% exactly, so the "SUM" column below reflects the full energy balance
% (q_in, q_out, q_rad, q_cond, and storage together) rather than just
% q_in - q_out.
A_cond = 2*(P.chamber_L + P.W) * P.t_wall;
L_cond = (P.hgap_narrow + P.hgap_wide) / 2;
R_link = (L_cond(1:end-1)/2 + L_cond(2:end)/2) ./ (P.k_wall * A_cond);   % [K/W], per link

Twall_iE  = results.Twall(iE,:).';

% ---- Floor's own wall-strip conduction (a GAIN unique to Twall(Ns), the
% wall wrapping the floor's edge) -- same formula as rhs()'s k==Ns branch
% and wall_rhs()'s Q_floor_to_wallstrip argument. ----
Tbase_iE = results.Tbase(iE);
Q_floor_to_wallstrip_iE = (Tbase_iE - Twall_iE(P.Ns)) * P.UA_floor_wall;

% Link-based, matching wall_rhs exactly: each link carries one heat flow
% seen with opposite sign by the two nodes it joins, so the axial terms
% sum to zero across the network.
q_cond = zeros(Np,1);   % net conduction INTO each wall node, [W]
for j = 1:(Np-1)
    q_link = (Twall_iE(j) - Twall_iE(j+1)) / R_link(j);
    q_cond(j)   = q_cond(j)   - q_link;
    q_cond(j+1) = q_cond(j+1) + q_link;
end
q_cond(P.Ns) = q_cond(P.Ns) + Q_floor_to_wallstrip_iE;

% Exact instantaneous storage rate for each wall node, obtained by
% differentiating the actual RHS (same technique used for the floor
% below), rather than assuming steady state (SUM = 0).
d_Twall_true = dY_end(P.idx.Twall);

% ---- Exact local mass-flow vector at final time ----
mdot_local_vec_iE = P.mdot_da + P.mdot_vapor_in + ...
    flipud(cumsum(flipud(results.mevap_plate(iE,:).')));

for j = 1:Np
    Tv_j    = results.Tv(iE,j);
    Twall_j = Twall_iE(j);

    % Same forced-convection coefficient as film_gap_coeffs, consistent
    % with vapor_gap_rhs and wall_rhs (P.h_in is not used here).
    hc_wall_j = mean(film_gap_coeffs(j, S_iE, P, mdot_local_vec_iE(j)));
    q_in_j      = hc_wall_j * (Tv_j - Twall_j) * Awall_vec(j);          % vapor -> wall
    q_out_j     = (Twall_j - P.Ta) / P.R_wall_out * Awall_vec(j);        % wall -> ambient
    q_rad_j     = results.q_rad_wall_plate(iE,j);                        % wall -> film+top enclosure (net loss)
    q_storage_j = P.mCp_wall_areal * Awall_vec(j) * d_Twall_true(j);     % + = wall node heating up

    row = q_in_j - q_out_j - q_rad_j + q_cond(j) - q_storage_j;

    fprintf(' %3d | %7.2f   %8.2f   %7.2f  | %+14.4f    %+14.4f      %+13.4f    %+10.4f   %+10.4f   %+10.4f\n', ...
        j, Tv_j, Twall_j, P.Ta, q_in_j, q_out_j, q_rad_j, q_cond(j), q_storage_j, row);
end



% ---- Floor: full breakdown, storage-corrected ----
% Recomputes the SAME q_solar_p / q_p2w already shown in the generic
% per-stage table above (so the numbers are traceable), but reports the
% back-side loss as two SEPARATE conduction terms (ground, wall-strip)
% rather than as one lumped "q_rad_gain" entry:
% there is no radiative term here at all; the water
% film sitting on the floor has no bearing on this path, since it's the
% floor's SOLID back side (not the wetted top) that loses heat this way.
q_solar_p_floor = results.q_solar_p_plate(iE, P.Ns);
q_p2w_floor      = -results.q_pw_plate(iE, P.Ns);   % sign convention matches the table above (loss FROM the floor solid)

Q_floor_to_ground_iE    = (Tbase_iE - P.Tground) / P.R_floor_ground * P.Ap_floor;
Q_floor_to_wallstrip_iE = (Tbase_iE - Twall_iE(P.Ns)) * P.UA_floor_wall;

d_Tp_true = reshape(dY_end(P.idx.Tp), P.Ns, P.Nx);   % dY_end computed once at the top of this function
d_Tbase_mean_true = mean(d_Tp_true(P.Ns,:));
q_storage_floor = P.rho_base*P.Cp_base*P.t_base*P.Ap_floor * d_Tbase_mean_true;   % + = floor heating up

sum_floor = q_solar_p_floor + q_p2w_floor - Q_floor_to_ground_iE - Q_floor_to_wallstrip_iE - q_storage_floor;

fprintf('\n--- Bottom base (floor) energy balance, stage %d ---\n', P.Ns);
fprintf('  q_solar_p    (absorbed, after water film)  : %+10.3f W\n', q_solar_p_floor);
fprintf('  q_p2w        (into the water film above)   : %+10.3f W\n', q_p2w_floor);
fprintf('  q_to_ground  (conduction, floor->ground)    : %+10.3f W\n', -Q_floor_to_ground_iE);
fprintf('  q_to_wallstrip (conduction, floor->wall)    : %+10.3f W\n', -Q_floor_to_wallstrip_iE);
fprintf('  q_storage    (floor heating up)             : %+10.3f W\n', q_storage_floor);
fprintf('  SUM (should be ~0)                          : %+10.3f W\n', sum_floor);




% ===================================================================
% INTERFACE RECIPROCITY CHECK
% Each convective pair must represent a single power seen identically by
% both adjoining control volumes. Under the nodal formulation
% (P.opt.interface_consistent = true) these residuals vanish identically;
% under the lumped alternative they quantify the resulting imbalance.
% ===================================================================
fprintf('\n--- Interface reciprocity check (film/plate <-> vapor gap) ---\n');
fprintf('  Each pair must be ONE number seen by both control volumes.\n');
fprintf('  pair                     |   film/plate side |     gap side |   residual\n');
max_resid = 0;
for j = 1:Np
    a = results.q_conv_wv_plate(iE,j);
    if P.opt.interface_consistent
        b = a;
    else
        b = results.hc_wv(iE,j) * (mean(results.Tw(j,:,iE)) - results.Tv(iE,j)) * P.W * P.L_stage(j);
    end
    r = a - b;  max_resid = max(max_resid, abs(r));
    fprintf('  film %2d -> gap %2d        | %+17.6f | %+12.6f | %+10.3e\n', j, j, a, b, r);
end
for j = 2:Np
    a = results.q_conv_vp_plate(iE,j-1);
    if P.opt.interface_consistent
        b = a;
    else
        b = results.hc_wv(iE,j) * (results.Tv(iE,j) - mean(results.Tp(j-1,:,iE))) * P.W * P.L_stage(j);
    end
    r = a - b;  max_resid = max(max_resid, abs(r));
    fprintf('  gap %2d  -> plate %2d      | %+17.6f | %+12.6f | %+10.3e\n', j, j-1, a, b, r);
end
fprintf('  MAX |residual| = %.4e W', max_resid);
if max_resid < 1e-9
    fprintf('   [PASS]\n');
else
    fprintf('   [FAIL -- interfaces are not conservative]\n');
end

% ===================================================================
% SOLAR CLOSURE
% Accounts for the disposal of the entire incident beam: glazing
% reflection and absorption, transmission through the aperture,
% spillage beyond the plate footprint, absorption within the films and
% solid surfaces, reflection from those surfaces, and any residual
% transmitted past the terminal stage.
% ===================================================================
Ig_iE = results.I(iE);
fprintf('\n--- Solar closure (final time) ---\n');
E_admit   = P.tau_g * Ig_iE * P.Ag;
E_target  = P.tau_g * Ig_iE * P.Ap;
E_spill   = E_admit - E_target;
E_films   = sum(results.q_solar_w_plate(iE,:));
E_plates  = sum(results.q_solar_p_plate(iE,:));
E_refl    = sum(results.q_solar_refl_plate(iE,:));
E_thru    = results.q_solar_terminal(iE);
fprintf('  incident on aperture   Ig*Ag           : %+10.3f W\n', Ig_iE*P.Ag);
fprintf('  reflected by glass     RF_g*Ig*Ag      : %+10.3f W\n', P.RF_g*Ig_iE*P.Ag);
fprintf('  absorbed by glass      alpha_g*Ig*Ag   : %+10.3f W\n', P.alpha_g*Ig_iE*P.Ag);
fprintf('  admitted (transmitted) tau_g*Ig*Ag     : %+10.3f W\n', E_admit);
fprintf('    -> targeted on plate 1 (tau_g*Ig*Ap) : %+10.3f W\n', E_target);
fprintf('    -> APERTURE SPILL (no target stage)  : %+10.3f W   (%.1f%% of admitted)\n', ...
    E_spill, 100*E_spill/max(E_admit,eps));
if P.opt.route_solar_spill
    fprintf('       [routed to wall node 1, alpha = %.2f -]\n', P.opt.alpha_wall_solar);
else
    fprintf('       (accounted only; set P.opt.route_solar_spill to deposit on the wall)\n');
end
fprintf('  of the targeted stream:\n');
fprintf('    absorbed in films                    : %+10.3f W\n', E_films);
fprintf('    absorbed in plates/floor             : %+10.3f W\n', E_plates);
fprintf('    reflected off solid surfaces (lost)  : %+10.3f W\n', E_refl);
fprintf('    transmitted past terminal stage      : %+10.3f W\n', E_thru);
fprintf('    closure residual (must be ~0)        : %+10.3e W\n', ...
    E_target - E_films - E_plates - E_refl - E_thru);

end % print_energy_balance

function check_steady_state(results, k, window_sec, tol)
t   = results.t;
P   = results.P;
if k == P.Ns
    x = results.x_base;   % floor stage is floor_L-long (= L+clearance), not L-long -- see extract_results
else
    x = results.x;
end

win = t >= (t(end) - window_sec);
nW  = nnz(win);
if k == P.Ns
    hdr(sprintf('QUASI-STEADY CHECK, STAGE %d (THE FLOOR)', k), ...
        'State drift against forcing drift over a trailing window.');
else
    hdr(sprintf('QUASI-STEADY CHECK, PLATE %d', k), ...
        'State drift against forcing drift over a trailing window.');
end
fprintf('Window: last %.1f s of run   (%.2f h, %d samples)\n', ...
        window_sec, window_sec/3600, nW);
fprintf('Cascade residence time: %.1f s -- the window MUST exceed this,\n', ...
        results.retention_time_total);
fprintf('or the slowest mode (the salt field) is invisible to this test.\n');
if window_sec < 2*results.retention_time_total
    fprintf(['** WARNING: window is shorter than 2 residence times. This\n' ...
             '   check cannot resolve the salt field and will pass\n' ...
             '   regardless of whether the cascade has converged. **\n']);
end
if nW < 20
    fprintf(['** WARNING: only %d samples in the window. Increase\n' ...
             '   FC.stepSize before trusting the drift statistics. **\n'], nW);
end

% ---- Criterion selection: ABSOLUTE drift vs. FORCING-TRACKING ratio ----
% These are two different tests and only one of them is valid at a time.
%
%   Constant irradiance (autonomous RHS): a true fixed point exists, so
%   the meaningful test is whether each state has stopped moving --
%   drift_rel < tol, an ABSOLUTE criterion.
%
%   Measured irradiance (non-autonomous RHS): no state can ever stop
%   moving, because every one tracks the diurnal forcing. The meaningful
%   test is instead whether the state follows the forcing without lag,
%   i.e. whether its fractional drift is commensurate with the
%   fractional drift of Ig over the same window.
%
% The tracking ratio R = drift_rel/I_drift_rel is therefore applied ONLY
% in the non-autonomous case. Under constant irradiance I_drift_rel = 0
% exactly, so R = drift_rel/eps -> O(1e16) and every field would be
% reported LAGGING however well converged the run is -- the ratio carries
% no information about an autonomous system.
use_absolute = P.opt.constant_irradiance;

tw_I  = t(win); tw_I = tw_I(:);
Iw    = results.I(win); Iw = Iw(:);
pI    = polyfit(tw_I - tw_I(1), Iw, 1);
I_drift_rel = abs(pI(1)) * window_sec / max(abs(mean(Iw)), eps);
R_tol = 3;    % state may drift up to 3x the fractional forcing drift

if use_absolute
    fprintf('Forcing: CONSTANT (Ig = %.1f W/m^2); residual forcing drift %.2e\n', ...
        P.G_const, I_drift_rel);
    fprintf('Criterion: ABSOLUTE -- steady if state drift < %.1e  [-]\n\n', tol);
else
    fprintf('Forcing (Ig) relative drift over window: %.2e  [-]\n', I_drift_rel);
    fprintf('Quasi-steady if state drift < %.0fx forcing drift  [-]\n\n', R_tol);
end

fields = {'Tw','delta','C','mevap'};
units  = {'K','m','kg/m^3','kg/m^2 s'};

is_steady_all = true;
mean_store    = struct();
err_store     = struct();

for f = 1:numel(fields)
    F   = results.(fields{f});
    Fk  = squeeze(F(k,:,win));
    Fm  = mean(Fk, 2);
    Fe  = max(abs(Fk - repmat(Fm, 1, size(Fk,2))), [], 2);
    rel = Fe ./ max(abs(Fm), eps);

    tw  = t(win); tw = tw(:);
    slope = zeros(numel(Fm),1);
    for j = 1:numel(Fm)
        p = polyfit(tw - tw(1), Fk(j,:).', 1);
        slope(j) = p(1);
    end
    drift_rel = abs(slope) * window_sec ./ max(abs(Fm), eps);

    % ---- Apply the criterion selected above ----
    if use_absolute
        % Autonomous system: test convergence to the fixed point directly.
        quasi = max(drift_rel) < tol;
        fprintf('  %-6s : max rel dev = %.2e   max rel drift = %.2e   (tol %.1e)   %s\n', ...
            fields{f}, max(rel), max(drift_rel), tol, ...
            ternary(quasi, 'STEADY', 'NOT CONVERGED'));
    else
        % Non-autonomous system: the tracking ratio
        %     R = (relative drift of the state) / (relative drift of Ig)
        % is O(1) when the response is quasi-steady and grows large when
        % thermal inertia causes the state to lag the forcing.
        R = drift_rel ./ max(I_drift_rel, eps);
        quasi = max(R) < R_tol;
        fprintf('  %-6s : max rel dev = %.2e   max rel drift = %.2e   drift/forcing = %6.2f   %s\n', ...
            fields{f}, max(rel), max(drift_rel), max(R), ...
            ternary(quasi, 'QUASI-STEADY', 'LAGGING'));
    end
    is_steady_all = is_steady_all && quasi;

    mean_store.(fields{f}) = Fm;
    err_store.(fields{f})  = Fe;
end

if use_absolute
    fprintf('\nOverall: %s\n', ternary(is_steady_all, ...
        'converged to steady state (all fields stationary within tolerance)', ...
        'NOT converged -- increase P.t_sim (slowest mode is the terminal-stage salt field)'));
else
    fprintf('\nOverall: %s\n', ternary(is_steady_all, ...
        'response is quasi-steady (state tracks the forcing without lag)', ...
        'one or more variables lag the forcing'));
end
foot();

if k == P.Ns
    figure('Name', sprintf('Quasi steady-state profiles - Stage %d (Floor)', k), 'Color', 'w');
else
    figure('Name', sprintf('Quasi steady-state profiles - Plate %d', k), 'Color', 'w');
end
for f = 1:numel(fields)
    subplot(2, 2, f);
    Fm = mean_store.(fields{f});
    Fe = err_store.(fields{f});
    errorbar(x, Fm, Fe, 'LineWidth', 1.4, 'CapSize', 3);
    xlabel('x [m]');
    ylabel(sprintf('%s [%s]', fields{f}, units{f}));
    title(sprintf('%s: mean \\pm max dev (last %.0f s)', fields{f}, window_sec));
    grid on;
end
try
    if k == P.Ns
        sgtitle(sprintf('Stage %d (Floor) - Quasi steady-state profiles', k));
    else
        sgtitle(sprintf('Plate %d - Quasi steady-state profiles', k));
    end
catch
end

end % check_steady_state

function steady_state_power_balance(results)
steady_check=results.steady_check ;
P = results.P;
hdr('QUASI-STEADY POWER BALANCE', 'Closure verification, not a reported result.');
fprintf('  Q_solar, incident on aperture     = %8.2f W\n', steady_check.Q_solar_aperture);
fprintf('  Q_solar, absorbed by the system   = %8.2f W\n', steady_check.Q_solar_absorbed);
fprintf('  Q_solar, not absorbed             = %8.2f W   (glazing reflection, aperture spill,\n', ...
    steady_check.Q_solar_unused);
fprintf('                                                  internal surface reflection)\n');
fprintf('  Q_fan_thermal, intake air         = %8.2f W\n', steady_check.Q_fan_thermal);
fprintf('  Q_feed_thermal                    = %8.2f W\n', steady_check.Q_feed_thermal);
fprintf('  Q_exhaust, air enthalpy out       = %8.2f W\n', steady_check.Q_vapor_out);
fprintf('  Q_brine, liquid enthalpy out      = %8.2f W\n', steady_check.Q_feed_out);
fprintf('  Q_wall loss                       = %8.2f W\n', steady_check.Q_wall_loss);
fprintf('  Q_glazing reflection              = %8.2f W\n', steady_check.Q_reflect_glass);
fprintf('  Q_glazing convection to ambient   = %8.2f W\n', steady_check.Q_glass_conv_loss);
fprintf('  Q_glazing radiation to sky        = %8.2f W\n', steady_check.Q_glass_rad_loss);
fprintf('  Q_floor conduction to ground      = %8.2f W\n', steady_check.Q_ground_loss);
fprintf('  Q_wall-strip conduction           = %8.2f W\n', steady_check.Q_wallstrip_loss);
fprintf('  ---------------------------------------------------\n');
fprintf('  Boundary streams (positive = leaving, ambient datum):\n');
fprintf('  Q_in,  supplied (solar + air + feed) = %8.2f W\n', steady_check.Q_in_total);
fprintf('  Q_out, net across the boundary       = %8.2f W\n', steady_check.Q_out_total);
if steady_check.Q_gain_total > 0
    fprintf('    of which reversed (entering)      = %8.2f W   (cascade below ambient)\n', ...
        steady_check.Q_gain_total);
end
fprintf('  Q_evaporation (mdot*h_fg)         = %8.2f W   (mdot_evap = %.4f g/s)\n', ...
    steady_check.Q_evap, steady_check.mdot_evap_total*1000);
fprintf('  Q_storage, sensible heat rate     = %8.2f W   (all control volumes)\n', ...
    steady_check.Q_storage_total);
fprintf('  ---------------------------------------------------\n');
fprintf('  Closure: Q_in - Q_out - Q_evap - Q_storage\n');
fprintf('                                    = %8.2f W   (%.2f %% of Q_in)\n', ...
    steady_check.Q_balance_residual, 100*steady_check.Q_balance_rel);
fprintf('  Q_entering, total (supplied + reversed) = %8.2f W\n', steady_check.Q_entering_total);
fprintf('  (This block is a closure verification, not a reported result.\n');
fprintf('   The thermal efficiency quoted in the summary is the\n');
fprintf('   time-integrated E_evap / E_in over the operating period.)\n\n');

% ===================================================================
% AIR-STREAM ENTHALPY CONVENTION CHECK
% ---------------------------------------------------------------------
% The inlet and outlet enthalpies of the air stream must be constructed
% on the same basis. Two conventions are available and they are not
% interchangeable:
%
%   per-component : mdot_da*Cp_da*(T-Ta) + mdot_vap*Cp_vap*(T-Ta)
%   moist-mixture : mdot_moist * Cp_m * (T-Ta),  Cp_m = (Cp_da + w*Cp_v)/(1+w)
%
% They agree only when the moist-mixture form is paired with the MOIST
% mass flow. Pairing Cp_m with the dry-air mass flow understates the
% enthalpy by a factor (1+w). The check below evaluates the inlet stream
% both ways and reports the difference, so that any future divergence
% between the two ends of the air loop is detected at source rather than
% appearing indirectly as an unexplained term in the power balance.
% ===================================================================
% DATUM: T_coil, not Ta. The loop is closed, so the air arriving at the
% still came from the AWG coil, not from outside. See fan_energy().
w_in_chk    = P.wv_fan_in;
dT_in_chk   = P.T_air_in - P.T_coil;
Cp_da_chk   = 1005 + 0.05*(P.T_air_in - 273.15);
Cp_vap_chk  = 1860 + 0.12*(P.T_air_in - 273.15);
H_in_comp   = (P.mdot_da + P.mdot_da*w_in_chk*0) * Cp_da_chk * dT_in_chk ...
            + (P.mdot_da*w_in_chk) * Cp_vap_chk * dT_in_chk;      % per-component
H_in_moist  = (P.mdot_da*(1+w_in_chk)) * moist_air_cp(P.T_air_in, w_in_chk) * dT_in_chk;
H_in_wrong  = P.mdot_da * moist_air_cp(P.T_air_in, w_in_chk) * dT_in_chk;  % mismatched pairing

fprintf('\n  ---- air-stream enthalpy convention check ----\n');
fprintf('  inlet enthalpy, per-component     = %9.2f W\n', H_in_comp);
fprintf('  inlet enthalpy, moist-mixture     = %9.2f W\n', H_in_moist);
fprintf('  difference between conventions    = %9.3e W   %s\n', ...
    H_in_comp - H_in_moist, ...
    ternary(abs(H_in_comp - H_in_moist) < 1e-6, '[PASS]', '[FAIL]'));
fprintf('  (mismatched pairing would give %9.2f W, low by %.2f W)\n', ...
    H_in_wrong, H_in_comp - H_in_wrong);

fprintf('\n  ---- downstream chain (closed loop) ----\n');
if isfield(results,'awg') && isfield(results,'downstream')
    awg = results.awg;
    D   = results.downstream;

    fprintf('  [ELECTRIC HEATER]\n');
    fprintf('    brine in                        = %8.4f g/s at %.2f K, w_salt = %.4f\n', ...
        D.mdot_brine_in*1000, results.brine_out.T, D.w_salt_in);
    fprintf('    concentrated brine out          = %8.4f g/s, w_salt = %.4f\n', ...
        D.mdot_brine_out*1000, D.w_salt_out);
    fprintf('    vapour raised                   = %8.4f g/s at %.2f K\n', ...
        D.mdot_vap*1000, D.T_vap_heater);
    fprintf('    duty: sensible %.1f W + latent %.1f W = %.1f W\n', ...
        D.Q_heater_sensible, D.Q_heater_latent, D.Q_heater_useful);
    fprintf('    electrical draw                 = %8.2f W   (eta = %.2f, losses %.0f%%)\n', ...
        D.W_heater, results.P.bl.eta_heater, 100*results.P.bl.f_heatloss);

    fprintf('  [FEED PREHEAT HX]\n');
    fprintf('    duty absorbed by the feed       = %8.2f W   (%.2f K rise)\n', ...
        D.Q_feed_absorbed, D.dT_feed);
    fprintf('    feed  %.2f K -> %.2f K   (pinch ceiling %.2f K)\n', ...
        D.T_feed_in, D.T_feed_out, D.T_feed_cap);
    fprintf('    vapour desuperheat available    = %8.2f W\n', D.Q_desuperheat);
    fprintf('    vapour latent capacity          = %8.2f W\n', D.Q_cond_capacity);
    % THREE regimes, and the label comes from the solve. A two-regime
    % test would report PINCH-LIMITED whenever the unconstrained supply
    % exceeds the ceiling, including the usual case in which the DESIGN
    % TARGET is what sets the answer.
    fprintf('    condensed fraction f_cond       = %8.4f  (unclamped %.4f)   [%s-LIMITED]\n', ...
        D.f_cond, D.f_cond_raw, upper(D.binding_limit));
    fprintf('      FINAL INSTANT. PART A quotes the vapour-weighted window\n');
    fprintf('      mean (%.4f -), which is higher: the feed capacity rate is\n', ...
        pick(results,'wmean','f_cond', D.f_cond));
    fprintf('      fixed while vapour supply falls off with the sun, so the\n');
    fprintf('      feed condenses a LARGER share of a smaller stream as the\n');
    fprintf('      day ends. Both are correct on their own basis.\n');
    fprintf('    Surplus vapour  %8.4f g/s: %.4f g/s -> air preheat,\n', ...
        (D.mdot_cond_airheat + D.mdot_cond_awg_vapour)*1000, ...
        D.mdot_cond_airheat*1000);
    fprintf('                    %8.4f g/s -> AWG.\n', ...
        D.mdot_cond_awg_vapour*1000);
    fprintf('    Whatever the feed cannot absorb is REROUTED, not lost:\n');
    fprintf('    it condenses downstream and is collected as product.\n');
    fprintf('  [AIR-TO-AIR RECUPERATOR]  (passive, no power)\n');
    fprintf('    effectiveness eps               = %8.3f  -\n', D.eps_recup);
    fprintf('    hot  %.2f K -> %.2f K   (recovered %.2f K of %.2f K available)\n', ...
        awg.Tv_exit, D.T_recup_hot_out, D.dT_recup, D.dT_recup_max);
    fprintf('    cold %.2f K -> %.2f K   (%.2f K rise to the heater)\n', ...
        results.P.T_coil, D.T_air_recup_out, D.dT_air_recup);
    fprintf('    heat recovered: sensible %.1f W + latent %.1f W = %.1f W\n', ...
        D.Q_recup_sensible, D.Q_recup_latent, D.Q_recup_total);
    if D.mdot_cond_recup > 0
        fprintf('    hot side crossed its dew point: %.4f g/s drained here -> product\n', ...
            D.mdot_cond_recup*1000);
    else
        fprintf('    hot side stayed above its dew point: no condensate here\n');
    end
    fprintf('    air lift avoided                = %8.2f W of %.2f W gross\n', ...
        D.Q_air_gross - D.Q_air_full, D.Q_air_gross);
    fprintf('  [VAPOUR-FIRED AIR PREHEATER]  (duty-specified, T_int is the output)\n');
    % "of the lift" is ambiguous: f_air_recovered is referred to the
    % POST-RECUPERATOR lift (Q_air_full), while the thermal-input block
    % refers the same duty to the GROSS lift (Q_air_gross). The two
    % denominators give different percentages for one duty, so each is
    % named explicitly wherever it is printed.
    fprintf('    surplus steam recovered         = %8.2f W  (%.1f%% of post-recup\n', ...
        D.Q_air_cond, 100*D.f_air_recovered);
    fprintf('                                                lift, %.1f%% of gross lift)\n', ...
        100*D.Q_air_cond/max(D.Q_air_gross,eps));
    fprintf('    air  %.2f K -> %.2f K (T_int)   = %8.2f K rise\n', ...
        D.T_air_recup_out, D.T_int, D.dT_air_cond);
    fprintf('      inlet is the RECUPERATOR outlet, not the coil: this unit\n');
    fprintf('      is third in the air path (coil -> recuperator -> here).\n');
    fprintf('    condensate collected here       = %8.4f g/s  -> product\n', ...
        D.mdot_cond_airheat*1000);

    fprintf('  [ELECTRIC AIR HEATER]  (outlet-T-specified, duty is the output)\n');
    fprintf('    lift after recuperator %.2f K -> %.2f K = %8.2f W\n', ...
        D.T_air_recup_out, results.P.T_air_in, D.Q_air_full);
    fprintf('    trim duty  %.2f K -> %.2f K     = %8.2f W\n', ...
        D.T_int, results.P.T_air_in, D.Q_air_trim);
    fprintf('    electrical rating (eta = %.2f)   = %8.2f W\n', ...
        results.P.eta_air_heater, D.W_air_heater);
    fprintf('    NOTE: RATE the element at the FULL lift (%.2f W). The steam\n', D.Q_air_full);
    fprintf('          supply falls to zero in the supply-limited regime, so\n');
    fprintf('          the element must be able to carry the whole duty alone.\n');

    fprintf('  [AWG]  two inlet streams\n');
    fprintf('    A: surplus vapour from the HX   = %8.4f g/s  -> fully condensed\n', ...
        D.mdot_cond_awg_vapour*1000);
    fprintf('    B: cascade exhaust              = %8.2f K, w = %.5f kg/kg\n', ...
        awg.Tv_exit, awg.wv_exit);
    fprintf('       dew point of the exhaust     = %8.2f K\n', awg.T_dew_return);
    fprintf('       coil / reject state          = %8.2f K, w = %.5f kg/kg (saturated)\n', ...
        awg.T_coil, awg.wv_coil);
    fprintf('       condensate from the air      = %8.4f g/s\n', D.mdot_cond_air*1000);
    fprintf('    cooling duty: A %.1f W + B %.1f W = %.1f W\n', ...
        D.Q_awg_vapour, D.Q_awg_air, D.Q_awg_total);
    fprintf('    reject temperature              = %8.2f K   (Ta + %.1f K approach)\n', ...
        awg.T_cond, results.P.dT_cond_approach);
    fprintf('    heat-pump lift                  = %8.2f K\n', awg.lift);
    fprintf('    COP (eta_II = %.2f, lift-dep.)  = %8.2f\n', awg.eta_II, awg.COP);
    fprintf('    compressor work                 = %8.2f W\n', awg.W_compressor);

    fprintf('  [AIR-SIDE MOISTURE AUDIT]\n');
    fprintf('    evaporation in the cascade      = %8.4f g/s\n', awg.mdot_evap_still*1000);
    fprintf('    drained in the recuperator      = %8.4f g/s\n', D.mdot_cond_recup*1000);
    fprintf('    condensed on the coil           = %8.4f g/s\n', D.mdot_cond_air*1000);
    fprintf('    recovered from the air, TOTAL   = %8.4f g/s\n', ...
        (D.mdot_cond_recup + D.mdot_cond_air)*1000);
    fprintf('    rejected to atmosphere          = %8.4f g/s   (0: loop is closed)\n', ...
        awg.mdot_moisture_rejected*1000);
    fprintf('    humidity-ratio closure residual = %8.2e g/s   (must be ~0)\n', ...
        awg.loop_moisture_residual*1000);

    fprintf('\n  ---- energy crossing the modelled boundary (the still) ----\n');
    fprintf('  absorbed solar radiation          = %8.2f W\n', awg.Q_solar_absorbed);
    fprintf('  blower electrical work            = %8.2f W\n', awg.W_blower);
    fprintf('  total                             = %8.2f W\n', awg.E_still_inputs);
    fprintf('  The AWG electrical demand lies outside THIS boundary. Do not quote\n');
    fprintf('  the still-only SEC as a system figure -- see the SYSTEM ENERGY\n');
    fprintf('  block below, which encloses both components.\n');
end

fprintf('\n  ---- thermal input by source ----\n');
fprintf('  (air terms referenced to T_coil = %.2f K, the state at which\n', P.T_coil);
fprintf('   loop air leaves the AWG. Ambient is NOT a datum here: no\n');
fprintf('   ambient air enters the closed loop.)\n');
% The air heater's REAL duty, on the T_coil datum -- the energy actually
% injected into the loop. steady_check.Q_fan_thermal is the ambient-datum
% BALANCE term and is not this unit's duty; it is printed in the closure
% block above and must not be quoted as a load.
% DATUM. Q_air_heater is Q_air_full, the lift REMAINING after the
% recuperator, i.e. referred to T_air_recup_out. On the T_coil datum
% declared in the header the energy entering the loop air is
% Q_air_gross, with the recuperator as a third source alongside the steam
% and the element. Q_air_gross is therefore used here; Q_air_full under a
% T_coil label would omit the passive recovery from the numerator and
% inflate every share computed against it.
Q_ah_ext   = results.downstream.Q_air_gross;                       % [W] coil -> setpoint
Q_recup_in = results.downstream.Q_recup_total;                     % [W] passive, free
Q_src_ext  = max(steady_check.Q_solar_absorbed,0) + max(Q_ah_ext,0) ...
             + max(steady_check.Q_feed_thermal,0);

% All three rows divide by the SAME denominator, Q_src_ext. An inline
% denominator on the solar row alone (solar + Q_air_heater) would leave
% three rows on two denominators and the shares would not sum to 100%.
fprintf('  absorbed solar radiation          = %8.2f W   (%5.1f %%)\n', ...
    steady_check.Q_solar_absorbed, ...
    100*max(steady_check.Q_solar_absorbed,0)/max(Q_src_ext,eps));
% ONE denominator for all three rows, or the shares do not sum to 100%.
% TWO SOURCES supply the air lift, and only one of them is purchased.
% Q_air_heater is the TOTAL thermal delivered to the loop air, which is
% the right denominator term; but the steam-fired share is internal
% recycle (evaporator vapour raised from the plant's own brine, already
% charged to the electric heater in the SEC block) and must not be
% double-counted as purchased work.
fprintf('  air lift, TOTAL into loop air     = %8.2f W   (%5.1f %%)\n', ...
    Q_ah_ext, 100*max(Q_ah_ext,0)/max(Q_src_ext,eps));
fprintf('    (coil %.2f K -> still inlet %.2f K, ALL sources)\n', ...
    P.T_coil, P.T_air_in);
fprintf('    of which recuperated [PASSIVE]  = %8.2f W   (%5.1f %% of lift)\n', ...
    Q_recup_in, 100*max(Q_recup_in,0)/max(Q_ah_ext,eps));
fprintf('    of which steam-fired [RECYCLE]  = %8.2f W   (%5.1f %% of lift)\n', ...
    results.downstream.Q_air_cond, ...
    100*max(results.downstream.Q_air_cond,0)/max(Q_ah_ext,eps));
fprintf('    of which electric trim [PURCH.] = %8.2f W   (%5.1f %% of lift)\n', ...
    results.downstream.Q_air_trim, ...
    100*max(results.downstream.Q_air_trim,0)/max(Q_ah_ext,eps));
% The three shares must exhaust the lift; a gap means one of the units
% is being credited on a different datum from the others.
resid_lift = Q_ah_ext - (Q_recup_in + results.downstream.Q_air_cond ...
                         + results.downstream.Q_air_trim);
if abs(resid_lift) > 1e-6*max(Q_ah_ext,eps)
    fprintf('    ** LIFT SPLIT DOES NOT CLOSE: %+.3e W -- datum mismatch.\n', resid_lift);
end
fprintf('  feed preheat  [INTERNAL RECYCLE]  = %8.2f W   (%5.1f %%)\n', ...
    steady_check.Q_feed_thermal, ...
    100*max(steady_check.Q_feed_thermal,0)/max(Q_src_ext,eps));
fprintf('    Enthalpy of the feed AT THE CASCADE INLET, relative to\n');
fprintf('    ambient. It is a genuine boundary term of the still, but it\n');
fprintf('    is NOT purchased energy: the rise across the HX is\n');
fprintf('    evaporator vapour raised from the plant''s own brine and is\n');
fprintf('    already charged to the electric heater in the SEC block.\n');
% GOR DENOMINATOR: BOUNDARY SUPPLY, *NOT* Q_src_ext.
% Q_src_ext is a datum-shifted DESCRIPTION of the energy delivered to the
% loop air, and 36% of it is the passive recuperator -- heat recovered
% from the loop's own exhaust that never crossed the plant boundary.
% Putting it in a performance ratio is perverse: a BETTER recuperator
% would enlarge the denominator and lower the GOR, penalising the plant
% for recovering heat. GOR must divide by what the plant is SUPPLIED,
% which is Q_src_total (= Q_in in the closure block above).
% The decomposition rows above and this ratio therefore have different
% denominators BY DESIGN, and that is stated explicitly here.
fprintf('  GOR, absorbed-solar denominator   = %7.4f  -\n', steady_check.GOR_solar_only);
fprintf('  GOR, total thermal input          = %7.4f  -  (denominator: boundary\n', ...
    steady_check.GOR_total_input);
fprintf('                                                 supply %.2f W, NOT the\n', ...
    steady_check.Q_evap/max(steady_check.GOR_total_input,eps));
fprintf('                                                 %.2f W of the rows above --\n', Q_src_ext);
fprintf('                                                 those include passive\n');
fprintf('                                                 recuperation, which is not\n');
fprintf('                                                 a supplied input.)\n');
fprintf('    Quote the total-input GOR. The solar-only denominator is\n');
fprintf('    large only because absorbed solar is a small share of the\n');
fprintf('    thermal input; it is not a solar performance figure.\n');
% Use the SAME share printed above (Q_ah_ext basis), not the ambient-datum
% fan_share_pct, so that this sentence agrees with the table three lines
% earlier.
% ATTRIBUTION. Q_ah_ext is the WHOLE lift -- recuperator, steam and
% element together -- so it is not the element's own share: the
% resistance element supplies roughly half of the lift and a comparable
% fraction of the thermal input. The claim being made is that the input is
% PURCHASED rather than solar, so the electric trim is the right numerator.
share_ah    = 100*max(Q_ah_ext,0)/max(Q_src_ext,eps);
share_elec  = 100*max(results.downstream.Q_air_trim,0)/max(Q_src_ext,eps);
if share_ah > 50
    fprintf(['  The AIR LIFT supplies %.0f%% of the thermal input, of which the\n' ...
             '  ELECTRIC element carries %.0f%% of the total (the rest is passive\n' ...
             '  recuperation and steam recycle). Absorbed solar is %.1f%%. The\n' ...
             '  driver is purchased electricity, not insolation, so the\n' ...
             '  total-input basis is the only defensible performance measure.\n'], ...
             share_ah, share_elec, ...
             100*max(steady_check.Q_solar_absorbed,0)/max(Q_src_ext,eps));
end


end

function Td = dew_point_from_Pv(Pv)
% Dew-point temperature [K] from vapour partial pressure [Pa], by
% inversion of the Magnus-Tetens form used elsewhere in this model.
Pv = max(Pv, 1);
a  = 17.62; b = 243.12;              % Magnus coefficients, degC
g  = log(Pv/611.2);
Td = b*g ./ (a - g) + 273.15;
end % dew_point_from_Pv


function m = steady_check_mdot_evap_total(results)
% Total evaporation rate across the cascade at the final sample [kg/s].
P = results.P;
m = 0;
for k = 1:P.Ns
    m = m + sum(results.mevap(k,:,end)) * P.dx_stage(k) * P.W;
end
end % steady_check_mdot_evap_total


function Q_dot = instantaneous_storage_rate(results, P)
% Rate of change of stored sensible heat, summed over every control
% volume, evaluated at the final stored sample by backward difference.
%
%   Q_dot = sum_i ( m_i * Cp_i * dT_i/dt )    [W]
%
% Under time-varying irradiance this term is non-zero and must appear in
% the instantaneous power balance. Positive values indicate the device is
% accumulating sensible heat.
t = results.t;
if numel(t) < 2
    Q_dot = 0; return
end
% Time derivatives are taken as a least-squares slope over the last few
% stored samples rather than a two-point backward difference. Under a
% measured irradiance record the state is never stationary, so this term
% is genuinely non-zero and its accuracy limits how tightly the power
% balance can close; a two-point difference on a coarse output grid
% carries a first-order truncation error proportional to the local
% curvature of the forcing.
nfit = min(5, numel(t));
idxf = (numel(t)-nfit+1):numel(t);
tf   = t(idxf) - t(idxf(1));
if (t(end) - t(end-1)) <= 0
    Q_dot = 0; return
end
den  = sum((tf - mean(tf)).^2);
wfit = (tf - mean(tf)) / max(den, eps);        % slope operator: dphi/dt = sum(wfit.*phi)
slope3 = @(A) reshape(reshape(A(:,:,idxf), [], nfit) * wfit(:), size(A,1), size(A,2));
ddt = @(A) slope3(A);                          % 3-D fields (stage, node, time)

Q_dot = 0;

% ---- Glazing (P.mCp_g is areal) ----
Q_dot = Q_dot + P.mCp_g * P.Ag * (results.Tg(idxf(:)).' * wfit(:));

% ---- Water films: rho*Cp*delta integrated over each stage ----
dTw    = ddt(results.Tw);
delta_e = results.delta(:,:,end);
C_e     = results.C(:,:,end);
Tw_e    = results.Tw(:,:,end);
for k = 1:P.Ns
    rho_k = water_density(Tw_e(k,:), C_e(k,:));
    Cp_k  = water_cp(Tw_e(k,:), C_e(k,:), P);
    Q_dot = Q_dot + sum(rho_k .* Cp_k .* delta_e(k,:) .* dTw(k,:)) * P.dx_stage(k) * P.W;
end

% ---- Cascade plates and floor slab ----
dTp = ddt(results.Tp);
for k = 1:P.Ns
    mCp_areal = P.rho_solid_stage(k) * P.Cp_solid_stage(k) * P.t_solid_stage(k);
    Q_dot = Q_dot + mCp_areal * sum(dTp(k,:)) * P.dx_stage(k) * P.W;
end

% ---- Wall / insulation nodes ----
Awall_vec = wall_area_per_gap(P);
dTwall = (results.Twall(idxf,:).' * wfit(:));
Q_dot = Q_dot + sum(P.mCp_wall_areal .* Awall_vec(:) .* dTwall(:));

% ---- Vapour gaps (sensible content of the moist air) ----
dTv = (results.Tv(idxf,:).' * wfit(:));
V_gap = P.V_gap_vec;
for j = 1:P.Ns
    rho_j = moist_air_density(results.Tv(end,j), results.wv(end,j));
    Cp_j  = moist_air_cp(results.Tv(end,j), results.wv(end,j));
    Q_dot = Q_dot + rho_j * V_gap(j) * Cp_j * dTv(j);
end
end % instantaneous_storage_rate


function run_diagnostics(results)
P = results.P;
issues = {};

hdr('STATE DIAGNOSTICS', 'Physical plausibility of the solved fields.');

fields_to_check = {'Tg','Tv','wv','Twall','Tw','delta','C','Tp','u','mevap','mfw'};
for i = 1:numel(fields_to_check)
    fname = fields_to_check{i};
    val = results.(fname);
    n_nan = sum(isnan(val(:)));
    n_inf = sum(isinf(val(:)));
    if n_nan > 0 || n_inf > 0
        issues{end+1} = sprintf('FAIL: %s contains %d NaN and %d Inf values -- integration broke down.', ...
            fname, n_nan, n_inf); %#ok<AGROW>
    end
end
if isempty(issues)
    fprintf('[PASS] No NaN/Inf values found in any field.\n');
else
    for i = 1:numel(issues)
        fprintf('%s\n', issues{i});
    end
end

delta_final = results.delta(:,:,end);
min_delta = min(delta_final(:));
[k_dry, x_dry] = find(delta_final == min_delta, 1);
if min_delta < P.delta_dryout
    fprintf('[FAIL] Film dry-out: minimum thickness %.4g mm at Plate %d, node %d (final time)\n', ...
        min_delta*1000, k_dry, x_dry);
    fprintf('       is AT/BELOW the dry-out floor (%.4g mm) -- evaporation there is fully throttled.\n', ...
        P.delta_dryout*1000);
elseif min_delta < P.delta_dryout_full
    pct = 100*(min_delta - P.delta_dryout)/(P.delta_dryout_full - P.delta_dryout);
    fprintf('[WARN] Film approaching dry-out: minimum thickness %.4g mm at Plate %d, node %d\n', ...
        min_delta*1000, k_dry, x_dry);
    fprintf('       is only %.0f%% of the way from dry-out (%.4g mm) to fully-unthrottled (%.4g mm).\n', ...
        pct, P.delta_dryout*1000, P.delta_dryout_full*1000);
else
    fprintf('[PASS] Film thickness stays above the dry-out throttle range everywhere (min = %.4g mm).\n', ...
        min_delta*1000);
end

C_final = results.C(:,:,end);
max_C_all_time = squeeze(max(max(results.C,[],1),[],2));
[max_C, i_worst] = max(max_C_all_time);
if max_C >= P.C_saturation
    fprintf('[FAIL] TDS reached %.1f kg/m3 at t = %.2f s (cap = %.1f)\n', ...
        max_C, results.t(i_worst), P.C_saturation);
elseif max_C >= P.C_saturation_ramp_start
    fprintf('[WARN] TDS approached %.1f kg/m3 at t = %.2f s (ramp-start = %.1f)\n', ...
        max_C, results.t(i_worst), P.C_saturation_ramp_start);
end

Tw_final = results.Tw(:,:,end);
rho_check =water_density(Tw_final, C_final);
mu_check  = water_viscosity(Tw_final, C_final);
u_final   = results.u(:,:,end);

bad_rho = (rho_check < 900) | (rho_check > 1500) | isnan(rho_check);
bad_mu  = (mu_check < 1e-5) | (mu_check > 1e-1) | isnan(mu_check);
bad_u   = (u_final < 0) | (u_final > 5) | isnan(u_final);
bad_T   = (Tw_final < 273) | (Tw_final > 373);

if any(bad_rho(:))
    fprintf('[FAIL] Water density outside plausible range (900-1500 kg/m3) at %d node(s).\n', sum(bad_rho(:)));
else
    fprintf('[PASS] Water density stays within a physically plausible range.\n');
end
if any(bad_mu(:))
    fprintf('[FAIL] Water viscosity outside plausible range (1e-5 to 1e-1 Pa.s) at %d node(s).\n', sum(bad_mu(:)));
else
    fprintf('[PASS] Water viscosity stays within a physically plausible range.\n');
end
if any(bad_u(:))
    fprintf('[WARN] Film velocity outside plausible range (0-5 m/s) at %d node(s) -- check for a units/formula error.\n', sum(bad_u(:)));
else
    fprintf('[PASS] Film velocity stays within a physically plausible range.\n');
end
if any(bad_T(:))
    fprintf('[FAIL] Film temperature outside liquid-water range (0-100 C) at %d node(s).\n', sum(bad_T(:)));
else
    fprintf('[PASS] Film temperature stays within the liquid-water range.\n');
end

rt = results.retention_time_plate;
if any(rt <= 0) || any(isnan(rt)) || any(isinf(rt))
    fprintf('[FAIL] Retention time is non-positive/NaN/Inf for at least one plate -- check film velocity.\n');
elseif any(rt > 3600)
    fprintf('[WARN] At least one plate has retention time over 1 hour (%.1f s) -- film may be nearly stalled.\n', max(rt));
else
    fprintf('[PASS] Retention times are all positive and within a plausible range.\n');
end

foot();

end % run_diagnostics

function print_vapor_pressure_diagnostics(results)
% Per-plate vapor-pressure / evaporation-driving-force table at final
% time. Recomputes Psat_w, Pv, dP, hm, hfg, mevap fresh from the state
% (same formulas used inside rhs()/extract_results()), so the numbers
% here are guaranteed consistent with what the solver actually used --
% nothing here is a separate/parallel calculation path.

P  = results.P;
iE = numel(results.t);

S.Tg    = results.Tg(iE);
S.Tv    = results.Tv(iE,:).';
S.wv    = results.wv(iE,:).';
S.Tw    = results.Tw(:,:,iE);
S.delta = results.delta(:,:,iE);
S.C     = results.C(:,:,iE);
S.Tp    = results.Tp(:,:,iE);

hdr('EVAPORATION DRIVING FORCE, PER STAGE   (final time)', ...
    'Saturation and bulk vapour pressures, transfer coefficient, latent heat.');

fprintf(' k  |  Tw_avg  C_avg  |   Psat_w      Pv       dP    |   hm        hfg      |  mevap    Q_evap\n');
fprintf('    |  [K]    [kg/m3]|   [Pa]        [Pa]     [Pa]   |  [m/s]    [J/kg]     |  [g/s]     [W]\n');
prule();

Q_evap_check = 0;
mdot_check   = 0;

% ---- Exact local mass-flow vector at final time ----
mdot_local_vec_iE = P.mdot_da + P.mdot_vapor_in + ...
    flipud(cumsum(flipud(results.mevap_plate(iE,:).')));

for k = 1:P.Ns
    Tw_k = S.Tw(k,:);
    C_k  = S.C(k,:);
    Tv_k = S.Tv(k);
    wv_k = S.wv(k);

    [~, hm_wv] = film_gap_coeffs(k, S, P, mdot_local_vec_iE(k));
    Psat_w = psat_saline(Tw_k, C_k, P);
    Pv_k   = vapor_partial_pressure(wv_k, Tv_k, P);
    dP     = evap_driving_dP(Psat_w, Pv_k, P);
    mevap_k = hm_wv .* dP .* P.Mwater ./ (P.Rgas .* Tw_k);
    mdot_k  = sum(mevap_k) * P.dx_stage(k) * P.W;         % [kg/s], this stage
    hfg_k   = mean(latent_heat(Tw_k, C_k, P));    % LOCAL hfg, this stage
    Q_evap_k = mdot_k * hfg_k;                    % [W], this plate

    Q_evap_check = Q_evap_check + Q_evap_k;
    mdot_check   = mdot_check + mdot_k;

    fprintf(' %2d  | %6.2f  %6.2f | %8.2f  %8.2f  %6.2f | %.4e  %8.1f    | %7.4f  %8.3f\n', ...
        k, mean(Tw_k), mean(C_k), ...
        mean(Psat_w), Pv_k, mean(dP), ...
        mean(hm_wv), hfg_k, ...
        mdot_k*1000, Q_evap_k);
end
foot();
end % print_vapor_pressure_diagnostics

function Print_AWG_Duty(results)
       
    a = results.awg;
    hdr('AWG DUTY   (from the still)', 'Humid-air state delivered and the duty required.');
    fprintf('  Humid air delivered to the AWG:\n');
    fprintf('    dry-air mass flow             = %8.4f kg/s\n', results.P.mdot_da);
    fprintf('    temperature                   = %8.2f K\n',    a.Tv_exit);
    fprintf('    humidity ratio                = %8.5f kg/kg\n', a.wv_exit);
    fprintf('    dew point                     = %8.2f K\n',    a.T_dew_return);
    fprintf('  Required of the AWG:\n');
    fprintf('    coil temperature (assumed)    = %8.2f K\n',    a.T_coil);
    fprintf('    condensate on the COIL        = %8.4f g/s   (= %.2f L/day over %g h)\n', ...
        a.mdot_condensate*1000, a.mdot_condensate*results.P.t_operating*3600, results.P.t_operating);
    fprintf('    total cooling duty            = %8.2f W   (both inlet streams)\n', a.Q_evaporator);
    fprintf('    NOTE: the unit rejects to AMBIENT. There is no reheat duty in\n');
    fprintf('          the closed loop -- no air, and no moisture, leaves it.\n');
    fprintf('  Air-side moisture audit:\n');
    fprintf('    evaporation in the cascade    = %8.4f g/s\n', a.mdot_evap_still*1000);
    fprintf('    drained in the recuperator    = %8.4f g/s\n', ...
        results.downstream.mdot_cond_recup*1000);
    fprintf('    condensed on the coil         = %8.4f g/s\n', a.mdot_condensate*1000);
    fprintf('    recovered from the air, TOTAL = %8.4f g/s\n', ...
        (results.downstream.mdot_cond_recup + a.mdot_condensate)*1000);
    fprintf('    rejected to atmosphere        = %8.4f g/s\n', a.mdot_moisture_rejected*1000);
    fprintf('    closure residual              = %8.2e g/s\n', a.loop_moisture_residual*1000);
    foot();

end

function print_procurement_spec(results)
% PART C -- duty specifications for procurement, one block per unit in
% process order.
%
% Basis is the final instant: equipment is rated for the duty it carries at
% the close of the window, not for the daily average. Condensate yields in
% the plasma block and the plant water balance are window means and are
% tagged as such. Stream identifiers refer to PART B.
P      = results.P;
t_op_s = P.t_operating*3600;
a      = results.awg;
b      = results.brine_out;
D      = results.downstream;
kgday  = @(m) m*t_op_s;

m_evap_still  = a.mdot_evap_still;                                 % [kg/s]
m_atm_harvest = P.mdot_da*(P.wv_fan_in - a.wv_coil);               % [kg/s]
m_cond_air    = a.mdot_condensate;                                 % [kg/s]
resid_air     = (m_cond_air + D.mdot_cond_recup) - (m_evap_still + m_atm_harvest);

hdr('PART C -- UNIT SPECIFICATIONS', ...
    'Basis: final instant unless tagged [win]. Stream IDs refer to PART B.');

% ------------------------------------------------------- cascade discharge
sec('CASCADE BRINE OUTLET   (stream, not a purchase)');
prow('Brine mass flow',      '%.4f', b.mdot*1000,        'g/s',   sprintf('%.2f kg/day', kgday(b.mdot)));
prow('Brine TDS',            '%.1f', b.C,                'kg/m3', '');
prow('Concentration factor', '%.2f', b.CF,               '-',     '');
prow('Brine temperature',    '%.2f', b.T,                'K',     '');
prow('Salt mass fraction',   '%.4f', b.w_salt,           '-',     '');
note('This stream is the crystalliser feed. ZLD closure occurs there, not');
note('in the still, and the cascade does not reach salt saturation.');

% ---------------------------------------------------------- crystalliser
sec(sprintf('CRYSTALLISER / ELECTRIC HEATER   (to %.1f wt%% salt)', 100*D.w_salt_out));
prow('Brine inlet, S4',      '%.4f', D.mdot_brine_in*1000,'g/s',  sprintf('%.2f kg/day at %.2f K', kgday(D.mdot_brine_in), b.T));
prow('  of which water',     '%.4f', D.mdot_water_in*1000,'g/s',  sprintf('%.2f kg/day', kgday(D.mdot_water_in)));
prow('Inlet salt fraction',  '%.4f', D.w_salt_in,        '-',     'TDS / rho');
prow('Salt, conserved',      '%.4f', D.mdot_salt*1000,   'g/s',   sprintf('%.2f kg/day', kgday(D.mdot_salt)));
prow('Salt slurry outlet',   '%.4f', D.mdot_brine_out*1000,'g/s', sprintf('%.2f kg/day', kgday(D.mdot_brine_out)));
prow('  residual water',     '%.4f', (D.mdot_brine_out-D.mdot_salt)*1000,'g/s', ...
     sprintf('%.2f kg/day', kgday(D.mdot_brine_out-D.mdot_salt)));
prow('Vapour raised, S5',    '%.4f', D.mdot_vap*1000,    'g/s',   sprintf('%.2f kg/day at %.2f K', kgday(D.mdot_vap), D.T_vap_heater));
prow('Mass closure',         '%.2e', (D.mdot_brine_in-D.mdot_vap-D.mdot_brine_out)*1000,'g/s','in - vap - slurry');
prow('Water closure',        '%.2e', D.water_residual*1000,'g/s', 'water_in - vap - slurry water');
prow('Sensible duty',        '%.3f', D.Q_heater_sensible/1000,'kW','');
prow('Latent duty',          '%.3f', D.Q_heater_latent/1000,'kW', '');
prow('DUTY',                 '%.3f', D.Q_heater_useful/1000,'kW', 'quote against this');
prow('Electrical rating',    '%.3f', D.W_heater/1000,    'kW',    'includes vessel losses');

% -------------------------------------------------------------- feed HX
sec('FEED PREHEAT HX   (crystalliser vapour -> seawater feed)');
prow('Feed inlet',           '%.2f', D.T_feed_in,        'K',     'supplied seawater');
prow('Feed outlet',          '%.2f', D.T_feed_out,       'K',     'solved in the RHS');
prow('Feed capacity rate',   '%.3f', D.Cden_feed,        'W/K',   '');
prow('Duty absorbed',        '%.3f', D.Q_feed_absorbed/1000,'kW', '');
prow('Vapour inlet',         '%.4f', D.mdot_vap*1000,    'g/s',   sprintf('at %.2f K liquor', D.T_boil));
prow('Condensing at',        '%.2f', D.T_sat_vap,        'K',     sprintf('BPE %.1f K stays upstream', D.BPE));
prow('Pinch ceiling',        '%.2f', D.T_feed_cap,       'K',     '');
prow('Desuperheat duty',     '%.3f', D.Q_desuperheat/1000,'kW',   '');
prow('Latent capacity',      '%.3f', D.Q_cond_capacity/1000,'kW', '');
prow('Condensed fraction',   '%.4f', D.f_cond,           '-',     sprintf('%s-limited', D.binding_limit));
note(sprintf('Vapour split, %.4f g/s total:', D.mdot_vap*1000));
prow('  condensed here',     '%.4f', D.mdot_cond_preheater*1000,'g/s', sprintf('%.2f kg/day, product', kgday(D.mdot_cond_preheater)));
prow('  to air preheat',     '%.4f', D.mdot_cond_airheat*1000,'g/s',   sprintf('%.2f kg/day, product', kgday(D.mdot_cond_airheat)));
prow('  to AWG',             '%.4f', D.mdot_cond_awg_vapour*1000,'g/s',sprintf('%.2f kg/day, stream A', kgday(D.mdot_cond_awg_vapour)));

% ------------------------------------------------------------------ AWG
sec('AWG / DEHUMIDIFIER   (two inlet streams)');
prow('Air volumetric flow',  '%.3f', P.Qfan,             'm3/s',  '');
prow('Dry-air mass flow',    '%.4f', P.mdot_da,          'kg/s',  '');
prow('Inlet air temperature','%.2f', a.Tv_exit,          'K',     'cascade exhaust');
prow('Inlet humidity ratio', '%.5f', a.wv_exit,          'kg/kg', '');
prow('Loop return humidity', '%.5f', P.wv_fan_in,        'kg/kg', '= w_sat(T_coil)');
prow('Coil saturation ratio','%.5f', a.wv_coil,          'kg/kg', '');
prow('Inlet dew point',      '%.2f', a.T_dew_return,     'K',     '');
prow('Coil temperature',     '%.2f', a.T_coil,           'K',     '');
prow('Dew-point margin',     '%.2f', a.T_dew_return-a.T_coil,'K', ...
     ternary(a.T_dew_return-a.T_coil > 3, 'ok', 'BINDING'));
prow('COOLING DUTY',         '%.3f', a.Q_evaporator/1000,'kW',    'quote against this');
prow('Compressor power',     '%.3f', a.W_compressor/1000,'kW',    sprintf('eta_II %.2f', a.eta_II));

sub('stream A, pure vapour from the feed HX');
prow('Vapour mass flow',     '%.4f', D.mdot_cond_awg_vapour*1000,'g/s', sprintf('%.2f kg/day', kgday(D.mdot_cond_awg_vapour)));
prow('Duty on this stream',  '%.3f', D.Q_awg_vapour/1000,'kW',    'condenses fully');

sub('stream B, humid cascade exhaust');
prow('Condensate recovered', '%.4f', m_cond_air*1000,    'g/s',   sprintf('%.2f kg/day', kgday(m_cond_air)));
prow('Duty on this stream',  '%.3f', D.Q_awg_air/1000,   'kW',    '');
if m_cond_air <= eps
    note(sprintf('No recovery: the coil at %.2f K sits at or above the exhaust', a.T_coil));
    note(sprintf('dew point of %.2f K, so nothing condenses on this stream.', a.T_dew_return));
else
    prow('  S3 evaporation',     '%.4f', m_evap_still*1000,'g/s', sprintf('%.2f kg/day', kgday(m_evap_still)));
    prow('  S7c recuperator drain','%.4f',-D.mdot_cond_recup*1000,'g/s', sprintf('-%.2f kg/day', kgday(D.mdot_cond_recup)));
    prow('  closure',            '%.2e', resid_air*1000,'g/s',    ternary(abs(resid_air) < 1e-6,'ok','CHECK'));
end
note('The reject goes to ambient; reheat comes from the separate air heater.');

% --------------------------------------------------------------- plasma
m_prod_total = results.product.mdot_total;
sec('PLASMA TREATMENT CHAMBER   (specified by throughput)');
prow('Feed-HX condensate',   '%.4f', results.product.mdot_preheater*1000,'g/s', sprintf('%.2f kg/day', kgday(results.product.mdot_preheater)));
prow('Air-preheat condensate','%.4f',results.product.mdot_airheat*1000,'g/s',   sprintf('%.2f kg/day', kgday(results.product.mdot_airheat)));
prow('Recuperator drain',    '%.4f', results.product.mdot_recup*1000,'g/s',     sprintf('%.2f kg/day', kgday(results.product.mdot_recup)));
prow('AWG stream A',         '%.4f', results.product.mdot_awg_vapour*1000,'g/s',sprintf('%.2f kg/day', kgday(results.product.mdot_awg_vapour)));
prow('AWG stream B',         '%.4f', results.product.mdot_awg_air*1000,'g/s',   sprintf('%.2f kg/day', kgday(results.product.mdot_awg_air)));
prule();
prow('TOTAL FRESHWATER',     '%.4f', m_prod_total*1000,  'g/s',   sprintf('%.2f kg/day, %.2f L/h', kgday(m_prod_total), m_prod_total*3600));
prow('Feed temperature',     '%.2f', a.T_coil,           'K',     'approximate');
note('DBD units are quoted on throughput. No plasma chemistry is modelled');
note('and none is required for sizing.');

% ------------------------------------------------------- electrical load
sec('ELECTRICAL LOAD   (for PV and battery sizing)');
note('Rate the hardware on these figures. The window-mean loads behind SEC');
note('in PART A are lower and would under-rate the equipment.');
prow('Blower',               '%.3f', a.W_blower/1000,    'kW',    '');
prow('Crystalliser heater',  '%.3f', D.W_heater/1000,    'kW',    '');
prow('AWG compressor',       '%.3f', a.W_compressor/1000,'kW',    '');
prow('Electric air heater',  '%.3f', D.W_air_heater/1000,'kW',    '');
prow('Recuperator',          '%.3f', 0,                  'kW',    'passive');
prule();
prow('TOTAL',                '%.3f', (a.W_blower+D.W_heater+a.W_compressor+D.W_air_heater)/1000,'kW','');

% ---------------------------------------------------- plant water balance
m_feed_water  = P.Vfeed*(P.rhow_in - P.TDSfeed)/1000;              % [kg/day]
if isfield(results,'wmean') && isfield(results.wmean,'valid') && results.wmean.valid
    w_slurry_day = kgday(results.wmean.mdot_brine_out - results.wmean.mdot_salt);
else
    w_slurry_day = kgday(D.mdot_brine_out - D.mdot_salt);
end
prod_day     = kgday(m_prod_total);
res_seawater = m_feed_water - prod_day - w_slurry_day;

sec(sprintf('PLANT WATER BALANCE   (over %.1f h of operation)', P.t_operating));
prow('Feed water in',        '%.2f', m_feed_water,       'kg/day','[win]');
prow('Freshwater product',   '%.2f', prod_day,           'kg/day',sprintf('[win] %s', pickstr(results,'product','basis','final instant')));
prow('Water in the slurry',  '%.2f', w_slurry_day,       'kg/day','[win]');
prule();
prow('SEAWATER RECOVERY',    '%.2f', 100*prod_day/max(m_feed_water,eps),'%','[win] quote this figure');
prow('Residual',             '%.3f', res_seawater,       'kg/day',sprintf('%.2f %% of feed water', 100*abs(res_seawater)/max(m_feed_water,eps)));
note('Two terms make up the residual: the product is a window integral');
note('against a final-state slurry rate, and film holdup is still filling at');
note('t_end. Both are real; a growing residual points at the cascade-to-');
note('crystalliser handoff.');
foot();
end


function v = pickstr(R, grp, fld, fallback)
% String-valued twin of pick(), for basis labels.
v = fallback;
if isfield(R, grp) && isfield(R.(grp), fld) && ischar(R.(grp).(fld))
    v = R.(grp).(fld);
end
end

function v = pick(R, grp, fld, fallback)
% Return R.(grp).(fld) when the group exists and is flagged valid,
% otherwise the fallback. Used so that a figure printed in several
% places resolves to ONE object rather than to whichever local variable
% happened to be in scope at each print site.
v = fallback;
if isfield(R, grp)
    g = R.(grp);
    if isfield(g,'valid') && ~g.valid, return; end
    if isfield(g, fld) && isfinite(g.(fld)), v = g.(fld); end
end
end

function s = ternary(cond, a, b)
% Small helper: returns a if cond is true, else b (MATLAB has no
% built-in ternary/inline-if operator).
if cond, s = a; else, s = b; end
end


% =========================================================================
function print_system_energy(results)
% Specific energy on the plant control volume: still, crystalliser, feed HX
% and AWG together, with only purchased work and absorbed solar crossing the
% boundary. The cascade-only SEC in PART A sits on a smaller control volume
% and the two must not be added.
if ~isfield(results,'system'), return; end
sysE = results.system;
a  = results.awg;
pr = results.product;

hdr('PLANT ENERGY   (still + crystalliser + feed HX + AWG)', ...
    'Crossing the boundary: blower, heater, compressor, air heater, solar.');

sec('PRODUCT AND WORK');
prow('Freshwater collected',  '%.4f', pr.mdot_total*1000, 'g/s',   sprintf('%.2f kg/day', pr.kg_per_day));
prow('Cascade evaporation',   '%.4f', results.awg.mdot_evap_still*1000,'g/s','not product until condensed');
prule();
prow('Blower',                '%.2f', sysE.W_blower,         'W',     '');
prow('Crystalliser heater',   '%.2f', sysE.W_heater,         'W',     '');
prow('AWG compressor',        '%.2f', sysE.W_compressor,     'W',     sprintf('COP %.2f, lift %.1f K', a.COP, a.lift));
prow('Electric air heater',   '%.2f', sysE.W_air_heater,     'W',     '');
prule();
prow('Total purchased work',  '%.2f', sysE.W_external,       'W',     '');
prow('Absorbed solar, free',  '%.2f', sysE.Q_solar_absorbed, 'W',     '');
prow('Solar share of input',  '%.2f', 100*sysE.Q_solar_absorbed/max(sysE.Q_solar_absorbed+sysE.W_external,eps), ...
                                                          '%',     '');

sec('SPECIFIC ENERGY');
prow('Blower',                '%.2f', sysE.SEC_blower,       'kWh/m3','');
prow('Crystalliser heater',   '%.2f', sysE.SEC_heater,       'kWh/m3','');
prow('AWG compressor',        '%.2f', sysE.SEC_compressor,   'kWh/m3','');
prow('Electric air heater',   '%.2f', sysE.SEC_air_heater,   'kWh/m3','');
prule();
prow('SEC, PLANT',            '%.2f', sysE.SEC_external,     'kWh/m3','quote this figure');
prow('  band, eta_II 0.55',   '%.2f', min(sysE.SEC_external_band),'kWh/m3', ...
     sprintf('to %.2f kWh/m3 at eta_II 0.35', max(sysE.SEC_external_band)));
foot();
end

% =========================================================================
function check_baseline_validity(results)
% Hard screens a reported design point must satisfy. Each line carries a
% verdict, the measured value with its unit, and the limit it was tested
% against. Explanatory text prints only where a screen fails, so a clean run
% reports one line per screen.
P = results.P;
hdr('VALIDITY SCREENS', 'Verdict, measured value, limit.');
fprintf('\n');

% ---- run covered the operating window -------------------------------
t_end = results.t(end);
ok    = t_end >= 0.99*P.t_operating*3600;
scr(ok, 'Run completion',        '%.0f', t_end,        's', ...
        sprintf('of %.0f s required', P.t_operating*3600));
if ~ok
    note('Stopped early on a terminal event. mfw_daily has fallen back to');
    note('rate extrapolation and is not a design figure.');
end
note(sprintf('mfw_daily basis: %s', results.mfw_daily_basis));

% ---- cascade evaporation cannot exceed the water fed ----------------
m_water_feed = P.Vfeed*(P.rhow_in - P.TDSfeed)/1000;               % [kg/day]
ok = results.water_recovery_pct <= 100;
scr(ok, 'Cascade evaporation',   '%.2f', results.water_recovery_pct, '%', ...
        sprintf('limit 100 %%, %.2f of %.2f kg/day fed', results.mfw_daily, m_water_feed));
if ~ok
    note('Impossible. The excess is initial film holdup draining, not');
    note('product: the storage rate is negative in the mass-balance block.');
end

% ---- plant product against the slurry ceiling -----------------------
pr_v = results.product;
if isfield(pr_v,'int_valid') && pr_v.int_valid
    m_prod_chk = pr_v.kg_per_day_int;  basis_chk = '[win]';
else
    m_prod_chk = pr_v.kg_per_day;      basis_chk = '[end x t_op]';
end
m_salt_chk  = P.Vfeed*P.TDSfeed/1000;
m_resid_chk = (m_salt_chk/P.bl.w_salt_target)*(1 - P.bl.w_salt_target);
rec_ceil    = 100*(m_water_feed - m_resid_chk)/max(m_water_feed,eps);
rec_prod    = 100*m_prod_chk/max(m_water_feed,eps);
ok = rec_prod <= rec_ceil + 1e-6;
scr(ok, 'Plant product',         '%.2f', rec_prod,     '%', ...
        sprintf('ceiling %.2f %%, %s', rec_ceil, basis_chk));
if ~ok
    note('Product exceeds what the feed can supply. The loop admits no');
    note('make-up, so there is no second water source: this is a basis or');
    note('storage-drainage artefact, not yield.');
end

% ---- TDS excursion over the whole trajectory ------------------------
% Tested on the unclamped field. results.C_out is capped at P.C_saturation,
% so a test against it saturates at the ceiling and reports the ceiling.
% A mid-run excursion invalidates the run even if the endpoint relaxes,
% because the model carries no crystallisation closure.
if isfield(results,'C_raw')
    C_traj_max = max(max(results.C_raw(:,end,:),[],1),[],3);
    ok = C_traj_max <= P.C_saturation;
    scr(ok, 'TDS trajectory max', '%.1f', C_traj_max,  'kg/m3', ...
            sprintf('limit %.1f kg/m3, final state %.1f kg/m3', P.C_saturation, results.brine_out.C));
    if ~ok
        note('The run passed through a state the model cannot represent.');
        note('A relaxed endpoint does not repair it.');
    end
end

% ---- terminal brine below NaCl saturation ---------------------------
if isfield(results,'C_raw')
    TDS_out = results.C_raw(end,end,end);
else
    TDS_out = results.TDS_out_per_plate(end);
end
ok = TDS_out <= P.C_saturation;
scr(ok, 'Brine concentration',   '%.1f', TDS_out,      'kg/m3', ...
        sprintf('limit %.1f kg/m3, solution basis', P.C_saturation));
if ~ok
    note('Above saturation the film concentration is clamped for property');
    note('evaluation only. Salt is conserved, but properties are evaluated');
    note('at a concentration the film does not have and no crystallisation');
    note('is modelled. The run is out of domain.');
end

% ---- coil below the dew point AT THE COIL INLET ---------------------
% Measured after the recuperator, not at the cascade exhaust. If the hot
% side crosses its dew point in the recuperator the air arrives saturated,
% a materially smaller margin than the raw exhaust would suggest.
if isfield(results,'awg')
    if isfield(results.downstream,'mdot_cond_recup') && results.downstream.mdot_cond_recup > 0
        T_dew_coil_in = results.downstream.T_recup_hot_out;
    else
        T_dew_coil_in = results.awg.T_dew_return;
    end
    marg = T_dew_coil_in - results.awg.T_coil;
    if isfield(P,'coil_mode') && strcmp(P.coil_mode,'auto') ...
            && ~(isfield(P,'T_coil_frosted') && P.T_coil_frosted)
        info('Dew-point margin', '%.2f', marg, 'K', 'auto coil, identity not a test');
    else
        ok = marg > 3;
        scr(ok, 'Dew-point margin',  '%.2f', marg,     'K',  'limit 3 K, at the coil inlet');
        if ~ok
            note('The coil sits too close to the dew point, so the AWG recovers');
            note('little from the cascade air. In a closed loop the remedy is not');
            note('coil_mode = auto, which is rejected because T_coil sets the still');
            note('inlet humidity and deriving it from the solved exhaust is');
            note('circular. Lower FC.T_coil, or raise FC.T_air_in.');
        end
    end
    ok = abs(results.awg.loop_moisture_residual) < 1e-5;
    scr(ok, 'Air-side moisture closure','%.2e', results.awg.loop_moisture_residual,'kg/s','limit 1e-5 kg/s');
end

% ---- feed preheat HX ------------------------------------------------
if isfield(results,'downstream')
    D = results.downstream;
    if isfield(P,'tfeed_dynamic') && P.tfeed_dynamic
        % The exchanger cannot be infeasible in this configuration: the pinch
        % caps the feed temperature and surplus vapour is rerouted rather than
        % demanded and not supplied. Which limit binds is what matters.
        if isfield(D,'binding_limit')
            info('Feed HX limit', '%.2f', D.dT_preheat_actual, 'K', ...
                 sprintf('%s-limited, %.2f K available', D.binding_limit, ...
                         min(D.dT_preheat_max, D.dT_pinch_max)));
            note(sprintf('Surplus %.4f g/s to air preheat, %.4f g/s to AWG. Rerouted,', ...
                 D.mdot_cond_airheat*1000, D.mdot_cond_awg_vapour*1000));
            note('not lost: both condense as product.');
        end
    else
        ok = D.preheater_feasible;
        scr(ok, 'Preheater feasibility','%.3f', D.f_cond_raw,'-','limit 1');
        if ~ok
            note(sprintf('Tfan_in = %.1f K is unattainable from the vapour at', P.Tfan_in));
            note(sprintf('Qfan = %.2f m3/s and Ta = %.1f K. Shortfall %.0f W.', ...
                 P.Qfan, P.Ta, D.Q_preheater_deficit));
        end
    end
    ok = D.w_salt_in < P.bl.w_salt_target;
    scr(ok, 'Crystalliser step',     '%.4f', D.w_salt_in,'-', ...
            sprintf('to %.4f salt mass fraction', D.w_salt_out));
    if ~ok
        note('Brine is already at or above the target. The crystalliser has no');
        note('duty and the vapour stream is empty.');
    end
end

% ---- buildability and residence time --------------------------------
ok = P.stack_back <= 2.5;
scr(ok, 'Buildability',          '%.3f', P.stack_back,'m',  'limit 2.5 m, rear stack');

rt_hr = results.retention_time_total/3600;
ok    = rt_hr <= P.t_operating;
scr(ok, 'Residence time',        '%.2f', rt_hr,       'h',  sprintf('window %.2f h', P.t_operating));
if ~ok
    note('Deep stages never reach steady state within the day; they behave');
    note('as batch concentrators rather than as a cascade.');
end

% ---- terminal salt field, reported not screened ---------------------
% The cascade has no process TDS limit: whatever leaves the last stage goes
% to the crystalliser. The terminal concentration is an output to report,
% not a quantity that must converge before the flowsheet is valid, and a
% diurnal plant is not expected to reach a stationary salt distribution
% within one operating window. Computed on the unclamped field where
% available, since differencing two clamped values measures nothing.
if isfield(results.brine_out,'C_win')
    if isfield(results.brine_out,'C_raw') && isfield(results.brine_out,'C_win_raw')
        drift = 100*abs(results.brine_out.C_raw - results.brine_out.C_win_raw) ...
                / max(results.brine_out.C_raw, eps);
    else
        drift = 100*abs(results.brine_out.C - results.brine_out.C_win) ...
                / max(results.brine_out.C, eps);
    end
    results.salt_field_drift_pct = drift;
    info('Salt field drift', '%.2f', drift, '%', 'window mean vs final state, not screened');
end
foot();
end

function info(label, fmt, value, unit, note_str)
% Reported-but-not-screened quantity. Shares the verdict column with scr()
% so the screen block reads as one aligned list.
fprintf('%s\n', deblank(sprintf('  %-6s %-28s %12s  %-7s %s', ...
        '[INFO]', label, sprintf(fmt, value), unit, note_str)));
end

function scr(ok, label, fmt, value, unit, limit_str)
% One validity screen per line: verdict, measured value with unit, limit.
fprintf('%s\n', deblank(sprintf('  %s %-28s %12s  %-7s %s', ...
        ternary(ok,'[PASS]','[FAIL]'), label, sprintf(fmt, value), unit, limit_str)));
end



% ---- report formatting helpers --------------------------------------
function hdr(title_str, subtitle_str)
% Report block header. One rule width is used everywhere so the parts of
% the report line up when read end to end.
fprintf('\n');
fprintf('%s\n', '==============================================================================');
fprintf('  %s\n', title_str);
if ~isempty(subtitle_str), fprintf('  %s\n', subtitle_str); end
fprintf('%s\n', '==============================================================================');
end

function foot()
fprintf('%s\n', '==============================================================================');
end

function sec(name)
% Section heading inside a report block.
fprintf('\n  %s\n', name);
fprintf('  %s\n', '----------------------------------------------------------------------------');
end

function sub(name)
% Sub-heading inside a section. Sits one indent level left of the rows it
% introduces, so a heading is never mistaken for a reported quantity.
fprintf('\n  %s\n', name);
end

function prule()
fprintf('  %s\n', '----------------------------------------------------------------------------');
end

function prow(label, fmt, value, unit, note_str)
% One reported quantity per line: label, value, unit, optional annotation.
% The unit column is mandatory, so a printed number cannot appear without
% its dimension. Dimensionless quantities carry '-' and ratios carry '%'.
fprintf('%s\n', deblank(sprintf('    %-30s %12s  %-7s %s', ...
        label, sprintf(fmt, value), unit, note_str)));
end

function note(text)
% Annotation line, indented under the rows it qualifies.
fprintf('      %s\n', text);
end

function idn(identity_str, meaning_str)
% Identity heading in the stream table.
fprintf('    %-34s %s\n', identity_str, meaning_str);
end

function ideq(lhs, rhs_terms, residual)
% Prints an identity as evaluated numbers plus its residual, both in g/s.
txt_ideq = sprintf('%.4f =', lhs*1000);
for i_ideq = 1:numel(rhs_terms)
    if i_ideq == 1
        txt_ideq = sprintf('%s %.4f', txt_ideq, rhs_terms(i_ideq)*1000);
    else
        txt_ideq = sprintf('%s + %.4f', txt_ideq, rhs_terms(i_ideq)*1000);
    end
end
fprintf('      %s g/s\n', txt_ideq);
fprintf('      residual %+.2e g/s\n', residual*1000);
end
end
