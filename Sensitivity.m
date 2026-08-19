function T = Sensitivity(varargin)
%==========================================================================
%  ZLDD CASCADE -- LOCAL SENSITIVITY (OAT) AND FACTOR DECLARATION
%
%  Model  : ZLDD_complete_modeling.m
%  Author : Tuhin Mahmud
%
%  Everything the parametric study needs is in this one file: baseline, factor
%  list, runner, metric extraction, validity screens and the figure suite.
%
%  IT IS ALSO THE FACTOR DECLARATION FOR THE WHOLE PROJECT. Morris_Global reads
%  the ranges and setters from here through Sensitivity('factors') rather than
%  re-declaring them, so the OAT and Morris studies cannot describe different
%  factors while appearing to describe the same ones.
%
%  USAGE
%    T = Sensitivity('only',{'Qfan'});      % one factor only
%    T = Sensitivity();                     % full sweep
%    T = Sensitivity('block','uncertainty');% the SEC band block
%    T = Sensitivity('dryrun',true);        % list the runs, solve nothing
%    T = Sensitivity('plot',true);          % sweep, then draw the figure suite
%    T = Sensitivity('load');               % results table from OAT.mat, no re-solve
%    T = Sensitivity('load','my.mat');      % from a checkpoint under another name
%    Sensitivity('figs', T);                % figures from a table already in hand
%    Sensitivity('figs', T, 'F5');          % one panel only
%    Sensitivity('figs', T, {'F2a','F4'});  % a named subset
%    E = Sensitivity('elast', T);           % elasticity ranking table only
%    Sensitivity('tornado', T, 'R_still_pct');   % tornado on any response column
%    F = Sensitivity('factors');            % the factor declaration itself
%    Sensitivity('pinch', T);               % SEC against cascade recovery
%    Sensitivity('validity', T);            % validity envelope, on request only
%    Sensitivity('export','svg',figs);   % export a named set of figures
%    Sensitivity('export', 'eps');          % every open figure -> vector
%    T = Sensitivity('resume',false);       % ignore any checkpoint, re-solve
%    T = Sensitivity('ckpt','none');        % write nothing at all
%
%  OUTPUT AND CHECKPOINTING. The results table is returned in memory AND the
%  sweep is checkpointed to OAT.mat after every solve, holding four variables:
%
%    rows   N x 1 cell,    rows{r} = the metric struct for run r
%    done   N x 1 logical, done(r) = true once run r has been written
%    runs   1 x N struct,  the run list (id, factor, label, unit, block,
%                          value, is_base)
%    sig    the run-list signature the resume check matches against
%
%  rows is PREALLOCATED to the full N and filled in place, so a sweep stopped
%  part way leaves the completed entries populated, the remainder empty ([])
%  and done true only for those written. Re-running skips the completed
%  entries and continues at the first false. Table assembly uses rows(done),
%  so a partial file loads and plots.
%
%  ONE FILE, ONE NAME. Every selection writes OAT.mat. A single-factor call
%  and the full sweep therefore share it, and the later run replaces the
%  earlier -- copy OAT.mat aside by hand before switching if those results
%  still matter. Resuming requires an exact run-list match (same factors,
%  same levels, same order); anything else starts at run 1.
%  Pass 'ckpt','myfile.mat' to override the name, or 'none' to write nothing.
%
%  Baseline struct, factor list and run metadata are attached as
%  T.Properties.UserData. Export from the caller if a CSV is wanted:
%      T = Sensitivity();  writetable(T,'results.csv');
%
%  The figure and elasticity routines are LOCAL to this file, so they are
%  not callable as Sensitivity_figs(T) from the command line -- go through the
%  Sensitivity('figs',...) entry point above.
%
%  WHY Qfan FIRST. At baseline the air heater and compressor together are
%  72 % of purchased work and both scale with mdot_da, while cascade
%  evaporation falls only sub-linearly until the air approaches saturation.
%  SEC is therefore expected to FALL toward low flow and turn up only when
%  the loop air saturates mid-cascade. If the curve is still falling at
%  0.05 m^3/s the range does not contain the minimum and must be extended
%  DOWNWARD -- discovering that after the rest of the sweep wastes a day.
%
%  METRICS. Three are plotted:
%    SEC_external    [kWh/m3]  PLANT boundary: blower + crystallizer +
%                              compressor + air heater, over collected
%                              product. The only energy figure every factor
%                              can move -- the still-only SEC sits on a
%                              control volume that the bl.* and AWG factors
%                              cannot reach, so their bars would be
%                              structurally zero.
%    R_still_pct     [%]       cascade evaporation / feed water. R_plant is
%                              saturated at ~99.6 % against a 99.96 %
%                              ceiling, so it can only move down and by
%                              0.14 % -- table column, never a tornado bar.
%    prod_L_m2_day   [L/m2/d]  the only metric that PENALISES area, so it
%                              opposes SEC in the area sweep and makes the
%                              collector trade-off legible.
%
%  BASIS. Every metric is on the window-mean basis published by the model,
%  so numerator and denominator share one basis across all runs.
%==========================================================================

% ---- DISPATCH. Everything lives in this one file, so the figure and
% elasticity routines are LOCAL functions and are not callable from the
% command line directly. They are reached through this entry point instead:
%     Sensitivity('figs', T)        figures from an existing results table
%     Sensitivity('elast', T)       returns the elasticity table, no figures
% SINGLE-ARGUMENT VERBS ARE ALLOWED. The dispatch must not demand nargin >= 2:
% Sensitivity('factors') and Sensitivity('export') would then fall through to
% the option parser below, which reads varargin{i+1} for every odd i, so a
% lone verb gives an index error rather than a usable message. Option NAMES
% ('only', 'block', 'dryrun', 'plot') are not verbs, so they match no case
% and still fall through as intended.
if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1}))
    verb = lower(char(varargin{1}));
    needsT = {'figs','elast','tornado','pinch','validity'};
    if strcmp(verb,'extract') && numel(varargin) < 3
        error('Sensitivity:needsArgs', ...
              'Sensitivity(''extract'', results, FC) needs both arguments.');
    end
    if ismember(verb, needsT) && numel(varargin) < 2
        error('Sensitivity:needsTable', ...
              'Sensitivity(''%s'', T) requires the results table.', verb);
    end
    switch verb
        case 'figs'
            Tin = varargin{2};
            % Optional third argument names the panels to draw, e.g.
            % Sensitivity('figs',T,'F5') or Sensitivity('figs',T,{'F2a','F4'}).
            if numel(varargin) >= 3, sel = varargin{3}; else, sel = []; end
            Tin = refresh_labels(Tin);
            Sensitivity_figs(Tin, sel); T = Tin; return
        case 'elast'
            T = Sensitivity_elasticity(refresh_labels(varargin{2})); disp(T); return
        case 'tornado'
            T = varargin{2};
            if numel(varargin) >= 3, rsp = varargin{3}; else, rsp = 'SEC_external'; end
            T = refresh_labels(T);
            zldd_tornado(T, rsp); return
        case 'extract'
            % Expose the metric extractor so a second study records the SAME
            % ~50 quantities per solve as the sweep does. Without this the
            % Morris run would keep only the responses it happened to score,
            % and answering a new question later would mean re-solving.
            T = local_extract(varargin{2}, varargin{3}); return
        case 'factors'
            % Expose the ONE factor declaration so other studies (Morris)
            % consume it rather than re-declaring ranges that would then
            % drift out of step with the OAT sweep.
            T = local_factors(); return
        case 'load'
            % Rebuild the results table from a checkpoint without re-solving.
            % load('OAT.mat') alone puts rows/done/runs/sig in the workspace,
            % not T, and every figure entry point needs T. This performs the
            % same assembly the sweep performs at its end.
            if numel(varargin) >= 2, ck = char(varargin{2}); else, ck = 'OAT.mat'; end
            T = local_table_from_ckpt(ck); return
        case 'pinch'
            T = refresh_labels(varargin{2});
            zldd_pinch_plot(T); return
        case 'validity'
            T = refresh_labels(varargin{2});
            zldd_validity_envelope(T); return
        case 'export'
            if numel(varargin) >= 2, fmt = varargin{2}; else, fmt = 'eps'; end
            if numel(varargin) >= 3, figs = varargin{3}; else, figs = []; end
            export_figs(fmt, figs); T = table(); return
    end
end

AUTHOR = 'Tuhin Mahmud';   % stamped on the console header and the provenance

% ckpt   = '' uses the selection-derived default (see below). A explicit name
%          overrides it; '' after parsing cannot be requested, so pass 'none'
%          to disable checkpointing entirely.
% resume = read an existing matching checkpoint and skip the solves it holds.
opt = struct('only',{{}},'block','','dryrun',false,'plot',false, ...
             'ckpt','','resume',true);

% UNKNOWN OPTIONS ARE AN ERROR. Assigning straight into opt accepts any name,
% so a misspelling or an option from a newer version of this file is stored
% and then ignored, and the sweep runs with the default the caller thought
% they had overridden.
known = fieldnames(opt);
if mod(numel(varargin),2) ~= 0
    error('Sensitivity:pairs','Options must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    nm = varargin{i};
    if ~(ischar(nm) || isstring(nm)) || ~ismember(char(nm), known)
        error('Sensitivity:unknownOption', ...
              'Unknown option "%s". Known options: %s.', ...
              char(string(nm)), strjoin(known.', ', '));
    end
    opt.(char(nm)) = varargin{i+1};
end

F = local_factors();

% ---- BASE/STRUCT CONSISTENCY, CHECKED BEFORE THE FIRST SOLVE.
% Every declared factor base must equal the corresponding field of the
% baseline struct. A mismatch anchors every normalized curve and every
% elasticity to a run 1 that does not sit at the declared base, and does so
% silently. Only flat single-field factors can be checked this way; A_plate
% (sets L and W) and the bl.* factors are handled below.
FCchk = local_baseline();
for f = 1:numel(F)
    if isfield(FCchk, F(f).name)
        d = abs(FCchk.(F(f).name) - F(f).base);
        assert(d <= 1e-12*max(1,abs(F(f).base)), 'Sensitivity:baseMismatch', ...
            ['Factor %s: declared base %g disagrees with baseline struct ' ...
             '%g. Fix local_factors() or the baseline before sweeping.'], ...
            F(f).name, F(f).base, FCchk.(F(f).name));
    end
end
assert(abs(FCchk.L*FCchk.W - 1.0) <= 1e-12, 'Sensitivity:baseMismatch', ...
    'A_plate base is declared 1.0 m^2 but L*W = %g.', FCchk.L*FCchk.W);
for nm = {'f_heatloss','BPE'}
    i = find(strcmp({F.name}, nm{1}),1);
    if ~isempty(i) && isfield(FCchk,'bl') && isfield(FCchk.bl, nm{1})
        assert(abs(FCchk.bl.(nm{1}) - F(i).base) <= 1e-12*max(1,abs(F(i).base)), ...
            'Sensitivity:baseMismatch', 'Factor %s: base %g vs FC.bl %g.', ...
            nm{1}, F(i).base, FCchk.bl.(nm{1}));
    end
end
clear FCchk

keep = true(1,numel(F));
if ~isempty(opt.only),  keep = ismember({F.name}, opt.only);        end
if ~isempty(opt.block), keep = keep & strcmp({F.block}, opt.block); end
F = F(keep);
if isempty(F), error('Sensitivity:noFactors','No factors selected.'); end

% ---- run list. Strict OAT: exactly ONE factor moves from a common
% baseline, and the baseline appears ONCE as run 1, shared by every factor.
FC0  = local_baseline();
runs = struct('id',{},'factor',{},'label',{},'unit',{},'block',{}, ...
              'value',{},'is_base',{});
runs(1) = struct('id',1,'factor','BASELINE','label','Baseline','unit','-', ...
                 'block','baseline','value',NaN,'is_base',true);
n = 1;
for f = 1:numel(F)
    for L = F(f).levels(:)'
        if abs(L - F(f).base) <= 1e-12*max(1,abs(F(f).base)), continue; end
        n = n+1;
        runs(n) = struct('id',n,'factor',F(f).name,'label',F(f).label, ...
                         'unit',F(f).unit,'block',F(f).block, ...
                         'value',L,'is_base',false);
    end
end
N = numel(runs);

fprintf('\n=================================================================\n');
fprintf('  ZLDD PARAMETRIC / OAT SWEEP -- %d factors, %d runs\n', numel(F), N);
fprintf('  Author : %s\n', AUTHOR);
fprintf('  Run    : %s   |   Nx = %d (sweep grid)\n', ...
        datestr(now,'yyyy-mm-dd HH:MM:SS'), FC0.Nx);
fprintf('=================================================================\n');
for f = 1:numel(F)
    fprintf('  %-14s %-30s %2d levels  [%s]\n', F(f).name, F(f).label, ...
            numel(F(f).levels), F(f).block);
end
fprintf('-----------------------------------------------------------------\n');
if opt.dryrun, T = struct2table(runs); disp(T); return; end

rows = cell(N,1); done = false(N,1);

% ---- CHECKPOINT FILE. ONE FIXED NAME FOR EVERY SELECTION.
% ONE NAME MEANS ONE FILE. A single-factor call and the full sweep both write
% OAT.mat, so whichever ran last is what is on disk -- see the note under
% USAGE. The signature check below only decides whether the file on disk can
% be RESUMED; it never preserves it.
if isempty(opt.ckpt), opt.ckpt = 'OAT.mat'; end
if strcmpi(opt.ckpt,'none'), opt.ckpt = ''; end

% A CHECKPOINT MUST MATCH THE RUN LIST, NOT JUST ITS LENGTH. Two sweeps with
% the same number of runs but different levels are different studies, and
% splicing one onto the other gives a table that looks valid and is not.
% THE FIELD IS 'value'. The run list carries id/factor/label/unit/block/
% value/is_base; 'level' is the name that same number takes later on the
% metric struct (M.level = runs(r).value), not on runs itself.
sig = struct('factor',{{runs.factor}},'value',[runs.value], ...
             'is_base',[runs.is_base],'N',N);

if opt.resume && ~isempty(opt.ckpt) && exist(opt.ckpt,'file')
    S = load(opt.ckpt);
    ok_ckpt = isfield(S,'rows') && isfield(S,'done') && numel(S.done) == N;
    if ok_ckpt && isfield(S,'sig')
        ok_ckpt = isequaln(S.sig, sig);
    elseif ok_ckpt && isfield(S,'runs')
        % A FILE CARRYING NO SIGNATURE is matched on the run list itself, so
        % solves already on disk are not thrown away.
        ok_ckpt = isequal({S.runs.factor},{runs.factor}) && ...
                  isequaln([S.runs.value],[runs.value]);
    else
        ok_ckpt = false;
    end
    if ok_ckpt
        rows = S.rows(:); done = logical(S.done(:));
        fprintf('  Resuming from %s -- %d of %d runs already complete.\n', ...
                opt.ckpt, sum(done), N);
    else
        % A DIFFERENT RUN LIST IS A DIFFERENT STUDY. It cannot be resumed onto
        % this one, so the solves start at 1 and the first write replaces the
        % file. Copy OAT.mat aside by hand before switching selections if
        % those results still matter.
        fprintf('  %s holds a different run list -- starting fresh (it will be replaced).\n', ...
                opt.ckpt);
    end
end
if ~isempty(opt.ckpt)
    fprintf('  Checkpoint: %s (written after every run)\n', opt.ckpt);
end

vis0 = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');
cleanupObj = onCleanup(@() set(0,'DefaultFigureVisible',vis0)); 

t_start = tic;
for r = 1:N
    if done(r), continue; end
    FC = FC0;
    FC.verbose = false;                 % batch mode: suppress the per-run audit
    if ~runs(r).is_base
        idx = find(strcmp({F.name}, runs(r).factor),1);
        FC  = F(idx).setter(FC, runs(r).value);
    end

    fprintf('[%3d/%3d] %-14s = %-9.4g ', r, N, runs(r).factor, runs(r).value);
    t0 = tic; err = '';
    try
        evalc('res = ZLDD_complete_modeling(FC);');   % swallow residual prints
        M = local_extract(res, FC); clear res
    catch ME
        M = local_extract([], FC); err = ME.message;
    end
    close all force

    M.run_id  = runs(r).id;        M.factor  = string(runs(r).factor);
    M.fac_lbl = string(runs(r).label);
    M.block   = string(runs(r).block);
    M.level   = runs(r).value;     M.is_base = runs(r).is_base;
    M.wall_s  = toc(t0);           M.error   = string(err);

    rows{r} = M; done(r) = true;

    % WRITTEN INSIDE THE LOOP, NOT AFTER IT. That is the entire point: a
    % sweep killed at run 60 resumes at run 60. -v7 matches the format of the
    % existing zldd_OAT_all.mat and costs ~90 kB and a few ms per run.
    if ~isempty(opt.ckpt)
        save(opt.ckpt,'rows','runs','done','sig','-v7');
    end

    if isempty(err)
        fprintf('SEC %7.2f  R_still %6.2f%%  prod %6.2f kg/d  %s (%.0fs)\n', ...
            M.SEC_external, M.R_still_pct, M.prod_kg_day, local_flags(M), M.wall_s);
    else
        fprintf('FAILED: %s\n', err);
    end
end
set(0,'DefaultFigureVisible',vis0);

% TABLE ASSEMBLY. local_extract() returns one fixed schema on every path,
% including the failure path, so a plain vertcat is safe. orderfields()
% against the prototype guards the field ORDER as well, since vertcat is
% order-sensitive on older releases.
proto = local_blank();
kept  = rows(done);
for i = 1:numel(kept), kept{i} = orderfields(kept{i}, proto); end
T = struct2table(vertcat(kept{:}));
T = movevars(T,{'run_id','factor','fac_lbl','block','level','is_base'},'Before',1);

% PROVENANCE. The table alone does not say which baseline, which grid or which
% MATLAB produced it, so it is carried with the table rather than alongside it.
T.Properties.UserData = struct( ...
    'author', AUTHOR, 'baseline_FC', FC0, 'factors', {F}, 'runs', {runs}, ...
    'matlab_version', version, 'script', mfilename('fullpath'), 'run_utc', ...
    datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS'));

% THE ASSEMBLED TABLE IS WRITTEN INTO THE CHECKPOINT. rows and done are the
% resume state and are written after every solve; T is the product and is
% written once, here, when the sweep completes. A later session then gets the
% table straight from the file without repeating the assembly, and T carries
% the provenance struct that rows alone does not hold.
%
% Appended rather than re-saved, so the resume state written during the sweep
% is left untouched. A partial run writes no T, which is what distinguishes a
% completed checkpoint from an interrupted one.
if ~isempty(opt.ckpt) && all(done)
    save(opt.ckpt,'T','-append');
    fprintf('  Results table T written to %s.\n', opt.ckpt);
end

fprintf('-----------------------------------------------------------------\n');
fprintf('  %d runs in %.1f min.\n', height(T), toc(t_start)/60);
nf = sum(~T.ok_all);
fprintf('  %d run(s) tripped a validity screen.\n', nf);
if nf > 0
    fprintf('  Exclude these from the response curves and report them as the\n');
    fprintf('  validity envelope:\n');
    bad = T(~T.ok_all,:);
    for i = 1:height(bad)
        fprintf('    %-14s = %-9.4g %s\n', bad.factor(i), bad.level(i), ...
                local_flags(table2struct(bad(i,:))));
    end
end
if sum(~T.ok_reynolds) > 0
    fprintf('  %d run(s) below turbulent gap-Re: still plotted, but the text\n', ...
            sum(~T.ok_reynolds));
    fprintf('  must say the closure is used outside its fitted range.\n');
end
if sum(~T.ok_salinity) > 0
    fprintf('  %d run(s) above the property-correlation salinity limit.\n', ...
            sum(~T.ok_salinity));
end
if all(T.cons_basis == "absent")
    fprintf('  Conservation screen inactive: no invariant field was found in\n');
    fprintf('  the model output on any run.\n');
elseif any(T.cons_basis ~= "absent")
    fprintf('  Max conservation residual over the sweep: %.3e (%d run(s) fail).\n', ...
            max(T.cons_resid_max), sum(~T.ok_conserve));
end
if any(T.stack_basis == "fallback_lower_bound")
    fprintf('  Stack height is a lower bound on %d run(s): the gap-volume\n', ...
            sum(T.stack_basis == "fallback_lower_bound"));
    fprintf('  field was unavailable, so the top compartment is not counted.\n');
end
if sum(~T.ok_stack) > 0
    fprintf('  %d run(s) exceed the buildable stack height (%.2f m max seen).\n', ...
            sum(~T.ok_stack), max(T.stack_h_m));
end
if any(T.Nx ~= 40)
    fprintf('  Grid: sweep run at Nx = %g, not the Nx = 40 of the baseline.\n', ...
            mode(T.Nx));
    fprintf('  Read results RELATIVE to run 1; absolute SEC values from this\n');
    fprintf('  table are not comparable with the baseline figures.\n');
end
fprintf('=================================================================\n\n');

if opt.plot
    Sensitivity_figs(T);
    disp(Sensitivity_elasticity(T));
end
end


%%=========================================================================
%%  BASELINE
%%=========================================================================
function FC = local_baseline()
% ONE DECLARATION, AND ONLY ONE. The baseline is read from
% Baseline_Input_Variable.m, which main.m and the TEA also read. No duplicate
% hard-coded copy is kept here as a fallback: a fallback is precisely how two
% baselines drift apart, because
% the sweep runs happily against the stale copy and nothing is raised. If the
% declaration is missing, that is an error, not something to work around.
%
% B3. THE GRID IS IMPOSED HERE, NOT INHERITED. The shared declaration carries
% the PUBLISHED grid, Nx = 40. Returning it unmodified would leave the
% deliberate Nx = 20 sweep economy documented in local_factors() silently
% unapplied, quadrupling the cost of the campaign with no message. The shared file
% remains the single source of truth for the PHYSICS; the grid is a SWEEP
% decision and is imposed here, explicitly, in one place. M.Nx records what
% was actually used on every row, so the provenance survives into the table.
NX_SWEEP = 20;                 % see the Nx note in local_factors()

if exist('Baseline_Input_Variable','file') ~= 2
    error('Sensitivity:noBaseline', ...
         ['Baseline_Input_Variable.m is not on the path. The sweep will not ' ...
          'fall back to a private copy of the baseline -- a duplicate ' ...
          'declaration drifts out of step silently and invalidates every ' ...
          'row without raising anything. Add it to the path and re-run.']);
end

FC = Baseline_Input_Variable();
FC.Nx      = NX_SWEEP;   % sweep economy, NOT the published grid
FC.verbose = false;      % 80 verbose audit reports are unreadable
end


%%=========================================================================
%%  FACTORS
%%=========================================================================
function F = local_factors()
% LEVELS RESTRICTED WHERE THE BALANCE DOES NOT CLOSE. Some combinations
% report plant recovery above 100 %, which is not a design being infeasible
% but a water balance that does not close: T_coil >= 298, T_air_in = 330,
% Vfeed = 80, Tfeed = 290 and Ta <= 305. Those levels are excluded here so
% the study reports only points the model can account for.
%
% THE RESTRICTION HIDES THE SYMPTOM, IT DOES NOT FIX THE CAUSE. The baseline
% itself sits at 99.43 %, half a percentage point below the same limit, and
% the conservation screen that would test closure directly does not yet run --
% cons_basis reports "absent" on every case, so ok_conserve passes by
% construction. Until the field names in conservation_residual() are wired
% and the residual is checked, treat closure as unverified rather than good.
%
% LEVEL DENSITY FOLLOWS MEASURED CURVATURE, not a uniform count. The measure
% is each factor's departure from a straight line as a fraction of its own
% SEC range: largest for A_plate, Qfan and T_air_in, mid-range for Vfeed,
% theta, hcomp and Np, small for T_coil, Tfeed_target and TDSfeed, and
% negligible for eps_recup. Points are spent where the curve bends and saved
% where it does not -- a dense grid on a response that is straight to
% 0.07 kWh/m3 buys nothing.
%
% Qfan and A_plate carry the densest grids because both have an interior
% optimum, and Qfan's spacing is DELIBERATELY NON-UNIFORM: 0.01 steps through
% 0.10-0.20 around the minimum near 0.145, coarser outside it. Uniform
% spacing over 0.03-0.30 would straddle the minimum and round it off.
%
% Dew-point failures were NOT removed on the same grounds: those margins are
% 1.9-2.9 K, positive but under a ~3 K design guard, which is the guard doing
% its job rather than the model failing.
% Setters are function handles, not field-name strings, because two factors
% are not single fields: plate area sets L and W together as sqrt(A), and
% bl.* cannot be reached by flat assignment.
%
% BLOCKS
%   design       levers you specify. Their tornado ranks DESIGN LEVERAGE,
%                not epistemic uncertainty -- say so in the section opening.
%   scenario     feedwater you are handed rather than choose.
%   uncertainty  2-level, reported as a band on SEC, NOT in the design
%                tornado. With kappa_w/Ta/Vwind out of the design block this
%                is the ONLY uncertainty quantification behind the headline.
%
% EXCLUDED, and why (this text belongs in the paper):
%   beta    must equal the PVGIS slope; perturbing it alone decouples the
%           glazing from the irradiance series.
%   Nx      grid parameter, not physics. HELD AT 20 FOR THE SWEEP, not at the
%           40 used for the published baseline, to keep the sweep affordable.
%           This is a deliberate trade and must be declared, because the
%           GCI study was performed at the finer grid and does NOT cover
%           these runs.
%           WHAT IT COSTS. Advection uses first-order upwinding, so the
%           leading discretisation error goes as dx, i.e. as 1/Nx: halving
%           Nx from 40 to 20 roughly DOUBLES it. Scaling the reported
%           GCI_80 = 0.24 % by that argument puts the Nx = 20 error near
%           1 % on absolute quantities. CONFIRM this against your own grid
%           study rather than taking the scaling on trust.
%           WHY IT IS STILL DEFENSIBLE. The sweep's conclusions are
%           DIFFERENCES between runs sharing one grid, and a systematic
%           discretisation bias largely cancels in a ratio Y/Y0. So the
%           RANKING and the SHAPE of the response curves survive; the
%           ABSOLUTE SEC of any single run does not, and must not be quoted
%           against the Nx = 40 figures in section 5.
%           REPORT RELATIVE, NOT ABSOLUTE. Normalize every sweep result to
%           the Nx = 20 baseline (run 1), which F2 already does.
%           VERIFY THE RANKING IS GRID-INDEPENDENT before publishing it:
%           re-run the two endpoint levels of the three highest-ranked
%           factors at Nx = 40 -- a handful of solves -- and check the ordering
%           does not change. If it does, the sweep must be repeated at 40.
%   eta_heater / eta_air_heater   resistance elements, ~1.0 by construction.
%   h_cryst thermochemical constant (~66 kJ/kg NaCl), not an estimate.
%   eta_II  already propagated as a band inside the model; sweeping it too
%           would double-count the compressor term.
%
% LEVEL SPACING is deliberately non-uniform:
%   A     geometric (/4 to x4) -- linear levels would cluster in the top
%         half and under-resolve the small-area end where curvature is.
%   Qfan  denser below baseline -- range is x10 and blower work goes roughly
%         as Q^3, so the turnover sits in the bottom third.
k = 0;

k=k+1; F(k) = mk('A_plate','Plate area','m^2','design',1.0, ...
    [0.16 0.2025 0.25 0.3025 0.36 0.4225 0.49 0.5625 0.64 0.7225 1.0 ...
     0.81 0.9025 1.1025 1.21 1.3225 1.44 1.5625 1.69 1.8225 1.96], @(FC,v) setLW(FC,v), ...
    ['L = W = sqrt(A), aspect ratio 1. Vfeed HELD CONSTANT, so this reads ' ...
     'as specific-area intensification -- scaling feed with area would be ' ...
     'a two-factor block move, not OAT. h_mean = hcomp + tan(theta)*L ' ...
     'rides along, so the bar carries a gap-diffusion term as well as area.']);

k=k+1; F(k) = mk('Vfeed','Feed flow','L/day','design',100, ...
    [90 95 100 110 120 130 150], @(FC,v) sf(FC,'Vfeed',v), ...
    ['Sets both the recovery denominator and the film loading. Spans the ' ...
     'same specific-loading group as A_plate, so the two bars will be near ' ...
     'mirror images -- report as ONE finding via the Gamma* collapse plot.']);

k=k+1; F(k) = mk('Qfan','Fan flow','m^3/s','design',0.25, ...
    [0.03 0.045 0.06 0.075 0.09 ...
     0.10 0.11 0.12 0.13 0.14 0.145 0.15 0.16 0.17 0.18 0.19 0.20 ...
     0.215 0.23 0.25 0.275 0.30], ...
    @(FC,v) sf(FC,'Qfan',v), ...
    ['RUN FIRST. Expect a minimum near 0.10-0.20. At the low end watch for ' ...
     'the loop air saturating mid-cascade: at baseline gap 1 already sits ' ...
     'at w = 0.0238 against w_sat ~0.0295, so a fifth of the flow may ' ...
     'saturate and the deep stages stop contributing. ' ...
     'LOW-END FLOOR IS SET BY THE CLOSURE, NOT BY PREFERENCE. Re scales ' ...
     'linearly with Qfan and baseline sits near Re = 2.8e4, so 0.03 m^3/s ' ...
     'holds every gap near Re = 3.4e3 -- above the 2300 blend threshold of ' ...
     'the Nu smoothstep. Going lower crosses onto the laminar branch ' ...
     'mid-sweep and the branch change, not the physics, would produce the ' ...
     'curvature. ok_reynolds is captured per run; if it trips at 0.03 the ' ...
     'floor must be raised, not the screen relaxed. ' ...
     'SPACING IS GEOMETRIC (ratio ~1.27 below baseline) because blower work ' ...
     'goes as Q^3 while heater and compressor go as mdot_da, so the minimum ' ...
     'is SHALLOW: three points across 0.10-0.20 cannot locate it. Six ' ...
     'levels now sit below 0.20 where the turnover is expected. ' ...
     'UPPER LIMIT 0.35, NOT 0.50. Above 0.35 the extra evaporation drives ' ...
     'the terminal brine to salt saturation and the run is screened out by ' ...
     'ok_brine, so levels beyond it buy failed runs rather than curve. ' ...
     '0.35 IS THE BOUNDARY ITSELF, so it is run as a probe: if it trips ' ...
     'ok_brine or ok_salinity, report it as the located saturation limit ' ...
     'rather than deleting it -- a screened run at a known constraint is a ' ...
     'result, and it is the only level that bounds the feasible range from ' ...
     'above. (The level was ARGUED FOR HERE BUT ABSENT FROM THE LIST, which ' ...
     'topped out at 0.30; it has been added.) The bar is ASYMMETRIC about ' ...
     'baseline (0.03-0.25 down, 0.25-0.35 up) and its raw length is NOT ' ...
     'comparable with symmetric factors -- rank on S_local, not bar length. ' ...
     'AND NOT ON S_fit EITHER: this factor is expected to be non-monotonic, ' ...
     'so the range-averaged exponent is meaningless for it by construction.']);

k=k+1; F(k) = mk('theta','Film tilt','deg','design',2, ...
    [1 2 3 4 5], @(FC,v) sf(FC,'theta',v), ...
    'The real gap-height lever; hcomp is only an offset.');

k=k+1; F(k) = mk('hcomp','Gap height offset','m','design',0.02, ...
    [0.005 0.010 0.015 0.020 0.030 0.040 0.050], @(FC,v) sf(FC,'hcomp',v), ...
    ['THE OTHER GAP-HEIGHT LEVER, and the more direct one. Mean gap height ' ...
     'is h_bar = hcomp + tan(theta)*L_ch, so at baseline (theta = 2 deg, ' ...
     'L_ch ~ 1.1 m) the taper contributes ~0.038 m against the 0.020 m ' ...
     'offset -- hcomp is about a third of h_bar and moving it 0.005-0.050 m ' ...
     'swings h_bar by roughly a factor of two. ' ...
     'MECHANISM. For a wide duct (h << W) the hydraulic diameter is ' ...
     'D_h ~ 2h, so h_c = Nu*k_a/D_h scales as 1/h: a TIGHTER gap raises the ' ...
     'convective and, through Chilton-Colburn, the mass-transfer ' ...
     'coefficient. Note that Re is nearly INVARIANT to h at fixed mass ' ...
     'flow, since Re = 2*mdot/(mu*(W+h)) -- so ok_reynolds will not bind ' ...
     'across this sweep and the response is pure transfer-coefficient, not ' ...
     'a regime change. ' ...
     'READ THE LOW END WITH SUSPICION. If blower work is computed from ' ...
     'Qfan alone rather than from a gap pressure drop, SEC will fall ' ...
     'MONOTONICALLY towards hcomp = 0.005 with no penalty, because the ' ...
     'model gets the higher h_c for free. Real duct loss goes roughly as ' ...
     'h^-3 at fixed volumetric flow, so a monotonic curve here is a MODEL ' ...
     'ARTEFACT, not an optimum, and must be reported as such. CHECK how ' ...
     'W_blower is formed before quoting any hcomp recommendation. ' ...
     'PARTIALLY REDUNDANT WITH theta: both set h_bar, so the two bars are ' ...
     'not independent and should be discussed together, not ranked against ' ...
     'each other. ' ...
     'UPPER LEVEL is bounded by buildable stack height, ' ...
     'roughly Np*(hcomp + t_plate); at 0.050 m and Np = 10 the stack is ' ...
     'already ~0.5 m before the tapered volume is counted.']);

k=k+1; F(k) = mk('T_coil','AWG coil temperature','K','design',295, ...
    [282 285 288 291 295], @(FC,v) sf(FC,'T_coil',v), ...
    ['Swept over its admissible range, NOT by percent -- absolute ' ...
     'temperature is not ratio-scaled (+/-50 % of 295 K is meaningless). ' ...
     'Acts twice: sets w_in = w_sat(T_coil) and the recuperator cold inlet, ' ...
     'and sets the heat-pump lift. Expect it to split SEC and R_still in ' ...
     'OPPOSITE directions -- the most interesting result in the block. ' ...
     '308 K may trip the dew-point guard; that is a result, not a bad run.']);

k=k+1; F(k) = mk('Np','Number of plates','-','design',10, ...
    [5 7 10 12 15], @(FC,v) sf(FC,'Np',v), ...
    ['Fills the [To be completed] placeholder in section 2 of the baseline ' ...
     'document. Np = 20 at low Vfeed is the likeliest breach of the ' ...
     'recovery ceiling.']);

k=k+1; F(k) = mk('T_air_in','Still inlet air temperature','K','design',325, ...
    [315 316 317 318 319 320 321 322 323 324 325], @(FC,v) sf(FC,'T_air_in',v), ...
    ['Absolute-range sweep. Drives the single largest load (60 % of ' ...
     'purchased work), so expect it at the top of the SEC tornado.']);

% k=k+1; F(k) = mk('Tfeed_target','Feed preheat target','K','design',350, ...
%     [330 335 340 345 350 355 360 365], @(FC,v) sf(FC,'Tfeed_target',v), ...
%     ['Capped near 368 K by dT_pinch_HX against T_sat_vap = 373.15, so 365 ' ...
%      'is the top safe level: above it the regime switches to pinch-limited ' ...
%      'and the bar flattens for the wrong reason. regime_feedHX is captured ' ...
%      'per run so this is visible rather than assumed.']);


k=k+1; F(k) = mk('Tfeed_target','Feed preheat target','K','design',320, ...
    [305 315 320 330 340 350], @(FC,v) sf(FC,'Tfeed_target',v), ...
    ['RANGE STRADDLES THE BASELINE: 305-350 K about a 320 K base. The ' ...
     'earlier 301-320 window sat entirely at or below baseline, so the ' ...
     'response curve never passed through run 1 and the bar was one-sided ' ...
     '-- not comparable with the two-sided bars of Qfan, T_coil and ' ...
     'T_air_in on the same chart. ' ...
     'FC.Tfeed_target IN Baseline_Input_Variable.m IS 320 AND MUST STAY IN ' ...
     'STEP with the base argument above. The base-consistency assertion at ' ...
     'the head of Sensitivity now refuses to run on a mismatch, so this can no ' ...
     'longer anchor every bar to the wrong run 1 in silence. ' ...
     'UPPER LIMIT is dT_pinch_HX below T_sat_vap = 373.15, i.e. 368 K, so ' ...
     '365 is the top safe level; above it the regime switches to ' ...
     'pinch-limited and the bar flattens for the wrong reason. ' ...
     'THE INERT UPPER BRANCH IS THE RESULT, NOT A FAILED RUN. If the feed ' ...
     'exchanger is capacity-limited the achieved Tfeed_casc saturates and ' ...
     'targets above it do nothing; regime_feedHX and the Tfeed_casc column ' ...
     'make that visible per run. Plot Tfeed_casc against Tfeed_target and ' ...
     'the departure from the 1:1 line IS the finding -- report the ' ...
     'binding-to-inert transition temperature explicitly.']);

k=k+1; F(k) = mk('eps_recup','Recuperator effectiveness','-','design',0.75, ...
    [0.70 0.75 0.85], @(FC,v) sf(FC,'eps_recup',v), ...
    ['eps = 0 disables the unit and quantifies its value directly -- likely ' ...
     'an abstract-worthy number. Turns the "both duties shrink together" ' ...
     'claim in the model comments into a figure.']);

k=k+1; F(k) = mk('Tfeed','Feed supply temperature','K','design',300, ...
    [300 305 310], @(FC,v) sf(FC,'Tfeed',v), ...
    ['Expected nearly inert while the HX is target-limited: the preheat ' ...
     'loop absorbs the change. A SHORT bar here is a positive result and ' ...
     'deserves one sentence -- it demonstrates the loop works.']);

k=k+1; F(k) = mk('TDSfeed','Feed salinity','kg/m^3','scenario',35, ...
    [5 20 35 55 75], @(FC,v) sf(FC,'TDSfeed',v), ...
    ['Denser through 35-45 where a real Gulf installation sits; coarser at ' ...
     'the ends where the point is range, not resolution. NOT a design lever ' ...
     '-- own subsection and figure, never the same tornado as Qfan. Two ' ...
     'cautions: recovery is defined against a denominator that SHRINKS with ' ...
     'salinity, so plot absolute distillate beside it; and confirm the ' ...
     'psat_saline fit limit at 75 (it is usually narrower than the density ' ...
     'fit, and is the binding one).']);

k=k+1; F(k) = mk('kappa_w','Film absorption coefficient','1/m','uncertainty',300, ...
    [150 225 300 425 600], @(FC,v) sf(FC,'kappa_w',v), ...
    ['Wide (/2 to x2): a lumped gray fit to a strongly spectral medium, ' ...
     'with a directional bias the model already concedes. Controls the ' ...
     'split of absorbed solar between film and plate. Expect a modest ' ...
     'effect on recovery and almost none on SEC -- if that holds, the known ' ...
     'bias does not threaten the headline number, which is the finding.']);

k=k+1; F(k) = mk('f_heatloss','Crystallizer vessel loss','-','uncertainty',0.15, ...
    [0.05 0.15 0.25], @(FC,v) ss(FC,'bl','f_heatloss',v), ...
    ['A pure engineering estimate with no stated basis, multiplying the ' ...
     'second-largest purchased load. Expect the largest SEC response here. ' ...
     'A parameter you invented needs a stated band.']);

k=k+1; F(k) = mk('BPE','Boiling-point elevation','K','uncertainty',8.5, ...
    [7.0 7.75 8.5 9.25 10.0], @(FC,v) ss(FC,'bl','BPE',v), ...
    ['Expect almost no propagation: the vapor condenses isothermally at ' ...
     'T_sat_vap so BPE never reaches the feed HX, and only lifts ' ...
     'T_vap_heater against a latent-dominated duty. DEMONSTRATING that is ' ...
     'worth more than the number -- it substantiates the constant-BPE ' ...
     'assumption a crystallizer reviewer will question.']);

k=k+1; F(k) = mk('Ta','Ambient temperature','K','uncertainty',308, ...
    [308 309 310 311.5 313], @(FC,v) sf(FC,'Ta',v), ...
    ['CONFIRMATION ONLY. Ta was insensitive on the once-through flowsheet, ' ...
     'but the closed loop changed the loss structure and Ta still sets the ' ...
     'AWG reject temperature, hence the lift. Two runs to confirm, then one ' ...
     'sentence excluding it -- cheaper than an assumption. ' ...
     'CONFIRMED BY THE BASELINE RUN: reject 313.00 K = Ta + 5.0 K approach, ' ...
     'lift 18.00 K, COP 7.38, compressor 774.6 W. Ta reaches SEC through the ' ...
     'compressor term and through the envelope loss, and through nothing ' ...
     'else -- no ambient air enters the loop, so RH_amb and w_amb do not ' ...
     'ride along with it. ' ...
     'ONE COUPLING IS BROKEN BY THE PERTURBATION: Tground is held at 305 K ' ...
     'while Ta moves +/-5 K, although soil and ambient are physically ' ...
     'linked. The floor conduction term is small (0.58 W at baseline), so ' ...
     'the distortion is second-order, but state it rather than let a ' ...
     'reviewer find it.']);
end

function s = mk(n,l,u,b,base,lev,set,note)
s = struct('name',n,'label',l,'unit',u,'block',b,'base',base, ...
           'levels',lev,'setter',set,'note',note);
end
function FC = sf(FC,f,v),        FC.(f) = v;            end
function FC = ss(FC,g,f,v),      FC.(g).(f) = v;        end
function FC = setLW(FC,A),       FC.L = sqrt(A); FC.W = sqrt(A); end


%%=========================================================================
%%  EXTRACTION AND VALIDITY SCREENS
%%=========================================================================
function M = local_extract(results, FC)
% One flat struct per run: metrics plus a flag for every screen the model
% prints. At x4 and /4 ranges a real fraction of runs WILL leave the
% validated envelope; capturing that at run time is the difference between
% a reported validity envelope and unexplained gaps in a curve.
%
% B1. ONE SCHEMA ON EVERY PATH. A failure branch returning a struct of two
% fields against roughly fifty on success is fatal: struct arrays cannot be
% concatenated across differing field sets, so a SINGLE failed run destroys
% the whole table -- at the end, after every solve has already been paid
% for. Since the sweep deliberately probes outside the
% validated envelope (x4 and /4 ranges, a saturation probe at Qfan = 0.35),
% failed runs are expected, not exceptional. The blank prototype below is
% returned unmodified on failure: NaN metrics, false screens, and a row that
% still carries its factor and level so the gap is visible in F6 rather than
% absent from the table.
M = local_blank();
if isempty(results) || ~isstruct(results), return; end

g   = @(s,f,d) gd(s,f,d);
sys = g(results,'system',struct());  pr = g(results,'product',struct());
awg = g(results,'awg',struct());     bo = g(results,'brine_out',struct());
fp  = g(results,'feedpre',struct()); D  = g(results,'downstream',struct());

% ---- primary ----
M.SEC_external  = g(sys,'SEC_external',NaN);
M.R_still_pct   = g(results,'water_recovery_pct',NaN);
M.prod_kg_day   = g(pr,'kg_per_day',NaN);
M.evap_kg_day   = g(results,'mfw_daily',NaN);
M.Nx            = g(FC,'Nx',NaN);   % provenance: the grid every metric sits on
A               = g(FC,'L',NaN)*g(FC,'W',NaN);
M.A_plate_m2    = A;
M.A_total_m2    = A*g(FC,'Np',NaN);
M.prod_L_m2_day = M.prod_kg_day / max(M.A_total_m2,eps);
M.Gamma_star    = g(FC,'Vfeed',NaN) / max(M.A_total_m2,eps);  % collapse group

% Mean gap height, the group hcomp and theta BOTH move. Captured so the two
% bars can be shown to be different routes to one variable rather than two
% independent levers -- plot SEC against h_bar and the hcomp and theta points
% should fall on a common curve. If they do NOT, the difference is the taper
% gradient itself and is worth a sentence.
% Per section 2.6.1 the gap spans hcomp + 2*tan(theta)*c to
% hcomp + 2*tan(theta)*(L_ch - c), so the mean is hcomp + tan(theta)*L_ch.
L_ch            = g(FC,'L',NaN) / max(cosd(g(FC,'beta',0)),eps);
M.h_bar_m       = g(FC,'hcomp',NaN) + tand(g(FC,'theta',NaN))*L_ch;
[M.stack_h_m, M.stack_basis] = stack_height(results, FC, M.h_bar_m);

% ---- secondary ----
M.R_plant_pct   = g(pr,'recovery_pct',NaN);
M.SEC_blower    = g(sys,'SEC_blower',NaN);
M.SEC_heater    = g(sys,'SEC_heater',NaN);
M.SEC_compr     = g(sys,'SEC_compressor',NaN);
M.SEC_airheat   = g(sys,'SEC_air_heater',NaN);
M.W_external    = g(sys,'W_external',NaN);
band            = g(sys,'SEC_external_band',[NaN NaN]);
M.SEC_band_lo   = min(band); M.SEC_band_hi = max(band);
M.brine_TDS     = g(bo,'C',NaN);      M.CF_cascade = g(bo,'CF',NaN);
M.Tfeed_casc    = g(fp,'Tfeed_out',NaN);
M.f_cond        = gn(results,{'wmean','f_cond'}, g(fp,'f_cond',NaN));
M.regime_feedHX = string(g(fp,'regime',''));
M.Tv_exit       = g(awg,'Tv_exit',NaN);
M.COP_awg       = g(awg,'COP',NaN);   M.lift_awg = g(awg,'lift',NaN);
M.residence_min = g(results,'retention_time_total',NaN)/60;

% ---- screens ----
te = g(results,'t',NaN); if ~isscalar(te), te = te(end); end
M.t_end   = te;
M.ok_run  = isfinite(te) && te >= 0.999*g(FC,'t',NaN);
M.ok_solver   = all(isfinite([M.SEC_external M.R_still_pct M.prod_kg_day]));
M.ok_massbal  = isfinite(M.R_still_pct) && M.R_still_pct <= 100;

% Conservation residual, per run. Tolerance 1e-8 on the normalized invariant
% is loose against the machine-precision closure reported at the design point
% and tight enough to catch a corner where the solve has degraded. When the
% model exposes no invariant vector the screen passes by construction and
% cons_basis records 'absent' -- do NOT claim invariant closure across the
% sweep in the paper unless cons_basis is something other than 'absent'.
[M.cons_resid_max, M.cons_basis] = conservation_residual(results);
M.ok_conserve = (M.cons_basis == "absent") || ...
                (isfinite(M.cons_resid_max) && M.cons_resid_max <= 1e-8);

mfw   = g(FC,'Vfeed',NaN)*g(results,'rhow_in',1000)/1000 ...
      - g(FC,'Vfeed',NaN)*g(FC,'TDSfeed',NaN)/1000;
wt    = gn(FC,{'bl','w_salt_target'},0.99);
msalt = g(FC,'Vfeed',NaN)*g(FC,'TDSfeed',NaN)/1000;
M.recovery_ceil = 100*(mfw - (msalt/max(wt,eps))*(1-wt))/max(mfw,eps);
M.ok_recovery   = isfinite(M.R_plant_pct) && M.R_plant_pct <= M.recovery_ceil+1e-6;

if isfield(results,'C_out') && ~isempty(results.C_out)
    M.TDS_traj_max = max(results.C_out(:,end));
else
    M.TDS_traj_max = NaN;
end
M.ok_brine = isfinite(M.TDS_traj_max) && M.TDS_traj_max <= 317;

% Dew-point margin at the COIL INLET. If the recuperator already drained
% condensate the air arrives saturated at its hot outlet, so screening the
% raw exhaust would report a margin the coil never sees.
if g(D,'mdot_cond_recup',0) > 0 && isfield(D,'T_recup_hot_out')
    Tdci = D.T_recup_hot_out;
else
    Tdci = g(awg,'T_dew_return',NaN);
end
M.dew_margin_K = Tdci - g(awg,'T_coil',NaN);
M.ok_dewpoint  = isfinite(M.dew_margin_K) && M.dew_margin_K > 3;

% Startup preheat cap must be RELEASED by t_end. If it still clamps, the
% feed temperature is (Tfeed + cap) -- a tuning constant, not a solved value.
M.dT_preheat    = M.Tfeed_casc - g(FC,'Tfeed',NaN);
M.ok_preheatcap = ~(isfinite(M.dT_preheat) && ...
                    M.dT_preheat >= 0.999*g(FC,'dT_preheat_cap',Inf));

M.ok_secband = ~isfinite(M.SEC_external) || ...
    (M.SEC_external >= M.SEC_band_lo-1e-9 && M.SEC_external <= M.SEC_band_hi+1e-9);

if isfield(results,'Re_gap') && ~isempty(results.Re_gap)
    M.Re_gap_min = min(results.Re_gap(end,:));
else
    M.Re_gap_min = NaN;
end
M.ok_reynolds = isfinite(M.Re_gap_min) && M.Re_gap_min >= 2300;

% UNIT BASIS -- THE TWO SCALES MUST NOT BE MIXED. TDS_traj_max is kg/m3 (the C
% carried by Ms = C*delta). The correlation limit is 150 g/kg SOLUTION (Sg),
% the binding member being viscosity per section 2.11. The two scales are NOT
% interchangeable -- exactly the Sg/Sm confusion the manuscript warns about.
% Converting at the local solution density, 0.150 * 1080 = 162 kg/m3. A
% threshold of 120 kg/m3 is over-strict by ~35 % and flags the BASELINE itself
% as invalid, since the reported terminal brine is ~125 kg/m3.
% CONFIRM before submission which correlation is binding at the top of the
% TDSfeed sweep -- psat_saline is often narrower than viscosity.
Sg_limit  = 0.150;                                 % kg salt / kg solution
rho_brine = g(results,'rho_brine_out',1080);       % kg/m3, fallback 1080
M.C_limit_kgm3 = Sg_limit * rho_brine;
M.ok_salinity  = isfinite(M.TDS_traj_max) && M.TDS_traj_max <= M.C_limit_kgm3;

% BUILDABLE STACK HEIGHT. Section 4 lists this among the validity screens for
% the design-case selection but gives no numeric bound, so nothing there can
% be enforced. 2.0 m is a placeholder standing for "a unit a
% technician can service without a platform" -- SET IT DELIBERATELY or delete
% the screen, because an unstated limit that silently passes everything is
% worse than no limit. Deliberately OUTSIDE ok_all: exceeding it is an
% engineering-practicality judgement, not a model failure, and such runs
% belong on the curves with a note rather than being dropped.
H_stack_max = 2.0;                                    % [m] PLACEHOLDER
M.ok_stack  = ~isfinite(M.stack_h_m) || M.stack_h_m <= H_stack_max;

M.ok_all = M.ok_run && M.ok_solver && M.ok_massbal && M.ok_conserve && ...
           M.ok_recovery && M.ok_brine && M.ok_dewpoint && M.ok_preheatcap && ...
           M.ok_secband;
% ok_reynolds and ok_salinity are deliberately OUTSIDE ok_all: they flag
% correlation-validity limits, not model failure. Such runs stay in the
% curves, but the text must say the closure is used outside its fit.
end

function M = local_blank()
% THE SCHEMA. Every field local_extract() can set, in the order it sets them,
% followed by the six identity fields and two bookkeeping fields the runner
% attaches afterwards. This is the single declaration of the results-table
% columns: add a metric here and in local_extract(), nowhere else.
%
% Numeric defaults are NaN, not 0. A zero in an SEC column is a measurement;
% a NaN is an absence, and the two must not be confusable in a printed table.
% Screens default to FALSE: a run that did not complete has not passed
% anything, and a screen that defaults true would quietly admit it.
n = NaN; f = false; s = "";
M = struct( ...
    'SEC_external',n, 'R_still_pct',n, 'prod_kg_day',n, 'evap_kg_day',n, ...
    'Nx',n, 'A_plate_m2',n, 'A_total_m2',n, 'prod_L_m2_day',n, ...
    'Gamma_star',n, 'h_bar_m',n, 'stack_h_m',n, 'stack_basis',s, ...
    'R_plant_pct',n, 'SEC_blower',n, 'SEC_heater',n, 'SEC_compr',n, ...
    'SEC_airheat',n, 'W_external',n, 'SEC_band_lo',n, 'SEC_band_hi',n, ...
    'brine_TDS',n, 'CF_cascade',n, 'Tfeed_casc',n, 'f_cond',n, ...
    'regime_feedHX',s, 'Tv_exit',n, 'COP_awg',n, 'lift_awg',n, ...
    'residence_min',n, 't_end',n, ...
    'ok_run',f, 'ok_solver',f, 'ok_massbal',f, ...
    'cons_resid_max',n, 'cons_basis',s, 'ok_conserve',f, ...
    'recovery_ceil',n, 'ok_recovery',f, 'TDS_traj_max',n, 'ok_brine',f, ...
    'dew_margin_K',n, 'ok_dewpoint',f, 'dT_preheat',n, 'ok_preheatcap',f, ...
    'ok_secband',f, 'Re_gap_min',n, 'ok_reynolds',f, 'C_limit_kgm3',n, ...
    'ok_salinity',f, 'ok_stack',f, 'ok_all',f, ...
    'run_id',n, 'factor',s, 'fac_lbl',s, 'block',s, 'level',n, ...
    'is_base',f, 'wall_s',n, 'error',s);
end


function [r, basis] = conservation_residual(results)
% C2. THE CONSERVATION SCREEN, PER RUN.
%
% TESTING R_still_pct <= 100 ALONE IS A BOUND, NOT A BALANCE: it cannot fail
% for any physically ordered solution and therefore screens nothing. The model already closes four invariants at machine
% precision at the design point; the question a sensitivity study has to
% answer is whether that survives the PERTURBED corners -- small A_plate, low
% Qfan, Np = 20 at low feed -- which is exactly where it is most at risk.
%
% The returned quantity is the largest NORMALIZED residual across whatever
% invariant vector the model exposes. FIELD NAMES ARE A GUESS and must be
% confirmed against ZLDD_complete_modeling.m; if none match, basis says
% 'absent' and the screen is not enforced, because a screen that fails every
% run because it cannot find its input is worse than no screen.
r = NaN; basis = "absent";
if ~isstruct(results), return; end
for nm = {'conservation','invariants','cons','balance_residuals'}
    if isfield(results, nm{1}) && ~isempty(results.(nm{1}))
        v = results.(nm{1});
        if isstruct(v)
            fn = fieldnames(v); acc = [];
            for i = 1:numel(fn)
                x = v.(fn{i});
                if isnumeric(x), acc = [acc; abs(x(:))]; end %#ok<AGROW>
            end
            if ~isempty(acc), r = max(acc); basis = string(nm{1}); return; end
        elseif isnumeric(v)
            r = max(abs(v(:))); basis = string(nm{1}); return
        end
    end
end
end


function [H, basis] = stack_height(results, FC, h_bar)
% TOTAL ENCLOSURE STACK HEIGHT [m].
%
%   H = sum_j (V_gap,j / A_p)  +  Np*t_plate  +  t_base + t_ins
%
% WHY THE GAP VOLUMES AND NOT Ns*h_bar. Gap 1, the top compartment beneath
% the glazing, is NOT an inter-plate gap: it is set by the glazing tilt over
% the chamber length, not by hcomp + taper. At the design case its air
% inventory is 0.406 kg against 0.066 kg in each interior gap, so its volume
% is about SIX TIMES an interior gap. Multiplying the interior mean height by
% Ns therefore understates the stack substantially, and understates it MOST
% at small hcomp -- exactly the corner of the hcomp sweep where a spurious
% optimum is most likely. The sum over actual gap volumes has no such bias.
%
% The model already carries V_gap,j because the gap storage term needs
% M_v,j = rho_m * V_gap,j, so this is a lookup, not a re-derivation.
%
% FIELD NAMES ARE A GUESS AND MUST BE CONFIRMED against the model. If none
% match, the fallback below is used and 'basis' says so -- treat any run
% flagged 'fallback' as an ESTIMATE and do not quote it as a buildability
% constraint. The fallback cannot see the tall top compartment and is
% therefore a LOWER BOUND on the true stack.
A = gd(FC,'L',NaN) * gd(FC,'W',NaN);
Vg = [];
for f = {'V_gap','Vgap','V_gap_j','gap_volume'}
    if isstruct(results) && isfield(results,f{1}) && ~isempty(results.(f{1}))
        Vg = results.(f{1})(:); break
    end
end

t_plate = firstfield(FC, {'ts','t_plate','t_s','plate_thickness'}, 0.001);
t_base  = firstfield(FC, {'t_base','tbase'},                      0.005);
t_ins   = firstfield(FC, {'t_ins','tins'},                        0.050);
Np      = gd(FC,'Np',NaN);

% Guard against a per-node or per-timestep array being picked up instead of
% the per-gap vector: Ns = Np+1 entries expected. A length mismatch means the
% field name matched something else, and summing it would give a silently
% wrong height rather than an error.
if ~isempty(Vg) && isfinite(Np) && numel(Vg) ~= Np+1
    Vg = [];
end

if ~isempty(Vg) && isfinite(A) && A > 0
    H = sum(Vg)/A + Np*t_plate + t_base + t_ins;
    basis = "gap_volumes";
else
    % Lower bound: every gap treated as an interior gap. Misses the tall top
    % compartment entirely.
    H = (Np+1)*h_bar + Np*t_plate + t_base + t_ins;
    basis = "fallback_lower_bound";
end
end

function v = firstfield(S, names, dflt)
v = dflt;
for i = 1:numel(names)
    if isstruct(S) && isfield(S,names{i}) && ~isempty(S.(names{i}))
        v = S.(names{i}); return
    end
end
end

function v = gd(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function v = gn(s,p,d)
v = d;
for i = 1:numel(p)
    if ~isstruct(s) || ~isfield(s,p{i}), return; end
    s = s.(p{i});
end
if ~isempty(s), v = s; end
end

function s = local_flags(M)
f = {'ok_run','ok_solver','ok_massbal','ok_conserve','ok_recovery','ok_brine', ...
     'ok_dewpoint','ok_preheatcap','ok_secband','ok_reynolds','ok_salinity', ...
     'ok_stack'};
bad = {};
for i = 1:numel(f)
    if isfield(M,f{i}) && ~M.(f{i}), bad{end+1} = f{i}(4:end); end %#ok<AGROW>
end
if isempty(bad), s = '[ok]'; else, s = ['[' strjoin(bad,',') ']']; end
end


%%=========================================================================
%%  POST-PROCESSING: LOG-ELASTICITY
%%=========================================================================
function E = Sensitivity_elasticity(T)
% Local sensitivity index, dimensionless and unit-free:
%
%     S_local = d(ln Y)/d(ln X)  at the baseline, by central difference on
%               the two levels bracketing it.
%     S_fit   = the same derivative regressed over the WHOLE admissible
%               range, i.e. a range-averaged power-law exponent.
%
% C1. THESE ARE NOT THE SAME QUANTITY AND MUST NOT BE CONFLATED. Reporting
% the global fit alone under a header promising the central difference hides
% the distinction. For a monotone power law they agree and nothing is
% lost. For a response with an interior optimum they can differ in SIGN. Qfan
% is exactly that case -- the entire premise of F1 is a minimum -- and a
% straight line through a falling and a rising branch returns a slope near
% zero, which would have ranked the study's most important design lever as
% inert. R2_fit exposes the failure: when it is low, S_fit is not an
% elasticity and must not be quoted as one.
%
% WHY THIS AND NOT BAR LENGTH. A tornado ranks factors by the SPREAD their
% levels produce, which conflates two different things: how responsive the
% output is, and how wide a range the analyst happened to choose. Qfan spans
% x11.7 and Tfeed a factor of 1.15, so Qfan wins the bar chart before any
% physics is consulted. The elasticity removes the range from the comparison
% and leaves the responsiveness. Report BOTH: S_local for mechanism, bar
% length for achievable design leverage over the admissible range. They
% answer different questions and disagreeing is informative, not a problem.
%
% Absolute-temperature factors (T_coil, T_air_in, Tfeed_target, Tfeed, Ta) are
% NOT ratio-scaled -- ln(325 K) has no physical meaning -- so a semi-elasticity
% d(ln Y)/dX in K^-1 is reported for those instead and must not be plotted on
% the same axis as the ratio-scaled ones.
absT = {'T_coil','T_air_in','Tfeed_target','Tfeed','Ta'};
base = T(T.is_base,:);
if isempty(base), error('Sensitivity:noBaseline','No baseline row in T.'); end
Y0   = base.SEC_external(1);
facs = unique(T.factor(~T.is_base),'stable');

vt = {'string','string','double','double','double','double','logical','double'};
vn = {'factor','scale','S_local','S_fit','R2_fit','range_pct','monotonic','n_ok'};
E  = table('Size',[numel(facs) numel(vn)],'VariableTypes',vt,'VariableNames',vn);
% A6. table('Size',...) fills doubles with ZERO, so an inapplicable or
% uncomputable entry prints as 0.000, indistinguishable from a measured null
% response -- the absolute-temperature factors would all read S_ln = 0.
% Absence must look like absence.
for v = {'S_local','S_fit','R2_fit','range_pct'}, E.(v{1})(:) = NaN; end

for i = 1:numel(facs)
    % A1: the baseline point belongs to this factor's curve and is required
    % for any index claiming to be local to the baseline.
    S = withbase(T, facs(i));
    S = S(S.ok_all,:);
    E.factor(i) = facs(i); E.n_ok(i) = height(S);
    if height(S) < 2, continue; end

    [x,ix] = sort(S.level); y = S.SEC_external(ix);
    k = isfinite(x) & isfinite(y); x = x(k); y = y(k);
    if numel(x) < 2, continue; end

    E.range_pct(i)  = 100*(max(y)-min(y))/Y0;
    dy              = diff(y);
    E.monotonic(i)  = all(dy >= -1e-12) || all(dy <= 1e-12);

    isratio = ~ismember(facs(i), absT);
    X0 = factor_base(facs(i));

    if isratio
        E.scale(i) = "ratio";
        % A4. eps_recup includes a level of exactly ZERO -- the disabled
        % recuperator. log(0) = -Inf propagates through polyfit and returned
        % NaN for the one factor whose value the abstract most wants to
        % quote. Non-positive levels cannot appear in a ratio-scaled index at
        % all; they are a discrete configuration change, not a perturbation,
        % and are reported as a separate delta in the text.
        p = x > 0 & y > 0;
        u = log(x(p)); w = log(y(p));
        [E.S_local(i), E.S_fit(i), E.R2_fit(i)] = ...
            local_and_fit(u, w, log(max(X0,realmin)));
    else
        E.scale(i) = "absolute";
        % Semi-elasticity d(ln Y)/dX in K^-1. ln(325 K) has no meaning, so
        % the abscissa stays linear and this column must NEVER be plotted on
        % the same axis as the ratio-scaled one.
        p = y > 0;
        [E.S_local(i), E.S_fit(i), E.R2_fit(i)] = ...
            local_and_fit(x(p), log(y(p)), X0);
    end
end
E = sortrows(E,'range_pct','descend');

% ---- READING THIS TABLE. S_local and S_fit answer different questions and
% are expected to disagree; the disagreement is the information.
%   S_local  central difference across the two levels bracketing the
%            baseline. This is the LOCAL index the section header refers to.
%   S_fit    single slope regressed over the whole admissible range: a
%            range-AVERAGED exponent, valid only if the response really is a
%            power law over that range.
%   R2_fit   how badly that assumption fails. A low R2 means S_fit must not
%            be quoted as an elasticity at all.
%   monotonic  false means a single scalar slope is meaningless by
%            construction. Qfan is the case in point: the entire premise of
%            F1 is an interior minimum, so a straight line through a falling
%            and a rising branch returns a slope near zero and would rank
%            the most important design lever in the study as inert. S_fit
%            alone reports exactly that number and nothing else.
if any(~E.monotonic & E.n_ok >= 2)
    nm = strjoin(cellstr(E.factor(~E.monotonic & E.n_ok >= 2)), ', ');
    fprintf(['  NON-MONOTONIC over the swept range: %s. For these, S_fit is\n' ...
             '  not an elasticity -- quote S_local and the response curve.\n'], nm);
end
end


function [S_local, S_fit, R2] = local_and_fit(u, w, u0)
% Slope of w against u, twice: locally at u0 and as a regression over all of
% it. The local value is a central difference on the two abscissae bracketing
% u0, which is what "evaluated at the baseline" means; it falls back to a
% one-sided difference when the swept range lies entirely on one side of the
% base, and that asymmetry is precisely why the fallback must be visible in
% the code rather than hidden by a global fit.
S_local = NaN; S_fit = NaN; R2 = NaN;
if numel(u) < 2, return; end
[u,ix] = sort(u); w = w(ix);

lo = find(u < u0 - eps(u0)*8, 1, 'last');
hi = find(u > u0 + eps(u0)*8, 1, 'first');
if ~isempty(lo) && ~isempty(hi)
    S_local = (w(hi) - w(lo)) / (u(hi) - u(lo));          % central
else
    j = find(abs(u - u0) <= eps(u0)*8, 1);                % the base point
    if ~isempty(j)
        if ~isempty(hi),     S_local = (w(hi)-w(j))/(u(hi)-u(j));  % forward
        elseif ~isempty(lo), S_local = (w(j)-w(lo))/(u(j)-u(lo));  % backward
        end
    end
end

p     = polyfit(u, w, 1);
S_fit = p(1);
res   = w - polyval(p, u);
sst   = sum((w - mean(w)).^2);
if sst > 0, R2 = 1 - sum(res.^2)/sst; else, R2 = NaN; end
end


function Sensitivity_figs(T, pick)
%==========================================================================
%  ZLDD OAT -- FIGURE SUITE
%
%  USAGE
%    Sensitivity_figs(T);        % reached via Sensitivity('figs',T)
%
%  A SUBSET IS DRAWN BY NAME. Passing a tag or a list of tags restricts the
%  suite to those panels, so a single figure can be re-cut after an edit
%  without re-drawing and re-closing the other five:
%
%    Sensitivity('figs', T, 'F5')            one panel
%    Sensitivity('figs', T, {'F2a','F4'})    several
%    Sensitivity('figs', T)                  all of them
%
%  Tags are F1, F2a, F2b, F3, F4, F5 and are matched case-insensitively. F6
%  (validity) and the pinch plot are separate entry points, not part of this
%  suite.
%
%  Figures are drawn on screen only; export them from the caller if needed.
%
%  WHAT TO PUT IN THE PAPER, and why. Seven figures are built; you do not
%  want all seven in the manuscript. The recommended set for a Desalination
%  submission is F1 (mechanism), F2a (locates the optimum), F4 (trade-off)
%  and F6 (validity). F2b, F3 and F5 are supporting material.
%
%    F1  SEC DECOMPOSITION vs Qfan -- STACKED AREA.
%        The single strongest figure available from this sweep. A response
%        curve shows THAT a minimum exists; the stacked decomposition shows
%        WHY: blower work climbing as Q^3 from the right, heater and
%        compressor falling as mdot_da from the left, and the minimum sitting
%        where the two cross. A reviewer who sees the mechanism does not have
%        to take the optimum on trust. Nothing else in the sweep explains
%        itself this directly.
%
%    F2a RESPONSE CURVES, ratio-scaled factors. Y/Y0 against X/X0 on log-log.
%        Slope IS the elasticity, so the ranking is readable off the figure
%        without a separate tornado.
%
%    F2b RESPONSE CURVES, absolute-temperature factors, linear in kelvin.
%        SEPARATE FIGURE, not a companion panel. Putting T on a ratio axis
%        would be a category error, and sharing a frame with F2a implied the
%        two abscissae were comparable when only the ordinate is. Splitting
%        them also doubles the width available to each, which the earlier
%        two-panel layout did not leave for eight overlapping series.
%
%    F3  TORNADO, split by block. Design, scenario and uncertainty MUST NOT
%        share a chart: the first ranks design leverage, the third ranks
%        epistemic uncertainty, and stacking them invites the reader to
%        compare a lever you choose with a number you guessed. Three small
%        panels, one per block, sharing an x-axis.
%
%    F4  TRADE-OFF LOCUS: SEC_external against prod_L_m2_day, each factor a
%        connected path, baseline marked. This is the figure that earns the
%        OAT its place -- it shows that the factors do not merely move the
%        outputs, they move them in DIFFERENT DIRECTIONS in the performance
%        plane. Factors whose paths run along the Pareto front are real
%        design levers; factors whose paths run perpendicular to it trade one
%        objective for the other and need a decision, not an optimum.
%
%    F5  REGIME DIAGNOSTIC: Tfeed_casc against Tfeed_target with the 1:1 line.
%        Departure from 1:1 locates the binding-to-inert transition of the
%        feed exchanger. Turns "the bar is short" into "the unit saturates
%        at T = ...", which is a result rather than an absence of one.
%
%    F6  VALIDITY ENVELOPE: every run as a marker in (factor, level) space,
%        colored by which screen it trips. Publishing where the model stops
%        being trustworthy is what makes the rest of it credible, and it is
%        far more convincing as a figure than as the prose list the runner
%        currently prints.
%
%  WHAT THIS SWEEP CANNOT TELL YOU. OAT explores only the axes through one
%  baseline point. It cannot see interactions, and in n dimensions it samples
%  a vanishing fraction of the space -- with 16 factors the fraction of the
%  hypercube within reach of the axes is negligible. Two consequences worth
%  stating in the paper rather than leaving for a reviewer:
%    (i)  the "optimum" Qfan is conditional on every other factor sitting at
%         baseline, and will move once Np or T_air_in move;
%    (ii) any claim that a factor is UNIMPORTANT is the weakest claim OAT can
%         make, because a factor inert along the baseline axis can still be
%         active through an interaction.
%  If a stronger statement is needed, Morris screening costs roughly
%  (k+1)*r solves for k factors and r trajectories -- and its sigma
%  detects interaction and nonlinearity directly, without the cost of a full
%  Sobol decomposition. That matters most for F4: the whole claim of that
%  figure is that factors move the two objectives in DIFFERENT directions,
%  and OAT cannot say whether those directions survive off the baseline axes.
%==========================================================================

% ---- CANVAS GEOMETRY AND TYPOGRAPHY ------------------------------------
% Taken from plot_baseline.m so the sweep figures and the Results-section
% figures form one set on the page: same frame width, same type sizes.
% Authored at final printed size, so nothing is rescaled between MATLAB and
% the typeset page.
FIG_W  = 19.0;   % [cm] Elsevier double-column maximum, all wide figures
FIG_H  = 11.0;   % [cm] standard height for multi-row layouts
FIG_H1 =  8.0;   % [cm] one-row layouts: shorter frame keeps tile proportions
FIG_WN = 12.0;   % [cm] narrow frame for single-axes figures

ST = local_style();
FS_tick = ST.tick;   % tick labels
FS_lbl  = ST.lbl;    % axis labels
FS_ttl  = ST.ttl;    % panel titles
FS_leg  = ST.leg;    % legends and in-axes annotation
MS      = ST.ms;     % marker size, matched to the type size

% 'Arial' not 'Helvetica': Helvetica is absent on most Windows and Linux
% MATLAB installs and falls back silently, changing text extents between the
% screen figure and the printed file.
old_defaults = get(0, {'DefaultAxesFontName','DefaultAxesFontSize', ...
                       'DefaultTextFontName','DefaultLineLineWidth', ...
                       'DefaultAxesLabelFontSizeMultiplier', ...
                       'DefaultAxesTitleFontSizeMultiplier'});
set(0,'DefaultAxesFontName',ST.font,'DefaultAxesFontSize',FS_tick, ...
      'DefaultTextFontName',ST.font,'DefaultLineLineWidth',1.0, ...
      'DefaultAxesLabelFontSizeMultiplier',1.0, ...
      'DefaultAxesTitleFontSizeMultiplier',1.0);
cleanupDefaults = onCleanup(@() set(0, ...
    {'DefaultAxesFontName','DefaultAxesFontSize','DefaultTextFontName', ...
     'DefaultLineLineWidth','DefaultAxesLabelFontSizeMultiplier', ...
     'DefaultAxesTitleFontSizeMultiplier'}, old_defaults)); 

% ---- WHICH PANELS TO DRAW. An empty or missing selection draws the suite.
% want() is the single gate every block goes through, so a tag matching
% nothing raises rather than silently producing the full set and leaving the
% caller to work out which window is which.
tags = ["F1","F2a","F2b","F3","F4","F5"];
if nargin < 2 || isempty(pick)
    pick = tags;
elseif ischar(pick)
    pick = string(pick);          % 'F5' is one tag, not the characters F and 5
else
    pick = string(pick(:)).';     % cellstr or string array, any orientation
end
bad = pick(~ismember(lower(pick), lower(tags)));
if ~isempty(bad)
    error('Sensitivity_figs:unknownFigure', ...
          'No figure named %s. Valid tags: %s.', ...
          strjoin(cellstr(bad),', '), strjoin(cellstr(tags),', '));
end
want = @(t) ismember(lower(string(t)), lower(pick));
fprintf('  Figure suite: drawing %s\n', strjoin(cellstr(pick), ', '));

base = T(T.is_base,:);
if isempty(base), error('Sensitivity_figs:noBaseline','No baseline row in T.'); end

C = struct('blower',[0.35 0.35 0.38],'heater',[0.80 0.42 0.20], ...
           'compr',[0.25 0.50 0.70],'airheat',[0.70 0.25 0.35], ...
           'base',[0.10 0.10 0.10],'bad',[0.75 0.15 0.15]);
absT = {'T_coil','T_air_in','Tfeed_target','Tfeed','Ta'};

%% ---- F1  SEC decomposition vs Qfan (stacked area) --------------------
% A1: inject the shared baseline, which IS the Qfan = 0.25 point and would
% otherwise be missing from this curve. A5: the minimum must be located over VALID runs
% only -- the low-Qfan end is exactly where ok_reynolds is expected to trip,
% and an optimum sitting on a screened run is not an optimum.
S = withbase(T,"Qfan"); S = S(S.ok_all,:);
if want("F1") && height(S) > 2
    f1 = newfig('F1_SEC_decomposition_Qfan', FIG_W, FIG_H1);
    comp = [S.SEC_blower S.SEC_compr S.SEC_heater S.SEC_airheat];
    comp(~isfinite(comp)) = 0;
    ar = area(S.level, comp, 'LineStyle','none'); hold on
    ar(1).FaceColor = C.blower; ar(2).FaceColor = C.compr;
    ar(3).FaceColor = C.heater; ar(4).FaceColor = C.airheat;
    plot(S.level, S.SEC_external,'k-','LineWidth',2)
    [~,im] = min(S.SEC_external);
    plot(S.level(im), S.SEC_external(im),'kv','MarkerFaceColor','w', ...
         'MarkerSize',9,'LineWidth',1.4)
    text(S.level(im), S.SEC_external(im), ...
         sprintf('  min %.0f kWh m^{-3} at %.3f m^3 s^{-1}', ...
                 S.SEC_external(im), S.level(im)),'VerticalAlignment','bottom')
    % A3. The anchor comes from the factor declaration, like every other base
    % in this file. base.level is NaN by construction for the baseline row,
    % and NaN propagates through IEEE arithmetic, so anchoring the line to
    % that row draws no marker at all.
    xline(factor_base("Qfan"),'k:','baseline','LabelOrientation','horizontal');
    xlabel('Fan volumetric flow  \itQ\rm_{fan}  [m^3 s^{-1}]','FontSize',FS_lbl)
    ylabel('SEC_{external}  [kWh m^{-3}]','FontSize',FS_lbl)
    legend({'Blower','Compressor','Crystallizer heater','Air heater','Total'}, ...
           'Location','northwest','Box','off','FontSize',FS_leg)
    % Caption sentence, not a title: the minimum is where vapor-fired preheat
    % stops covering the air-heating demand, so it is a crossing of two terms
    % rather than a smooth optimum.
    local_title('Purchased-work decomposition against fan volumetric flow', ...
                'FontSize',FS_ttl,'FontWeight','normal')
    grid on; box off; set(gca,'Layer','top','XScale','log')
end

%% ---- F2a  normalized response, ratio-scaled factors --------------------
% SPLIT FROM THE ABSOLUTE-TEMPERATURE PANEL INTO ITS OWN FIGURE. The two
% panels share a y-axis and nothing else: this one is log-log, so its slope
% IS the elasticity and the ranking is readable straight off the figure;
% the other is linear in kelvin, because a ratio T/T0 on an absolute scale
% is a category error -- 330/325 is not "1.5 % more temperature" in any
% sense the model responds to. Side by side at half width each curve is
% roughly 8 cm across with up to eight overlapping series, and the legend
% covers the lower-left data. At full width both are legible, and each can
% be placed in the manuscript where its argument is made.
Y0  = base.SEC_external(1);
des = unique(T.factor(T.block=="design" & ~T.is_base),'stable');

if want("F2a")
f2a = newfig('F2a_response_ratio_scaled', 15.0, FIG_H); hold on
for i = 1:numel(des)
    if ismember(des(i),absT), continue; end
    S = withbase(T, des(i)); S = S(S.ok_all,:);
    if height(S)<2, continue; end
    X0 = factor_base(des(i));               % one declaration, not two
    k  = S.level > 0;                       % eps_recup = 0 is not a point
    if ~any(k), continue; end               % on a log abscissa
    plot(S.level(k)/X0, S.SEC_external(k)/Y0,'-o','LineWidth',1.4, ...
         'MarkerSize',MS,'DisplayName',char(S.fac_lbl(1)))
end
plot(1,1,'kp','MarkerFaceColor','k','MarkerSize',12,'DisplayName','Baseline')
set(gca,'XScale','log','YScale','log'); grid on; box off
xlabel('Factor level  \itX\rm / \itX\rm_0  [-]','FontSize',FS_lbl)
ylabel('SEC / SEC_0  [-]','FontSize',FS_lbl)
local_title('Baseline-normalized SEC: ratio-scaled factors  (slope = elasticity)', ...
      'FontSize',FS_ttl,'FontWeight','normal')
legend('Location','eastoutside','Box','off','FontSize',FS_leg)
end

%% ---- F2b  normalized response, absolute-temperature factors ------------
% Linear abscissa in kelvin. These factors have no meaningful ratio to the
% baseline, so they cannot go on the log panel above and their slopes are
% not elasticities.
if want("F2b")
f2b = newfig('F2b_response_absolute_temperature', 15.0, FIG_H); hold on
for i = 1:numel(des)
    if ~ismember(des(i),absT), continue; end
    S = withbase(T, des(i)); S = S(S.ok_all,:);
    if height(S)<2, continue; end
    plot(S.level, S.SEC_external/Y0,'-o','LineWidth',1.4,'MarkerSize',MS, ...
         'DisplayName',char(S.fac_lbl(1)))
end
yline(1,'k:','HandleVisibility','off'); grid on; box off
xlabel('Factor level  \itT\rm  [K]','FontSize',FS_lbl)
ylabel('SEC / SEC_0  [-]','FontSize',FS_lbl)
local_title('Baseline-normalized SEC: absolute-temperature factors', ...
      'FontSize',FS_ttl,'FontWeight','normal')
legend('Location','eastoutside','Box','off','FontSize',FS_leg)
end

%% ---- F3  tornado ------------------------------------------------------
% Drawn by the local zldd_tornado so that F3 and the cascade-recovery
% tornado are ONE figure style. Splitting design, scenario and uncertainty
% across three panels with three different x-scales makes the blocks
% incomparable by construction; instead the block is named beside each
% factor and every bar sits on a common axis.
%
% Bars are the deviation at the LOWEST and HIGHEST swept level, not the min
% and max of the response. The distinction matters for a non-monotonic factor
% such as Qfan: a min/max range bar discards the direction of the response
% and hides which end of the sweep produced it.
if want("F3"), zldd_tornado(T, 'SEC_external'); end

%% ---- F4  trade-off locus ----------------------------------------------
if want("F4")
f4 = newfig('F4_tradeoff_plane', FIG_WN, FIG_H); hold on
for i = 1:numel(des)
    S = withbase(T, des(i)); S = S(S.ok_all,:);   % A1: paths were broken at
    if height(S)<2, continue; end                 % the baseline point
    plot(S.prod_L_m2_day, S.SEC_external,'-o','LineWidth',1.3, ...
         'MarkerSize',MS,'DisplayName',char(S.fac_lbl(1)))
end
plot(base.prod_L_m2_day(1), base.SEC_external(1),'kp', ...
     'MarkerFaceColor','k','MarkerSize',14,'DisplayName','Baseline')
xlabel('Areal productivity  [L m^{-2} d^{-1}]','FontSize',FS_lbl)
ylabel('SEC_{external}  [kWh m^{-3}]','FontSize',FS_lbl)
local_title({'Performance plane: which factors are levers and which are trades', ...
       'down-right is unambiguously better'},'FontSize',FS_ttl,'FontWeight','normal')
legend('Location','eastoutside','Box','off','FontSize',FS_leg); grid on; box off
end

%% ---- F5  feed-exchanger regime ----------------------------------------
% The exchanger delivers the target until its duty saturates, after which the
% achieved feed temperature falls behind and further increases in the target
% buy nothing. The upper panel locates that departure against the 1:1 line;
% the lower panel plots the shortfall directly, because on the 1:1 axes a
% departure of a few kelvin over a 50 K span is a line thickness.
%
% The shortfall axis holds a floor of +/-1 K so that an exchanger which never
% saturates reads as a flat line on a scale with real units, rather than
% autoscaling zero noise into an apparent trend.
S = withbase(T,"Tfeed_target"); S = S(S.ok_all,:);
if want("F5") && height(S) > 2
    f5 = newfig('F5_feedHX_regime', 13.0, 12.0);
    tl = tiledlayout(f5,2,1,'TileSpacing','compact','Padding','compact');
    dS  = S.level - S.Tfeed_casc;          % shortfall [K], zero while target binds
    ib  = find(dS > 1, 1);
    lim = [min(S.level)-5 max(S.level)+5];

    ax1 = nexttile(tl); hold(ax1,'on')
    plot(ax1,lim,lim,'k:','LineWidth',1.0,'DisplayName','1:1 line (target delivered)')
    plot(ax1,S.level,S.Tfeed_casc,'-o','LineWidth',1.6,'MarkerSize',MS+1, ...
         'DisplayName','Achieved  \itT\rm_{feed,casc}')
    if ~isempty(ib)
        xline(ax1,S.level(ib),'r--',sprintf('departs 1:1 near %.0f K',S.level(ib)), ...
              'LabelOrientation','horizontal','FontSize',FS_leg,'HandleVisibility','off')
    end
    ylabel(ax1,'Achieved  \itT\rm_{feed,casc}  [K]','FontSize',FS_lbl)
    legend(ax1,'Location','northwest','Box','off','FontSize',FS_leg)
    grid(ax1,'on'); box(ax1,'off'); xlim(ax1,lim); ylim(ax1,lim)
    set(ax1,'XTickLabel',[])

    ax2 = nexttile(tl); hold(ax2,'on')
    yline(ax2,0,'k-','LineWidth',0.8,'HandleVisibility','off')
    plot(ax2,S.level,dS,'-s','LineWidth',1.6,'MarkerSize',MS+1,'Color',[0.80 0.33 0.20])
    if ~isempty(ib)
        xline(ax2,S.level(ib),'r--','HandleVisibility','off')
    else
        text(ax2,mean(lim),0.45,'target delivered at every swept level: shortfall = 0', ...
             'HorizontalAlignment','center','FontSize',FS_leg,'FontAngle','italic', ...
             'Color',[0.45 0.45 0.45])
    end
    xlabel(ax2,'Feed preheat target  \itT\rm_{feed,target}  [K]','FontSize',FS_lbl)
    ylabel(ax2,'Shortfall  \itT\rm_{target} - \itT\rm_{casc}  [K]','FontSize',FS_lbl)
    grid(ax2,'on'); box(ax2,'off'); xlim(ax2,lim)
    ylim(ax2,[-1 max(1.2, 1.15*max(dS))])

    local_title(tl,{'Feed exchanger: does the preheat target bind?', ...
              'a shortfall above zero marks a duty-limited exchanger, not an inert factor'}, ...
          'FontSize',FS_ttl,'FontWeight','normal')
end


end

%% ------------------------------------------------------------------ util
function X0 = factor_base(name)
% A2. THE BASELINE ANCHOR, READ FROM THE ONE DECLARATION.
%
% The base is read from local_factors() and never re-declared here. A second
% list written by hand drifts out of step: omit a ratio-scaled factor such as
% hcomp and it falls to the geometric mean of its levels (0.0186 m) instead of
% its declared base (0.020 m), so the hcomp curve in F2 normalizes to the
% wrong anchor and does not pass through unity while every other curve does.
% A duplicated constant is a defect waiting for an edit; exactly one place
% holds a base.
F = local_factors();
i = find(strcmp({F.name}, char(name)), 1);
if isempty(i)
    error('Sensitivity:unknownFactor', ...
          'No factor named %s in local_factors().', char(name));
end
X0 = F(i).base;
end


function u = factor_unit(name)
% THE UNIT, READ FROM THE SAME ONE DECLARATION as the base.
%
% local_factors() carries a unit for every factor but local_extract() never
% copied it into the results table, so the level columns of the tornado and
% the validity envelope printed bare numbers -- "305 / 330" says nothing
% without kelvin. Reading it here rather than re-declaring it keeps the one
% source of truth that factor_base already established.
F = local_factors();
i = find(strcmp({F.name}, char(name)), 1);
if isempty(i), u = ''; else, u = F(i).unit; end
if strcmp(u,'-'), u = ''; end        % dimensionless: print nothing, not "-"
end


function s = factor_axislabel(name, label)
% "Plate area (m^2)" -- ONE LINE, name and unit, nothing else.
%
% The block a factor belongs to (design / scenario / uncertainty) is a
% property of the STUDY DESIGN, not of the quantity, and it meant nothing to
% a reader of the figure; it also forced a second text line per row, which
% pushed the label off the row center and made the tall rows read as
% misaligned. The block still governs which factors appear together and is
% stated in the text, where it belongs.
u = factor_unit(name);
if isempty(u), s = char(label); else, s = sprintf('%s (%s)', char(label), u); end
end


function T = refresh_labels(T)
%REFRESH_LABELS  Take fac_lbl from the factor declaration, not from the file.
%
% ONE DECLARATION, INCLUDING THE WORDING. local_extract() copies the label into
% every results row at solve time, so a checkpoint carries the wording that
% local_factors() holds at the moment of the sweep. Every figure reads fac_lbl
% from the table, which makes the stored label the thing the reader sees: an
% edit to local_factors() moves the declaration and leaves the artwork behind,
% and the two disagree with no error anywhere.
%
% The mapping is by factor NAME, which is a model identifier and does not
% change with wording, so this re-reads text only and touches no numbers. The
% baseline row carries no factor and is passed through.
if ~ismember('fac_lbl', T.Properties.VariableNames), return; end
F = local_factors();
for i = 1:numel(F)
    m = T.factor == string(F(i).name);
    if any(m), T.fac_lbl(m) = string(F(i).label); end
end
end


function S = withbase(T, name)
% A1. RE-ATTACH THE SHARED BASELINE TO A FACTOR'S OWN CURVE.
%
% Run construction skips any level equal to the factor base, because the
% baseline is solved ONCE as run 1 and shared -- correct, and it saves one
% solve per factor. But that row carries factor = "BASELINE", so a filter of
% the form T.factor == <name> excludes it, and the consequence is not
% cosmetic:
%   - F1's stacked area has a hole at Qfan = 0.25, the baseline itself;
%   - F2's normalized curves do not pass through (1,1);
%   - F4's paths break at the baseline point;
%   - the elasticity fit never sees (X0, Y0), so a "local" index at the
%     baseline is computed from data that excludes it.
% Every factor whose level list contains its own base is affected.
%
% The fix is injection, not re-solving: copy the baseline row, relabel it to
% this factor and place it at the factor's declared base level. The physics
% is identical by construction, since run 1 IS this factor at its base.
S = T(T.factor == string(name) & ~T.is_base, :);
b = T(T.is_base, :);
if isempty(b) || isempty(S), S = sortrows(S,'level'); return; end
b = b(1,:);
b.factor  = string(name);
b.fac_lbl = S.fac_lbl(1);
b.block   = S.block(1);
b.level   = factor_base(name);
S = sortrows([S; b], 'level');
end

function x = normlev(v)
if max(v)==min(v), x = 0.5*ones(size(v)); else, x = (v-min(v))/(max(v)-min(v)); end
end

function h = newfig(name, w_cm, h_cm)
% One creation point so every figure is authored at final printed size in
% centimeters. Position in PIXELS would make the printed size depend on the
% display, which is what breaks a figure set across machines.
h = figure('Name',name,'Color','w','Units','centimeters', ...
           'Position',[2 2 w_cm h_cm]);
end


function D = zldd_tornado(T, resp, bases, cutoff)
%ZLDD_TORNADO  Two-sided tornado on any response column of the results table.
%
%   Reached from the command line through the Sensitivity entry points:
%       Sensitivity('tornado', T)                    % SEC_external
%       Sensitivity('tornado', T, 'R_still_pct')     % cascade recovery
%
%   Each factor gets TWO bars from a common zero at the baseline response:
%   one for the lowest swept level, one for the highest. Bars are colored by
%   the DIRECTION THE FACTOR WAS MOVED, not by the sign of the response --
%   the bar's own side already carries the sign, so coloring by it wastes a
%   channel. Coloring by factor direction makes an inverse response visible
%   at a glance (an orange bar pointing left) and puts both bars of a
%   non-monotonic factor on the same side, where a range bar would hide it.
%
%   ONE FIGURE STYLE FOR BOTH RESPONSES. SEC and cascade recovery are drawn by
%   this same routine so F3 and the recovery tornado are directly comparable.
%   Splitting design, scenario and uncertainty over three panels with three
%   different x-scales makes the blocks incomparable by construction; the
%   block is named beside each factor instead.
%
%   THE FULL SWEPT RANGE IS SHOWN, screened runs included. Filtering on ok_all
%   shortens a factor's range to whatever happened to pass, which misreports
%   the sweep as narrower than it is. Screen status stays in the returned table
%   as ok_lo and ok_hi, NaN where that side holds no runs; the validity
%   envelope is F6's job.
%
%   A "--" LEVEL MARKS A SIDE THAT HOLDS NO RUNS. A level list whose extreme
%   entry is its own base has no runs beyond it, so that column is empty by
%   construction. Those rows carry § and are named on the console.
%
%   THE OPTIMUM SENSE IS DECLARED PER RESPONSE in the switch below: minimum for
%   SEC, maximum for recovery and productivity.
%
%   INTERIORITY IS A POSITION IN THE LEVEL LIST. The optimum is taken over the
%   swept levels together with the base, and is flagged when it falls strictly
%   between the first and last of them.
%
%   RANGES ARE NOT COMMENSURATE unless you make them so. Bar length is set by
%   how far the factor was swept, so this ranks leverage only if the level
%   lists represent comparable intervals. eps_recup swept to zero is a
%   no-recuperator case, not a design range, and will dominate any SEC tornado
%   for that reason alone. State the convention in the caption.

if nargin < 2 || isempty(resp), resp = 'R_still_pct'; end
if nargin < 4 || isempty(cutoff), cutoff = NaN; end   % NaN -> auto
if ~ismember(resp, T.Properties.VariableNames)
    error('zldd_tornado:noResponse','No column "%s" in the table.', resp);
end

% ---- factor baselines. The baseline level is not in the level list (runs at
% the base are skipped in construction), so it cannot be recovered from T and
% is read from local_factors(), the single declaration the rest of the file
% uses. xbase sets the low/high split, the deviation columns and the
% interior-optimum candidate set together, so it is not restated locally.
% Override by passing a struct with the same field names.
if nargin < 3 || isempty(bases)
    F_ = local_factors();  bases = struct();
    for q_ = 1:numel(F_), bases.(F_(q_).name) = F_(q_).base; end
end

base_row = T(T.is_base,:);
if height(base_row) ~= 1
    error('zldd_tornado:baseline','Table has %d baseline rows, need 1.', height(base_row));
end
y0 = base_row.(resp);

% THE FULL SWEPT RANGE IS SHOWN, screened runs included. Filtering on ok_all
% silently shortens a factor's range to whatever happened to pass, which
% misreports the sweep as narrower than it was -- T_air_in loses its 305 and
% 330 levels that way and appears one-sided when it is not.
V    = T(~T.is_base, :);
facs = unique(V.factor,'stable');

name=strings(0); dlo=[]; dhi=[]; llo=[]; lhi=[]; oklo=[]; okhi=[]; xbase=[];
blk = strings(0);  flb = strings(0);
for i = 1:numel(facs)
    f = facs(i);
    if ~isfield(bases, char(f))
        warning('zldd_tornado:noBase','No baseline declared for %s -- skipped.', f);
        continue
    end
    S  = sortrows(V(V.factor == f, :), 'level');
    xb = bases.(char(f));
    below = S(S.level < xb, :);   above = S(S.level > xb, :);
    if isempty(below) && isempty(above), continue, end

    % ok is NaN where a side holds no runs. A level list that does not straddle
    % its base has no extreme on one side, and NaN records the absence rather
    % than a screen verdict on a run that does not exist.
    if isempty(below), d1=NaN; l1=NaN; k1=NaN;
    else, d1 = below.(resp)(1)   - y0;  l1 = below.level(1);   k1 = below.ok_all(1);
    end
    if isempty(above), d2=NaN; l2=NaN; k2=NaN;
    else, d2 = above.(resp)(end) - y0;  l2 = above.level(end); k2 = above.ok_all(end);
    end

    name(end+1,1)=f;  xbase(end+1,1)=xb;                        %#ok<AGROW>
    blk(end+1,1) = S.block(1);   flb(end+1,1) = S.fac_lbl(1);   %#ok<AGROW>
    dlo(end+1,1)=d1;  llo(end+1,1)=l1;  oklo(end+1,1)=k1;       %#ok<AGROW>
    dhi(end+1,1)=d2;  lhi(end+1,1)=l2;  okhi(end+1,1)=k2;       %#ok<AGROW>
end

% DIRECTION OF IMPROVEMENT, declared per response. A minimum is the target for
% SEC and a maximum for recovery and productivity, and the two cannot be
% distinguished from the response values alone. An undeclared response warns.
switch resp
    case {'SEC_external','SEC_still','SEC'},              sense = 'min';
    case {'R_still_pct','R_plant_pct','prod_L_m2_day'},   sense = 'max';
    otherwise
        warning('zldd_tornado:sense', ...
            'No optimisation sense declared for "%s"; using min.', resp);
        sense = 'min';
end

% INTERIOR OPTIMUM, recorded per factor. A tornado bar is read as the range a
% factor can move the response over, which holds only where the response is
% monotone in that factor. A response that turns inside the swept window
% reaches its best value at neither extreme, so the bar pair understates the
% factor and the turning point is marked on the panel.
%
% The candidate set is the swept levels together with the base: the base is an
% operating point like any other, and without it the level adjacent to it sits
% at the edge of the set, where a monotone response places the optimum.
%
% Interiority is a property of position in the sorted level list, so it is
% tested there. A test on response values needs a tolerance in the units of
% whichever column is plotted, and returns true for any factor with a NaN at
% one extreme. At least three surviving points are needed for an interior
% point to exist.
%
% Candidates are the full swept range, not the ok_all subset, matching the
% runs the bars are drawn from.
optlev = nan(numel(name),1);  optval = optlev;
interior = false(numel(name),1);
for i = 1:numel(name)
    S   = V(V.factor == name(i), :);
    lev = [S.level;      xbase(i)];
    val = [S.(resp);     y0     ];
    keep = ~isnan(val) & ~isnan(lev);
    lev = lev(keep);  val = val(keep);
    [lev, ix] = sort(lev);  val = val(ix);
    [lev, iu] = unique(lev, 'stable');  val = val(iu);   % base may repeat a level
    if numel(lev) < 3, continue, end
    if strcmp(sense,'min'), [mv,mi] = min(val); else, [mv,mi] = max(val); end
    optlev(i) = lev(mi);  optval(i) = mv - y0;
    interior(i) = (mi > 1) && (mi < numel(lev));
end

% SPAN IS THE LARGER ONE-SIDED EXCURSION, not the low-to-high spread. The two
% coincide only for an opposite-sign pair. Rows are ranked on it, so the
% convention belongs in the caption.
span = max([abs(dlo) abs(dhi)],[],2,'omitnan');
D = table(name,flb,blk,xbase,llo,dlo,oklo,lhi,dhi,okhi,span,optlev,optval, ...
     interior, ...
     'VariableNames',{'factor','label','block','base','level_lo','d_lo', ...
                      'ok_lo','level_hi','d_hi','ok_hi','span', ...
                      'opt_level','opt_delta','interior'});
D = sortrows(D,'span','descend');
D.one_sided = isnan(D.d_lo) | isnan(D.d_hi);
D.nonmono   = ~D.one_sided & sign(D.d_lo) == sign(D.d_hi) & D.d_lo ~= 0;

% PERTURBATION FROM THE BASE. Bar length is a leverage ranking only where the
% factors are perturbed comparably, which they are not here: Qfan runs
% -88 %/+20 % about its base while Np runs -50 %/+50 %. The perturbation is
% printed beside each level, so a bar long only because its factor spans a
% wider range can be discounted.
%
% Absolute temperatures are reported in kelvin. A fraction of the base
% measures a perturbation only on a ratio scale, and the factor notes in
% local_factors() sweep T_coil and T_air_in over admissible ranges for that
% reason. The absT flag selects the format per factor; lo_pct and hi_pct carry
% NaN for those, and lo_dev and hi_dev carry the kelvin deviation.
isT = false(height(D),1);
for i = 1:height(D), isT(i) = strcmp(factor_unit(D.factor(i)),'K'); end
D.absT   = isT;
D.lo_pct = 100*(D.level_lo - D.base)./D.base;   D.lo_pct(isT) = NaN;
D.hi_pct = 100*(D.level_hi - D.base)./D.base;   D.hi_pct(isT) = NaN;
D.lo_dev = D.level_lo - D.base;                 % K for absT, raw units otherwise
D.hi_dev = D.level_hi - D.base;

% ---- INERT FACTORS GO TO A FOOTNOTE, NOT A ROW ------------------------
% A factor whose response does not move still needs reporting -- inertness is
% a result, and a reader who cannot find eps_recup will assume it was never
% swept. But five near-empty rows cost a third of the panel height and shrink
% every bar that does carry information, so they are named and quantified
% beneath the axis instead of drawn.
%
% The threshold is RELATIVE to the largest bar (2 % by default), because an
% absolute one in kWh/m3 would be meaningless on R_still and vice versa. Pass
% cutoff = 0 to draw every factor.
if isnan(cutoff), cutoff = 0.02 * max(D.span); end
inert = D.span < cutoff;
Dfull = D;                       % returned in full; only the plot is trimmed
D     = D(~inert,:);
Dout  = sortrows(Dfull(inert,:),'span','descend');

% ---- draw ---------------------------------------------------------------
% BOTH BARS SHARE ONE LINE and both start at zero, which is what makes this a
% tornado rather than a paired bar chart. The longer is drawn first so that a
% same-sign pair -- a non-monotonic factor -- nests instead of hiding.
%
% THREE LEVEL COLUMNS: low level at the left margin, high level at the right,
% baseline boxed on the zero line. Bar length alone does not say what the
% factor was set to, and a reader cannot check a range that is not printed.
%
% ALL BARS ARE DRAWN SOLID. Validity is reported by the envelope figure and
% by ok_all in the results table, so it is not annotated a second time here.
% The screen status of each extreme is still returned in D (ok_lo, ok_hi) for
% anything that needs it.
ST_ = local_style();          % one type preset for the whole file
n  = height(D);
yy = (n:-1:1).';
C_dn = [0.20 0.45 0.70];
C_up = [0.80 0.35 0.20];
C_op = [0.10 0.45 0.25];   % interior optimum marker
xm   = max(abs([D.d_lo; D.d_hi; D.opt_delta]),[],'omitnan');
eps_draw = 0.004 * xm;     % below this a bar is thinner than its own edge
LM   = 3.90;   % left margin, in units of xm. The descriptive labels are long
               % ('Still inlet air temperature'), and at 1.7 they overprint
               % the low-level column. The column is in axis units while the
               % type is in points, so the width needed here scales inversely
               % with the canvas: a narrower canvas needs a larger LM.

% The figure Name becomes the export filename, so it has to identify the
% response as well as the figure type.
switch resp
    case 'SEC_external', figname = 'F3_tornado_SEC';
    case 'R_still_pct',  figname = 'F7_tornado_cascade_recovery';
    otherwise,           figname = ['F_tornado_' resp];
end
% CANVAS AUTHORED AT PRINTED WIDTH. This figure carries more small type than
% any other in the file -- two level columns, a value on every bar and a
% footnote block -- so it is the one most exposed to placement scaling. The
% layout is entirely in axis units, so the canvas sets the type size relative
% to the artwork and nothing else: height follows width to hold the aspect
% ratio, and the bar values reach the page at their authored 7.5 pt.
TW = 17.0;                      % [cm] canvas width, near the printed width
figure('Name',figname,'Color','w', ...
       'Units','centimeters','Position',[2 2 TW (0.95*n+3.75)]);
ax = axes; hold(ax,'on'); box(ax,'off'); grid(ax,'on')
ax.YAxis.TickValues = []; ax.Toolbar = [];

% Alternating row bands. With sixteen unlabelled rows the eye loses which
% level column belongs to which bar; a faint band is cheaper than gridlines
% and does not compete with the bars.
for i = 1:2:n
    patch([-9e9 9e9 9e9 -9e9],[yy(i)-0.5 yy(i)-0.5 yy(i)+0.5 yy(i)+0.5], ...
          [0.955 0.955 0.955],'EdgeColor','none','HandleVisibility','off');
end

% SAME-SIGN PAIRS GET THEIR OWN HALF-ROW EACH; opposite-sign pairs share the
% full row. Drawing both from zero on one line is the tornado's whole point
% when they diverge, but when BOTH extremes move the response the same way
% -- a non-monotonic factor, or one whose admissible levels lie on one side
% of the base -- the shorter bar is drawn on top of the longer one and the
% figure silently reports one number where there are two. A_plate is the
% worst case here: +19.5 and +18.4 differ by 6 %, so the nested bar was
% invisible. Splitting the row costs a little height and hides nothing.
for i = 1:n
    dl = D.d_lo(i);  dh = D.d_hi(i);
    same = ~isnan(dl) && ~isnan(dh) && sign(dl) == sign(dh);
    if same
        dd = [dl dh];  cc = {C_dn, C_up};
        yc = [yy(i)+0.16, yy(i)-0.16];  hh = [0.14 0.14];
    else
        dd = [dl dh];  cc = {C_dn, C_up};
        yc = [yy(i) yy(i)];             hh = [0.30 0.30];
    end
    % Bars narrower than eps_draw are below the width at which a rectangle is
    % distinguishable from its own edge, and the row is labelled instead. The
    % threshold scales with the panel because the response carries different
    % units on the SEC and recovery tornadoes.
    if all(isnan(dd)) || max(abs(dd),[],'omitnan') < eps_draw
        text(0.02*xm, yy(i),'no effect','FontSize',ST_.ann, ...
             'Color',[.55 .55 .55],'FontAngle','italic');
        continue
    end
    % longest first so that, in the shared-row case, a shorter opposite bar
    % is never overdrawn
    [~,ord] = sort(abs(dd),'descend');
    for k = ord
        if isnan(dd(k)), continue, end
        rectangle('Position',[min(0,dd(k)) yc(k)-hh(k) abs(dd(k)) 2*hh(k)], ...
                  'FaceColor',cc{k},'EdgeColor',cc{k});
    end
    if D.interior(i)
        plot(D.opt_delta(i), yy(i),'d','MarkerSize',7, ...
             'MarkerFaceColor',C_op,'MarkerEdgeColor','w','LineWidth',0.9);
        text(D.opt_delta(i), yy(i)-0.36, sprintf('%g',D.opt_level(i)), ...
             'HorizontalAlignment','center','FontSize',ST_.ann,'Color',C_op, ...
             'FontWeight','bold');
    end
    for k = 1:2
        if isnan(dd(k)) || abs(dd(k)) < eps_draw, continue, end
        inside = abs(dd(k)) > 0.15*xm;
        if inside
            text(dd(k)-sign(dd(k))*0.012*xm, yc(k), sprintf('%+.1f',dd(k)), ...
                 'Color','w','FontWeight','bold','FontSize',ST_.ann, ...
                 'HorizontalAlignment',ternary(dd(k)>0,'right','left'));
        else
            text(dd(k)+sign(dd(k))*0.015*xm, yc(k), sprintf('%+.1f',dd(k)), ...
                 'Color',cc{k},'FontWeight','bold','FontSize',ST_.ann, ...
                 'HorizontalAlignment',ternary(dd(k)>0,'left','right'));
        end
    end
end
xline(0,'k-','LineWidth',1.2,'HandleVisibility','off');

for i = 1:n
    % THE BLOCK STAYS VISIBLE WITHOUT SPLITTING THE CHART. Design levers,
    % scenario conditions and uncertainty parameters answer different
    % questions and should not be ranked against one another as if they were
    % interchangeable -- but three separate panels made that point at the
    % cost of three inconsistent axes. One chart with the block named beside
    % each factor keeps the distinction and the common scale.
    % NO MATLAB IDENTIFIERS ON THE FIGURE. 'Qfan' and 'hcomp' are variable
    % names internal to the model; a reader of the paper has never seen them,
    % and printing them makes a published figure look like console output.
    % The slot is worth more as the UNIT, which the figure otherwise lacks
    % entirely -- the LOW and HIGH columns are bare numbers, so without it
    % "305 / 330" is unreadable. Traceability to the code belongs in the
    % supplementary table, not in the artwork.
    % White backing on the text: with the grid extending past the bars,
    % a gridline would otherwise strike through the labels.
    lbl_i = factor_axislabel(D.factor(i), D.label(i));
    if D.interior(i),   lbl_i = [lbl_i '  *'];  end    %#ok<AGROW>
    if D.one_sided(i),  lbl_i = [lbl_i '  §'];  end    %#ok<AGROW>
    text(-LM*xm, yy(i), lbl_i, ...
         'HorizontalAlignment','left','FontSize',ST_.lbl,'VerticalAlignment','middle', ...
         'BackgroundColor','w','Margin',0.5);
    text(-1.13*xm, yy(i), levcol(D,i,'lo'), ...
         'HorizontalAlignment','right','FontSize',ST_.tick,'Color',C_dn, ...
         'FontWeight','bold','Interpreter','none','BackgroundColor','w','Margin',0.5);
    text( 1.13*xm, yy(i), levcol(D,i,'hi'), ...
         'HorizontalAlignment','left','FontSize',ST_.tick,'Color',C_up, ...
         'FontWeight','bold','Interpreter','none','BackgroundColor','w','Margin',0.5);
end
% COLUMN HEADINGS IN SENTENCE CASE AND IN FULL. Capitals read as console
% output in a typeset figure, and an abbreviated heading spends the caption on
% defining the artwork instead of on the result.
text(-LM*xm, n+0.85,'Factor','HorizontalAlignment','left','FontSize',ST_.tick, ...
     'Color',[.3 .3 .3],'FontWeight','bold');
text(-1.13*xm, n+0.85,'Low level (deviation from baseline)', ...
     'HorizontalAlignment','right','FontSize',ST_.tick, ...
     'Color',C_dn,'FontWeight','bold');
text( 1.13*xm, n+0.85,'High level (deviation from baseline)', ...
     'HorizontalAlignment','left','FontSize',ST_.tick, ...
     'Color',C_up,'FontWeight','bold');

xlim([-(LM+0.06)*xm 2.00*xm]); ylim([0.35 n+1.4]);

% TICKS AND GRID ARE CLIPPED, BUT SYMMETRICALLY ABOUT ZERO. The label column
% occupies negative x INSIDE the axes -- it is not a margin -- so MATLAB
% otherwise keeps ticking out to the axis limit and draws a gridline under
% the factor names at an x no bar reaches. Clipping to the data span alone
% left an asymmetric grid, which on a diverging chart misreads as a shifted
% zero. The window is taken one round tick beyond the longest bar on BOTH
% sides, so the grid is symmetric and the reader can measure either way.
tk = get(gca,'XTick');
set(gca,'XTick', tk(abs(tk) <= 1.80*xm));

% The y-axis line is dropped with it. There is no y quantity here -- the rows
% are categories -- so a vertical rule at the far left reads as an axis that
% means something.
set(gca,'YColor','none');
xlabel(resp_axislabel(resp, y0),'FontSize',ST_.lbl);

% THE OMITTED FACTORS GO TO THE CONSOLE, NOT ONTO THE AXIS. Six names and
% six numbers set below the xlabel render at roughly 6 pt at print width,
% collide with the axis and read as a stray annotation. The caption has room
% for a proper sentence, so the list is printed for pasting there instead.
if ~isempty(Dout)
    parts = strings(height(Dout),1);
    for i = 1:height(Dout)
        parts(i) = sprintf('%s (%.2f)', Dout.label(i), Dout.span(i));
    end
    fprintf(['  FOR THE CAPTION -- omitted, |delta| < %.3g over their swept\n' ...
             '  range: %s.\n'], cutoff, strjoin(cellstr(parts), ', '));
end
% The marker key is listed only where a marker is drawn.
p1 = patch(NaN,NaN,C_dn,'EdgeColor','none');
p2 = patch(NaN,NaN,C_up,'EdgeColor','none');
hL = [p1 p2];
sL = {'Low level','High level'};
if any(D.interior)
    p3 = plot(NaN,NaN,'d','MarkerSize',7,'MarkerFaceColor',C_op, ...
              'MarkerEdgeColor','w','LineWidth',0.9);
    hL(end+1) = p3;
    sL{end+1} = sprintf('Interior %simum (level marked)', ...
                        ternary(strcmp(sense,'min'),'min','max'));
end
legend(hL, sL, 'Orientation','horizontal','Location','northoutside','Box','off', ...
       'FontSize',ST_.leg,'Interpreter','none');

% FOOTNOTES SIT IN A STRIP BELOW THE AXIS, DIRECTLY UNDER THE XLABEL. Inside
% the axes they overlap the bottom factor row, and a text object below ylim
% does not clip by default, so it lands on the xlabel instead.
%
% The block is anchored to the MEASURED BOTTOM OF THE XLABEL, not to the
% bottom of the frame. The axes already carries a bottom margin sized for the
% tick labels and the xlabel, so anchoring to the frame leaves that whole
% margin as white space between the xlabel and the footnotes.
fn = strings(0);
if any(D.interior)
    fn(end+1) = sprintf(['*  %simum lies strictly inside the swept range, so ' ...
        'neither extreme reaches it'], ternary(strcmp(sense,'min'),'min','max'));
end
if any(D.one_sided)
    fn(end+1) = "§  levels lie on one side of the base only, so this row is not a two-sided bar";
end
fn(end+1) = "deviations from baseline in parentheses: absolute temperatures in K, all other factors in %";

% ---- FOOTNOTE APPEARANCE. The four knobs are here rather than inline, so a
% change of look is one edit and not a hunt through the annotation call.
FN_COL = [0.20 0.20 0.20];   % text color; darker reads better in print
FN_FS  = ST_.ann;            % [pt] font size, from the shared preset
LH     = 0.40;               % [cm] line pitch; raise to open the block up
GAP    = 0.10;               % [cm] xlabel to the first footnote line; lower
                             % to move the whole block up, raise to drop it
DX     = 0.00;               % [cm] horizontal shift, + moves right

drawnow
fh = ancestor(ax,'figure');
fh.Units = 'centimeters';  ax.Units = 'centimeters';
p0 = ax.Position;                              % capture BEFORE the frame grows

% How far the block reaches below the xlabel, and how much of that the
% existing bottom margin already provides.
blk  = GAP + LH*numel(fn) + 0.15;
xl   = ax.XLabel;  xl.Units = 'centimeters';
e = get(xl,'Extent');                          % [cm], relative to the axes
if numel(e) == 4 && isfinite(e(2))
    xlab_bot = p0(2) + e(2);                   % [cm] from the frame bottom
else
    xlab_bot = p0(2) - 1.10;                   % fallback if Extent is unset
end
strip = max(0, blk - xlab_bot);                % grow only by what is missing

fh.Position(4) = fh.Position(4) + strip;
ax.Position    = [p0(1) p0(2)+strip p0(3) p0(4)];
drawnow
e = get(xl,'Extent');                          % re-measure after the move
if numel(e) == 4 && isfinite(e(2))
    top = p0(2) + strip + e(2);
else
    top = xlab_bot + strip;
end
fp = fh.Position;
for q = 1:numel(fn)
    yq = top - GAP - q*LH;                     % [cm] from the frame bottom
    xq = p0(1) + DX;                           % [cm] from the frame left edge
    annotation(fh,'textbox', [xq/fp(3) yq/fp(4) 1-xq/fp(3) LH/fp(4)], ...
        'String', char(fn(q)), 'EdgeColor','none', 'Margin',0, ...
        'FontName','Arial','FontSize',FN_FS,'FontAngle','italic', ...
        'Color',FN_COL,'VerticalAlignment','middle', ...
        'FitBoxToText','off','Interpreter','none');
end
xl.Units = 'data';  ax.Units = 'normalized';  fh.Units = 'centimeters';

set(ax,'Layer','top','FontName','Arial','FontSize',8.5);

D = Dfull;                       % return every factor, trimmed or not
fprintf('  Baseline %s = %.3f  |  %d of %d factors drawn (cutoff %.3g)\n', ...
        resp, y0, n, height(Dfull), cutoff);
% Structural notes on the sweep, one line each. interior is the flag the panel
% marks, so the listing and the figure name the same factors.
lst = @(t) strjoin(cellstr(t.label(:)), ', ');
if any(Dfull.interior)
    fprintf('  Interior %simum (marked): %s\n', sense, lst(Dfull(Dfull.interior,:)));
end
if any(Dfull.nonmono)
    fprintf('  Same-sign extremes: %s\n', lst(Dfull(Dfull.nonmono,:)));
end
if any(Dfull.one_sided)
    fprintf('  Levels do not straddle base: %s\n', lst(Dfull(Dfull.one_sided,:)));
end
end


function export_figs(fmt, figs)
%EXPORT_FIGS  Write every open figure at its authored size, as vector.
%
%   Sensitivity('export')          % EPS, the Elsevier vector format
%   Sensitivity('export','svg')    % SVG, for further editing
%
%   THE PAINTERS RENDERER IS FORCED. With OpenGL the output is a raster image
%   wrapped in a vector container -- the file carries a vector extension and
%   the content does not scale, which is the commonest way an export is not
%   actually vector. Painters is what emits paths.
%
%   PRINT, NOT EXPORTGRAPHICS. exportgraphics tight-crops to the drawn
%   content, so the saved aspect ratio comes from the bounding box rather than
%   the authored frame, and figures built on a common canvas land on the page
%   at different sizes. PaperSize is set equal to PaperPosition so there is no
%   page margin and no rescaling.
%
%   FILENAMES COME FROM THE FIGURE Name PROPERTY, not the figure number:
%   numbers depend on the order the figures happened to open, so a rerun in a
%   different order would overwrite the wrong file.
%
%   A stacked-area panel is the one case where vector can disappoint --
%   painters antialiases each patch path separately and can leave hairline
%   seams between the bands of F1. Check that one at high zoom; if the seams
%   show, re-export F1 alone as 600 dpi raster.
if nargin < 1 || isempty(fmt),  fmt = 'eps'; end
if nargin < 2, figs = []; end

% A SPECIFIC SET OF FIGURES CAN BE NAMED. Exporting everything open is the
% right default when the whole suite has just been drawn, but a Morris run
% that follows a sweep would otherwise re-export the sweep's figures and
% overwrite files that were already correct.
if isempty(figs), figs = findobj(groot,'Type','figure'); end
figs = figs(isgraphics(figs,'figure'));
if isempty(figs), warning('Sensitivity:export','No open figures.'); return, end
[~,ord] = sort([figs.Number]);  figs = figs(ord);

fprintf('  Exporting %d figure(s) as %s ...\n', numel(figs), upper(fmt));
for i = 1:numel(figs)
    h  = figs(i);
    nm = h.Name;  if isempty(nm), nm = sprintf('figure_%d', h.Number); end
    nm = regexprep(nm,'[^A-Za-z0-9]+','_');
    nm = regexprep(nm,'_+','_');
    nm = regexprep(nm,'^_|_$','');

    old = get(h,{'Units','PaperUnits','PaperPositionMode','PaperPosition', ...
                 'PaperSize','InvertHardcopy','Color','Renderer','RendererMode'});
    h.Units='centimeters';  pos = h.Position;
    h.PaperUnits='centimeters';  h.PaperPositionMode='manual';
    h.PaperPosition=[0 0 pos(3) pos(4)];  h.PaperSize=[pos(3) pos(4)];
    h.InvertHardcopy='off';  h.Color='w';
    h.Renderer='painters';   h.RendererMode='manual';
    set(findall(h,'Type','axes'),'Toolbar',[]);   % or it prints into the file
    drawnow expose                                 % let the layout settle

    switch lower(fmt)
        case 'eps', print(h,[nm '.eps'],'-depsc2','-vector');
        case 'svg', print(h,[nm '.svg'],'-dsvg','-vector');
        otherwise,  error('Sensitivity:fmt','Use eps or svg.');
    end
    set(h,{'Units','PaperUnits','PaperPositionMode','PaperPosition', ...
           'PaperSize','InvertHardcopy','Color','Renderer','RendererMode'},old);
    fprintf('    %s.%s   (%.1f x %.1f cm)\n', nm, lower(fmt), pos(3), pos(4));
end
fprintf('  Done -> %s\n', pwd);
end

function T = local_table_from_ckpt(ck)
%LOCAL_TABLE_FROM_CKPT  Results table from a checkpoint, without re-solving.
%
%   T = Sensitivity('load');            % OAT.mat in the current folder
%   T = Sensitivity('load','my.mat');   % a checkpoint saved under another name
%
% The checkpoint holds rows, done, runs and sig, and — once a sweep completes —
% the assembled table T. A stored T is returned as it stands, so the provenance
% attached at the end of the sweep survives. Older or partial checkpoints hold
% no T and are assembled here by the same path the sweep uses, so a table from
% either route is interchangeable with a swept one.
%
% A partial checkpoint loads: rows(done) keeps the completed entries and the
% table is short by the runs still outstanding. The count is reported so a
% partial file is not mistaken for a complete one.
%
% On the assembly path UserData carries the run list and the factor
% declaration, but not the baseline configuration struct: FC0 belongs to the
% solve and is not written to the checkpoint separately. It is present only
% when a stored T supplies it.
if ~exist(ck,'file')
    error('Sensitivity:noCheckpoint', ...
          'No checkpoint "%s" in %s.', ck, pwd);
end
S = load(ck);

% A CHECKPOINT WRITTEN BY A COMPLETED SWEEP ALREADY HOLDS T. Taking it as
% stored preserves the provenance struct, including the baseline configuration,
% which cannot be reconstructed from rows. Checkpoints written before T was
% stored, and partial ones, fall through to the assembly below.
if isfield(S,'T') && istable(S.T)
    T = S.T;
    fprintf('  Loaded %s: stored table, %d runs, %d factor(s).\n', ...
            ck, height(T), numel(unique(T.factor(~T.is_base))));
    return
end

need = {'rows','done'};
if ~all(isfield(S, need))
    error('Sensitivity:badCheckpoint', ...
          '"%s" holds %s; a sweep checkpoint holds rows, done, runs and sig.', ...
          ck, strjoin(fieldnames(S)', ', '));
end
kept = S.rows(logical(S.done));
if isempty(kept)
    error('Sensitivity:emptyCheckpoint','"%s" holds no completed runs.', ck);
end
proto = local_blank();
for i = 1:numel(kept), kept{i} = orderfields(kept{i}, proto); end
T = struct2table(vertcat(kept{:}));
T = movevars(T,{'run_id','factor','fac_lbl','block','level','is_base'},'Before',1);

ud = struct('source', ck, 'loaded_utc', ...
    datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS'));
if isfield(S,'runs'),  ud.runs    = S.runs;          end
ud.factors = local_factors();
T.Properties.UserData = ud;

fprintf('  Loaded %s: %d of %d runs, %d factor(s), baseline %s.\n', ...
        ck, height(T), numel(S.done), ...
        numel(unique(T.factor(~T.is_base))), ...
        ternary(any(T.is_base),'present','MISSING'));
if height(T) < numel(S.done)
    fprintf('  Partial checkpoint: %d run(s) outstanding.\n', ...
            numel(S.done) - height(T));
end
end


function S = local_style()
%LOCAL_STYLE  Type sizes and figure conventions for every figure in this file.
%
% The figure suite and the standalone figures reached through the dispatch
% (validity, pinch, tornado) are separate functions with separate workspaces,
% so the sizes are held here rather than in any one of them. Authored at final
% printed size; nothing is rescaled between MATLAB and the typeset page, so a
% figure MUST be placed in the manuscript at its native width. Type is fixed
% in points while the canvas is fixed in centimeters, so placing a canvas into
% a narrower column scales every label by the width ratio: a canvas half again
% wider than the column carries 9 pt type onto the page at 6 pt, below the
% publisher floor. Canvas widths therefore sit near the printed width.
%
% FIELDS
%   tick/lbl/ttl/leg/ann  [pt] type sizes at final printed size
%   ms                    marker size, matched to the type size
%   titles                false: in-plot titles are suppressed. The typesetter
%                         strips them and sets the caption instead, so a title
%                         on the axes is either duplicated or lost. Set true
%                         while working interactively.
%   font                  one family for the whole set
S = struct('tick',8, 'lbl',9, 'ttl',9, 'leg',8, 'ann',7.5, 'ms',3.5, ...
           'titles',false, 'font','Arial');
end


function varargout = local_title(varargin)
%LOCAL_TITLE  Title only when local_style().titles is on.
%
% Journal figures carry a caption, not a title. A title that states a finding
% duplicates the caption when both survive and loses the finding when the
% typesetter strips the artwork's own text. Every title call in this file goes
% through here, so the suite is switched from one place and the sentence stays
% in the code for whoever writes the caption.
% OUTPUT ONLY ON REQUEST. A function that always assigns its output echoes
% ans at the prompt on every call, so a suppressed title prints a line of
% console noise per figure.
ST = local_style();
if ST.titles, h = title(varargin{:}); else, h = gobjects(0); end
if nargout > 0, varargout{1} = h; end
end


function s = resp_axislabel(resp, y0)
%RESP_AXISLABEL  Axis label for a tornado on any response column.
%
% NO MATLAB IDENTIFIERS ON THE AXIS. A label built from the column name prints
% a variable the reader has never seen and an anchor value with no unit
% attached, so neither the quantity nor its scale is readable from the figure.
% The typeset symbol and its unit are declared here, once, beside the
% identifier they belong to.
%
% THE UNIT IS PRINTED ONCE, IN THE BRACKET. The bracket carries the unit of the
% plotted quantity and the parenthesis carries the baseline anchor as a bare
% number, read in that same unit. A unit symbol inside the parenthesis sets
% the same thing twice on one line.
switch resp
    case 'SEC_external'
        sym = '\Delta SEC_{external}';  unit = 'kWh m^{-3}';
        anch = sprintf('%.1f', y0);
    case 'R_still_pct'
        sym = '\Delta \itR\rm_{still}';  unit = 'percentage points';
        anch = sprintf('%.1f', y0);
    case 'prod_L_m2_day'
        sym = '\Delta areal productivity';  unit = 'L m^{-2} d^{-1}';
        anch = sprintf('%.2f', y0);
    otherwise
        sym = ['\Delta ' strrep(resp,'_','\_')];  unit = '';
        anch = sprintf('%.2f', y0);
end
if isempty(unit)
    s = sprintf('%s from baseline (%s)', sym, anch);
else
    s = sprintf('%s from baseline (%s)   [%s]', sym, anch, unit);
end
end


function s = fmtlev(x)
if isnan(x), s = '--'; else, s = sprintf('%g', x); end
end


function s = levcol(D, i, side)
%LEVCOL  One level column entry: the level and its deviation from the base.
%
% The deviation is a percentage on ratio-scaled factors and a kelvin
% difference on absolute temperatures, where a fraction of the base carries no
% physical meaning and is not comparable with the fractions on the other rows.
%
% A side holding no runs prints "--" with no deviation: nothing is computed
% there, as opposed to computed and undefined.
if strcmp(side,'lo')
    lv = D.level_lo(i);  pc = D.lo_pct(i);  dv = D.lo_dev(i);
else
    lv = D.level_hi(i);  pc = D.hi_pct(i);  dv = D.hi_dev(i);
end
if isnan(lv),      s = '--';
elseif D.absT(i),  s = sprintf('%s  (%+.0f K)', fmtlev(lv), dv);
else,              s = sprintf('%s  (%+.0f%%)',   fmtlev(lv), pc);
end
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end


function zldd_validity_envelope(T)
%ZLDD_VALIDITY_ENVELOPE  Which swept levels leave the model's valid range.
%
%   Sensitivity('validity', T)
%
%   NOT PART OF THE DEFAULT FIGURE SUITE. It is diagnostic rather than a
%   result, and the same information is in the results table as ok_all and
%   the individual ok_* columns, so it is drawn only when asked for.
%   Worth generating before submission even if it is not inserted: T_coil
%   fails a screen at every swept level, which is a statement the text
%   should make whether or not the figure carries it.
% THE ADMISSIBLE INTERVAL IS THE OBJECT, NOT THE INDIVIDUAL RUN. A marker per
% solved case puts a hundred dots on the panel, most of them passing and
% therefore carrying no information, while the quantity the reader wants --
% how far along its range each factor stays inside the model's fitted domain
% -- has to be inferred by scanning for the first colored dot. A shaded band
% over the failing portion of each row states it directly, and the swept
% levels are kept as small ticks so the sampling density remains visible.
%
% ONE COLOR PER SCREEN. Writing the failing screen name beside every marker
% produces a column of overlapping strings that run off the right edge at
% print size. Screens that fail nowhere are dropped from the legend, since a
% key for an absent symbol sends the reader hunting for it.
%
% BANDS ARE DRAWN PER CONTIGUOUS RUN OF FAILURES and per screen at reduced
% height, so a level that trips two screens shows both rather than one hiding
% the other, and a factor whose failures are split across both ends of its
% range is not drawn as a single band spanning the admissible middle.
%
% THE ONSET LEVEL IS PRINTED IN PHYSICAL UNITS at the edge of each band. The
% abscissa is normalized within each factor so that rows of different ranges
% are comparable, but a normalized position is not a number the reader can
% act on; the onset level is.
%
% ROWS ARE SORTED BY THE FRACTION OF LEVELS FAILING, most-constrained at the
% top. Factor declaration order carries no meaning here.
ST = local_style();
scr  = {'ok_recovery','ok_dewpoint','ok_salinity','ok_reynolds','ok_preheatcap'};
snm  = {'terminal recovery ceiling','glazing dew point', ...
        'salinity correlation','gap Reynolds number','preheat capacity'};
scol = [0.80 0.33 0.20; 0.35 0.30 0.62; 0.15 0.55 0.55; 0.85 0.60 0.15; 0.45 0.45 0.45];
scr  = scr(ismember(scr, T.Properties.VariableNames));
snm  = snm(1:numel(scr));  scol = scol(1:numel(scr),:);

facs = unique(T.factor(~T.is_base),'stable');
nfail = zeros(numel(facs),1); nlev = nfail;
for i = 1:numel(facs)
    S = T(T.factor==facs(i),:);
    bad = false(height(S),1);
    for s = 1:numel(scr), bad = bad | ~S.(scr{s}); end
    nfail(i) = sum(bad); nlev(i) = height(S);
end

% Screens with no failure anywhere are dropped entirely: they contribute no
% band, no legend key and no color.
active = false(1,numel(scr));
for s = 1:numel(scr), active(s) = any(~T.(scr{s})(~T.is_base)); end

[~,ord] = sortrows([nfail./nlev nfail]);    % least-constrained at the bottom
facs = facs(ord); nfail = nfail(ord); nlev = nlev(ord);

n  = numel(facs);
f6 = newfig('F6_validity_envelope', 17.0, 0.62*n + 3.6);
ax = axes; hold(ax,'on')

LM = 0.46;      % label column width, in normalized abscissa units
RM = 1.10;      % right edge, leaving room for the fail/levels column
for i = 1:n
    S  = sortrows(T(T.factor==facs(i),:),'level');
    xs = normlev(S.level);
    y  = i;                                    % most-constrained at the top,
                                               % rows filled from the bottom up

    if mod(i,2) == 1
        patch([-9 9 9 -9],[y-0.5 y-0.5 y+0.5 y+0.5],[0.957 0.957 0.957], ...
              'EdgeColor','none','HandleVisibility','off');
    end

    % the swept range, then the sampling as short ticks
    plot([0 1],[y y],'-','Color',[0.78 0.78 0.78],'LineWidth',1.1, ...
         'HandleVisibility','off')
    plot(xs, y*ones(size(xs)),'|','MarkerSize',4,'Color',[0.62 0.62 0.62], ...
         'LineWidth',0.7,'HandleVisibility','off')

    % half-spacing, used to extend a band to the midpoint of its neighbours
    if numel(xs) > 1, hstep = 0.5*median(diff(xs)); else, hstep = 0.04; end

    ns = nnz(active);  k = 0;
    for s = 1:numel(scr)
        if ~active(s), continue, end
        k = k + 1;
        m = ~S.(scr{s});
        if ~any(m), continue, end
        hh = ternary(ns > 1, 0.30/ns, 0.30);
        yc = ternary(ns > 1, y + 0.30 - (k-0.5)*2*hh, y);
        seg = local_runs(m(:).');
        for q = 1:size(seg,1)
            a = max(0, xs(seg(q,1)) - hstep);
            b = min(1, xs(seg(q,2)) + hstep);
            patch([a b b a],[yc-hh yc-hh yc+hh yc+hh], scol(s,:), ...
                  'EdgeColor','none','FaceAlpha',0.85,'HandleVisibility','off');
            % onset level, printed on the admissible side of the band edge
            % The boundary is printed as an inequality on the physical level,
            % not as a bare number: a bare "1" beside a band at the low end
            % does not say whether it is the last failing level or the first
            % admissible one. White backing keeps it off the sampling ticks.
            if seg(q,1) > 1
                text(a-0.010, yc, sprintf('>= %s', fmtlev(S.level(seg(q,1)))), ...
                     'HorizontalAlignment','right','VerticalAlignment','middle', ...
                     'FontSize',ST.ann,'Color',scol(s,:),'FontWeight','bold', ...
                     'BackgroundColor','w','Margin',0.5)
            elseif seg(q,2) < numel(xs)
                text(b+0.010, yc, sprintf('<= %s', fmtlev(S.level(seg(q,2)))), ...
                     'HorizontalAlignment','left','VerticalAlignment','middle', ...
                     'FontSize',ST.ann,'Color',scol(s,:),'FontWeight','bold', ...
                     'BackgroundColor','w','Margin',0.5)
            end
        end
    end

    text(-0.03, y, factor_axislabel(facs(i), S.fac_lbl(1)), ...
         'HorizontalAlignment','right','FontSize',ST.lbl, ...
         'VerticalAlignment','middle');
    if nfail(i) > 0, cc = [0.72 0.15 0.15]; wt = 'bold';
    else,            cc = [0.60 0.60 0.60]; wt = 'normal';
    end
    text(1.035, y, sprintf('%d/%d',nfail(i),nlev(i)),'HorizontalAlignment','left', ...
         'FontSize',ST.ann,'Color',cc,'FontWeight',wt);
end
text(1.035, n+0.78,'Failing / levels','FontSize',ST.ann, ...
     'Color',[0.35 0.35 0.35],'FontWeight','bold');

% legend proxies: the admissible range, then one key per screen that fails
h = plot(NaN,NaN,'-','Color',[0.78 0.78 0.78],'LineWidth',3, ...
         'DisplayName','Inside the fitted domain');
for s = 1:numel(scr)
    if ~active(s), continue, end
    h(end+1) = patch(NaN,NaN,scol(s,:),'EdgeColor','none', ...
                     'DisplayName',['Outside: ' snm{s}]);   %#ok<AGROW>
end
legend(h,'Orientation','horizontal','Location','northoutside','Box','off', ...
       'FontSize',ST.leg)

set(ax,'YTick',[],'XTick',[0 0.25 0.5 0.75 1], ...
       'XTickLabel',{'Lowest','','Mid','','Highest'},'FontSize',ST.tick)
xlabel('Swept level, normalized within each factor  [-]','FontSize',ST.lbl)
xlim([-LM RM]); ylim([0.4 n+1.0]); box off
set(ax,'YColor','none','Layer','top','FontName',ST.font)
tk = get(ax,'XTick'); set(ax,'XTick', tk(tk >= 0 & tk <= 1));

fprintf('  Validity envelope: %d of %d cases outside a fitted domain.\n', ...
        sum(nfail), sum(nlev));
end


function seg = local_runs(m)
%LOCAL_RUNS  Start and end indices of each contiguous true run in a logical row.
%
% A factor whose failures fall at both ends of its range must be drawn as two
% bands; spanning first to last failure would shade the admissible middle.
d = diff([false, logical(m), false]);
seg = [find(d == 1).', (find(d == -1) - 1).'];
end


function zldd_pinch_plot(T)
%ZLDD_PINCH_PLOT  SEC against cascade recovery, resolving the preheat pinch.
%
%   Sensitivity('pinch', T)
%
%   WHY THIS SUBSET. Only factors that alter the air-brine contacting duty at
%   FIXED WETTED AREA and FIXED recuperator effectiveness are fitted. The
%   restriction is not cosmetic:
%
%     Area-changing factors (A_plate, Np) retain 60-133 kWh/m3 of purchased
%     air heating at recoveries of 46-52 %, where the fitted factors need
%     none, because recovery and area set the vapor surplus independently.
%     They are drawn as open squares to show that they do NOT fall on the
%     curve -- omitting them would overstate the result.
%
%     eps_recup moves SEC by 73.5 kWh/m3 over its buildable range at recovery
%     constant to five decimal places, so it is orthogonal to this axis
%     entirely and is not plotted.
%
%   The claim the figure supports is therefore "within a fixed geometry, SEC
%   passes through a minimum near R_still = 50 %", not "SEC is a function of
%   R_still". The quadratic is a summary of that tendency, not a model: an rms
%   residual near 36 kWh/m3 is quoted on the axis so it cannot be read as one.
%
%   SCREENED RUNS ARE EXCLUDED FROM THE FIT and drawn hollow. A vertex set by
%   an inadmissible point locates the optimum where the model does not hold,
%   which is the failure this exclusion exists to prevent. It stays in force
%   whether or not any level trips a screen: with the ranges held inside the
%   admissible region nothing is excluded and the hollow markers report only
%   that status, and widening any range brings it back into play.

core = ["Qfan","T_air_in","Vfeed","T_coil","theta","hcomp"];
area = ["A_plate","Np"];

base = T(T.is_base,:);
V    = T(~T.is_base,:);
S    = V(ismember(V.factor,core) & V.ok_all, :);
if height(S) < 6
    error('Sensitivity:pinch','Too few admissible runs in the fitted subset.');
end

c   = polyfit(S.R_still_pct, S.SEC_external, 2);
vx  = -c(2)/(2*c(1));
yh  = polyval(c, S.R_still_pct);
R2  = 1 - sum((S.SEC_external-yh).^2)/sum((S.SEC_external-mean(S.SEC_external)).^2);
rms = sqrt(mean((S.SEC_external-yh).^2));

ST = local_style();
figure('Name','Fig4_SEC_vs_cascade_recovery','Color','w', ...
       'Units','centimeters','Position',[2 2 17 11]);
ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'off'); set(ax,'Layer','top')

xs = linspace(min(S.R_still_pct), max(S.R_still_pct), 200);
plot(xs, polyval(c,xs),'-','Color',[.55 .55 .55],'LineWidth',1.8, ...
     'DisplayName',sprintf('Quadratic fit  (\\itR\\rm^2 = %.2f, vertex %.1f %%)',R2,vx));
xline(vx,':','Color',[.55 .55 .55],'HandleVisibility','off');

co = lines(numel(core));
for i = 1:numel(core)
    g  = sortrows(V(V.factor==core(i),:),'R_still_pct');
    ok = g(g.ok_all,:);  bad = g(~g.ok_all,:);
    lb = ok.fac_lbl(1);
    plot(ok.R_still_pct, ok.SEC_external,'-o','Color',co(i,:),'MarkerSize',4.5, ...
         'MarkerFaceColor',co(i,:),'LineWidth',1.2,'DisplayName',char(lb));
    if ~isempty(bad)
        plot(bad.R_still_pct, bad.SEC_external,'o','Color',co(i,:), ...
             'MarkerSize',4.5,'LineWidth',1.1,'HandleVisibility','off');
    end
end
for j = 1:numel(area)
    g = sortrows(V(V.factor==area(j),:),'R_still_pct');
    h = plot(g.R_still_pct, g.SEC_external,'s','Color',[.45 .45 .45], ...
             'MarkerSize',5,'LineWidth',1.0);
    if j == 1, h.DisplayName = 'Area-changing factors (not fitted)';
    else,      h.HandleVisibility = 'off';
    end
end
plot(base.R_still_pct, base.SEC_external,'kp','MarkerSize',16, ...
     'MarkerFaceColor','k','DisplayName','Baseline');
plot(NaN,NaN,'o','Color',[.35 .35 .35],'MarkerSize',4.5,'LineWidth',1.1, ...
     'DisplayName','Outside the fitted salinity domain');

xlabel('Cascade recovery  \itR\rm_{still}  [%]','FontSize',ST.lbl);
ylabel('SEC_{external}  [kWh m^{-3}]','FontSize',ST.lbl);
% Three columns: the key holds nine entries, and each legend row costs
% height that comes out of the axes.
legend('Location','northoutside','NumColumns',3,'Box','off','FontSize',ST.leg);
set(ax,'FontName',ST.font,'FontSize',ST.tick);

fprintf(['  Pinch fit: n = %d admissible runs, vertex at R_still = %.1f %%,\n' ...
         '  R^2 = %.3f, rms residual %.1f kWh/m3. Baseline sits at %.1f %%.\n'], ...
        height(S), vx, R2, rms, base.R_still_pct);
end