function M = Morris_Global(varargin)
%ZLDD_MORRIS  Morris elementary-effects screening of the ZLDD cascade model.
%
%   M = Morris_Global('dryrun',true);        % plan and cost, solve nothing
%   M = Morris_Global('r',10);               % 10 trajectories, design block
%   M = Morris_Global('r',10,'plot',true);   % and draw the figures
%   M = Morris_Global('only',{'Qfan','Np','T_air_in','Vfeed'});
%   M = Morris_Global('resume',false);       % ignore an existing Global.mat
%   M = Morris_Global('file','Global_uncertainty.mat');   % a second study
%
%   SOLVE NOW, DRAW LATER. 'plot' is off by default; the figures are reachable
%   from a finished study without solving anything again:
%
%       M = Morris_Global('r',10);               % solve, no figures
%       load('Global.mat');                      % or come back to it later
%       h = Morris_Global('figs', M);            % draw
%       Morris_Global('figs', M, 'export','svg') % draw and save as vector
%
%   M reaches Global.mat only when the study completes, so a part-run file
%   carries its runs but not the statistics the figures plot.
%
%   load Global.mat                          % readable mid-run, not only at
%   sum(done), numel(runs)                   % the end: what is planned, what
%   rows{7}                                  % is finished, and every metric
%   runs(7)                                  % recorded for a finished solve
%
%   WHY THIS AND NOT ANOTHER OAT. The OAT study samples only the axes through
%   one baseline point, so two of its conclusions are conditional in a way it
%   cannot itself detect: the fan-flow optimum near 0.145 m^3/s holds only
%   with every other factor at baseline, and "this factor is unimportant" is
%   the weakest claim that design can support, because a factor inert along
%   the baseline axis may still be active through an interaction. Morris
%   resolves both at (k+1)*r solves rather than the ~1e4 a variance
%   decomposition needs. It is a SCREENING method: it ranks and it detects
%   interaction, but it does not apportion variance. Do not quote it as a
%   Sobol index.
%
%   WHAT THE STATISTICS MEAN
%     mu*   mean |elementary effect|. The importance ranking. Absolute values
%           are taken because this response is non-monotonic in fan flow,
%           inlet air temperature and plate area; plain mu cancels the two
%           limbs of those factors and reports near zero for a factor that
%           dominates the design. Report mu*.
%     mu    signed mean, kept as a diagnostic: mu << mu* is itself the
%           signature of a non-monotonic or interacting factor.
%     sigma spread of the effects across the sample. High sigma means the
%           influence depends on WHERE in the factor space it is measured --
%           interaction or nonlinearity. This is the statistic the OAT study
%           cannot produce, and the reason for running this at all.
%
%   INTERNAL CONTROL. eps_recup is linear to 0.07 kWh/m3 across the OAT
%   sweep, with cascade recovery constant to five decimals. It must therefore
%   return sigma ~ 0 here. If it does not, the implementation is wrong and no
%   other row in the table means anything. Check that row first.
%
%   SAMPLING. Standard Morris trajectory design (Morris 1991; Campolongo et
%   al. 2007) on the unit hypercube, so effects are dimensionless and
%   comparable across factors whose physical ranges differ by orders of
%   magnitude -- the defect that makes raw bar length a poor ranking in the
%   OAT tornado. With p levels and Delta = p/(2(p-1)) the design is symmetric and
%   each trajectory buys k effects for k+1 solves.
%
%   OPTIONS
%     'r'        trajectories, default 10.  Cost is (k+1)*r solves.
%     'p'        grid levels, default 4  (Delta = 2/3).
%     'only'     cellstr of factor names.  'block' selects a whole block.
%     'resp'     responses scored, default {'SEC_external','R_still_pct'}
%     'Nx'       grid, default 20, matching the OAT sweep.
%     'seed'     RNG seed, default 20260101. The design is random; a fixed
%                seed is what makes the study repeatable for a referee.
%     'screen'   true: effects whose endpoints fail a validity screen are
%                discarded rather than averaged in.
%     'timeout'  per-solve wall-clock cap in seconds, default 900 (15 min).
%     'file'     data file, default 'Global.mat'. '' disables all writing.
%     'resume'   true: reuse an existing file whose design matches.
%     'export'   'svg' | 'eps': write the Morris figures as vector after
%                drawing them. Only the Morris figures, not everything open.
%     'range'    struct narrowing a factor's range FOR THIS STUDY ONLY, e.g.
%                struct('T_coil',[282 295],'eps_recup',[0.70 0.85]).
%
%                WHY THIS EXISTS. The declared level lists are the right
%                ranges for a sweep, where an inadmissible level costs one run
%                and is reported as the validity envelope. In a trajectory
%                design the cost is quite different: T_coil >= 298 K breaches
%                the recovery ceiling whatever the other factors do, so
%                every point drawn there dies, and a dead point destroys BOTH
%                elementary effects adjoining it. Sampling 282-308 K therefore
%                loses roughly three quarters of the effects and leaves two or
%                three per factor -- too few to separate anything.
%
%                It is a RESTRICTION, never an extension: a range wider than
%                the declaration is rejected, because that would sample the
%                model where the sweep never checked it and would put a second
%                source of truth beside local_factors(). Restrictions are
%                echoed in the header and stored in the design signature, so a
%                checkpoint from a differently restricted study will not
%                resume onto this one.
%     'dryrun'   plan only.     'plot'  draw the figures.
%
%   WHAT IS SAVED, AND WHEN. One file, 'Global.mat', laid out like the OAT
%   sweep's: RUNS is the complete plan, DONE marks what has finished, ROWS
%   holds one metric struct per finished solve. The file is WRITTEN BEFORE THE
%   FIRST SOLVE and rewritten after every solve, so at any instant the file on
%   disk describes the whole study and carries every run completed so far.
%   Ctrl-C, a crash or a power cut costs the solve in flight, nothing else.
%
%   The plan is complete up front because every applied value is a
%   deterministic function of the seed, the factor set and the ranges --
%   nothing about run 87 waits on run 86. So RUNS is filled in at design time,
%   which is what lets a half-finished file be read as a study rather than as
%   a heap of numbers: the rows not yet done are already named.
%
%   The same ~50 quantities the sweep records are captured for EVERY solve,
%   not only the responses being scored, through Sensitivity('extract'). A
%   study costing 3-4 h should not have to be re-run because a reviewer asks
%   about brine salinity or the dew-point margin; ROWS carries them all, and
%   elementary effects for any column can be recomputed from it without
%   solving anything.
%
%   Requires Sensitivity.m on the path: ranges and setters are read from its
%   declaration via Sensitivity('factors'), so the two studies cannot drift.

%% --------------------------------------------------------------- commands
% DRAW-LATER ENTRY POINT, taken before the option parser because 'figs' is a
% command and not a name-value pair. A study solved without 'plot' is not a
% study that has to be solved again to be looked at.
if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1})) && ...
        strcmpi(char(varargin{1}), 'figs')
    if nargin < 2 || ~isstruct(varargin{2})
        error('Morris_Global:figsNeedsM', ...
             ['Morris_Global(''figs'', M) needs the study struct. Either keep ' ...
              'the M returned by the solving call, or load it: ' ...
              'load(''Global.mat''), which carries M once the study completes.']);
    end
    M = Morris_Global_figs(varargin{2});    % figure handles, not the study
    if nargin >= 4 && strcmpi(char(varargin{3}), 'export')
        Sensitivity('export', varargin{4}, M);
    end
    return
end

%% ---------------------------------------------------------------- options
opt = struct('r',10,'p',4,'only',{{}},'block','design', ...
             'resp',{{'SEC_external','R_still_pct'}},'Nx',20, ...
             'seed',20260101,'screen',true,'timeout',900, ...
             'file','Global.mat','resume',true, ...
             'dryrun',false,'plot',false,'export','','range',struct());
% UNKNOWN OPTIONS ARE AN ERROR, NOT A NEW FIELD. A bare
% opt.(varargin{i}) = varargin{i+1} accepts anything: a misspelling, or an
% option belonging to a newer version of this file, is silently stored and
% then silently ignored, so the study runs to completion with the setting the
% caller believed they had changed still at its default. That failure is
% invisible and expensive at 3-4 h per run.
known = fieldnames(opt);
if mod(numel(varargin),2) ~= 0
    error('Morris_Global:pairs','Options must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    nm = varargin{i};
    if ~(ischar(nm) || isstring(nm)) || ~ismember(char(nm), known)
        error('Morris_Global:unknownOption', ...
              ['Unknown option "%s". Known options: %s.\n' ...
               'If this option exists in a newer version of the file, the ' ...
               'copy on the path is out of date.'], ...
              char(string(nm)), strjoin(known.', ', '));
    end
    opt.(char(nm)) = varargin{i+1};
end

if exist('Sensitivity','file') ~= 2
    error('Morris_Global:noSweep','Sensitivity.m must be on the path.');
end
F = Sensitivity('factors');

if ~isempty(opt.only)
    missing = setdiff(opt.only, {F.name});
    if ~isempty(missing)
        error('Morris_Global:unknownFactor','Not declared: %s', strjoin(missing,', '));
    end
    keep = ismember({F.name}, opt.only);
elseif ~strcmpi(opt.block,'all')
    keep = strcmpi({F.block}, opt.block);
else
    keep = true(1,numel(F));
end
F = F(keep);  k = numel(F);
if k < 2, error('Morris_Global:tooFew','Need at least two factors.'); end

% Ranges are the endpoints of the OAT level list. Taking them from the same
% declaration is what keeps the two studies comparable; a range typed here
% would be a second source of truth and would drift.
lo = arrayfun(@(f) min(f.levels), F);
hi = arrayfun(@(f) max(f.levels), F);

narrowed = strings(0);
rf = fieldnames(opt.range);
for a = 1:numel(rf)
    i = find(strcmp({F.name}, rf{a}), 1);
    if isempty(i)
        if ismember(rf{a}, opt.only) || isempty(opt.only)
            warning('Morris_Global:rangeUnused', ...
                    'Range given for %s, which is not in this factor set.', rf{a});
        end
        continue
    end
    rg = opt.range.(rf{a});
    if numel(rg) ~= 2 || rg(2) <= rg(1)
        error('Morris_Global:badRange','range.%s must be [lo hi].', rf{a});
    end
    if rg(1) < lo(i)-eps(100) || rg(2) > hi(i)+eps(100)
        error('Morris_Global:rangeWider', ...
             ['range.%s = [%g %g] falls outside the declared [%g %g]. ' ...
              'A study may restrict a range but not extend one -- the model ' ...
              'has not been checked there.'], rf{a}, rg(1), rg(2), lo(i), hi(i));
    end
    narrowed(end+1,1) = sprintf('%s [%g %g]', rf{a}, rg(1), rg(2)); %#ok<AGROW>
    lo(i) = rg(1);  hi(i) = rg(2);
end
if any(hi <= lo), error('Morris_Global:degenerate','A factor has zero range.'); end

% Integer factors are detected from the declaration rather than named here, so
% a new one is picked up without editing this file.
isint = arrayfun(@(f) all(abs(f.levels-round(f.levels)) < eps(100)) && ...
                      (max(f.levels)-min(f.levels)) >= 2 && strcmp(f.unit,'-'), F);

%% ------------------------------------------------------- trajectory design
rng(opt.seed,'twister');
p     = opt.p;
Delta = p/(2*(p-1));
grid_ = linspace(0, 1-Delta, p/2);

X = cell(opt.r,1);
for t = 1:opt.r
    B  = tril(ones(k+1,k),-1);
    Ds = diag(2*randi([0 1],k,1)-1);
    Ps = eye(k); Ps = Ps(randperm(k),:);
    xs = grid_(randi(numel(grid_),1,k));
    X{t} = min(max((ones(k+1,1)*xs + ...
           (Delta/2)*((2*B-ones(k+1,k))*Ds+ones(k+1,k)))*Ps, 0), 1);
end
N = opt.r*(k+1);

fprintf('\n=================================================================\n');
fprintf('  ZLDD MORRIS SCREENING\n');
fprintf('  Factors : %d   Trajectories : %d   Levels : %d   Delta = %.3f\n', ...
        k, opt.r, p, Delta);
fprintf('  Solves  : %d   Grid : Nx = %d   Seed : %d\n', N, opt.Nx, opt.seed);
fprintf('  Timeout : %g s per solve\n', opt.timeout);
fprintf('-----------------------------------------------------------------\n');
for i = 1:k
    dec = [min(F(i).levels) max(F(i).levels)];
    tag = '';
    if abs(lo(i)-dec(1)) > eps(100) || abs(hi(i)-dec(2)) > eps(100)
        tag = sprintf('   restricted from [%g %g]', dec(1), dec(2));
    end
    fprintf('   %-16s %10.4g  ->%10.4g   [%s]%s%s\n', F(i).name, lo(i), hi(i), ...
            F(i).unit, tern(isint(i),'  integer',''), tag);
end
if opt.dryrun
    fprintf('-----------------------------------------------------------------\n');
    fprintf('  DRY RUN -- nothing solved.\n');
    M = struct('factors',{F},'X',{X},'Delta',Delta,'lo',lo,'hi',hi,'isint',isint);
    return
end

%% --------------------------------------------------------------- timeout
% TRUE PREEMPTION IS NOT POSSIBLE FROM THE CALLER in plain MATLAB: once the
% ODE integrator is inside a step, no timer in this workspace can interrupt
% it. parfeval solves it properly -- the solve runs on a worker and CANCEL
% kills it -- so the toolbox is used when it is present and the cap is real.
%
% Without the toolbox the study still runs, but a pathological point holds the
% loop for as long as it takes. The checkpoint below is what limits the damage
% then: the run is resumable, so killing MATLAB by hand costs one solve rather
% than the study. The clean alternative is a wall-clock check inside the
% model's ODE OutputFcn, which belongs in ZLDD_complete_modeling, not here.
usePar = opt.timeout > 0 && ~isempty(ver('parallel')) && ...
         license('test','Distrib_Computing_Toolbox');
poolFailed = false;
if usePar
    % A POOL THAT WILL NOT START COSTS THE TIMEOUT, NOT THE STUDY. The licence
    % check above says the toolbox is installed; it cannot say the pool will
    % come up. A stale job, a blocked loopback port or an interrupted shutdown
    % all fail here, and an uncaught failure ends the whole study before the
    % first solve -- for want of a cap that only ever truncates stragglers.
    try
        if isempty(gcp('nocreate')), parpool('local',1); end
        fprintf('  Per-solve timeout active via parfeval (1 worker).\n');
    catch ME
        usePar = false; poolFailed = true;
        warning('Morris_Global:poolFailed', ...
            ['Parallel pool would not start (%s) -- continuing without the ' ...
             'per-solve timeout. The run is saved after every solve, so a ' ...
             'stalled point can be interrupted by hand and resumed. To clear ' ...
             'a stale pool: delete(gcp(''nocreate'')), then restart MATLAB.'], ...
            ME.message);
    end
end
if ~usePar && ~poolFailed && opt.timeout > 0
    warning('Morris_Global:noTimeout', ...
        ['Parallel Computing Toolbox not available -- the per-solve timeout ' ...
         'cannot be enforced. The run is checkpointed, so interrupting a ' ...
         'stalled solve by hand costs one point, not the study.']);
end

%% ----------------------------------------------------------- the run plan
% Every applied value follows from the design alone, so the entire plan is
% built here, before anything is solved. Integer factors are rounded at the
% same point: the coordinate actually applied is a property of the plan, not
% an outcome of the solve.
nr   = numel(opt.resp);
ID   = @(t,j) (t-1)*(k+1)+j;     % trajectory/point -> linear run id

V    = zeros(N,k);               % applied physical values
Xapp = zeros(N,k);               % applied unit coordinates, after rounding
for t = 1:opt.r
    for j = 1:k+1
        v = lo + X{t}(j,:).*(hi-lo);
        % THE APPLIED COORDINATE IS RECORDED, NOT THE REQUESTED ONE. Rounding
        % an integer factor moves the point off the Morris grid: Np steps
        % 5 -> 11.67, which the model receives as 12, so the perturbation
        % actually applied is 0.700 of the range where the design asked for
        % 0.667. Dividing by the nominal Delta would inflate that factor's
        % elementary effect by 5 % -- a bias in the ranking, not noise in it.
        v(isint)      = round(v(isint));
        V(ID(t,j),:)  = v;
        Xapp(ID(t,j),:) = (v-lo)./(hi-lo);
    end
end

tmpl = struct('id',0,'traj',0,'point',0,'factor','','label','','unit','', ...
              'block','','value',NaN,'step',NaN,'is_start',false, ...
              'x',zeros(1,k),'v',zeros(1,k));
runs = repmat(tmpl, N, 1);
for t = 1:opt.r
    for j = 1:k+1
        id = ID(t,j);
        runs(id).id    = id;
        runs(id).traj  = t;
        runs(id).point = j;
        runs(id).x     = Xapp(id,:);
        runs(id).v     = V(id,:);
        if j == 1
            % The anchor of the trajectory: no factor has moved yet.
            runs(id).is_start = true;
            runs(id).factor   = 'START';
            runs(id).label    = 'Trajectory start';
            runs(id).unit     = '-';
            runs(id).block    = 'start';
        else
            i = find(abs(X{t}(j,:) - X{t}(j-1,:)) > 1e-9, 1);
            runs(id).factor = F(i).name;
            runs(id).label  = F(i).label;
            runs(id).unit   = F(i).unit;
            runs(id).block  = F(i).block;
            runs(id).value  = V(id,i);
            runs(id).step   = Xapp(id,i) - Xapp(ID(t,j-1),i);
        end
    end
end

rows = cell(N,1);                % full metric struct per finished solve
done = false(N,1);
Y    = nan(N, nr);
OK   = false(N,1);
WALL = nan(N,1);
STAT = strings(N,1);
ERRM = strings(N,1);             % why a run failed, when it did

% Factor metadata WITHOUT the setters. Function handles saved into a .mat
% capture their workspace and reload as broken references if Sensitivity.m
% has changed; the names are what a reader of the file needs.
fac = struct('name',{F.name},'label',{F.label},'unit',{F.unit}, ...
             'block',{F.block},'lo',num2cell(lo),'hi',num2cell(hi), ...
             'isint',num2cell(isint));

% THE DESIGN, NOT THE COUNT, IS WHAT A RESUME MUST MATCH. Two studies with
% the same number of solves but a different seed, factor set, range or grid
% are different studies, and resuming one onto the other would silently mix
% them into a table that looks valid and is not.
sig = struct('names',{{F.name}},'lo',lo,'hi',hi,'r',opt.r,'p',opt.p, ...
             'seed',opt.seed,'Nx',opt.Nx,'resp',{opt.resp});

if opt.resume && ~isempty(opt.file) && exist(opt.file,'file')
    S = load(opt.file);
    % The layout is checked as well as the design. A file carrying a matching
    % sig but a different set of variables would pass the sig test and then
    % throw on the field reads below, so both are tested before anything is
    % taken from it.
    if isfield(S,'sig') && isequaln(S.sig, sig) && ...
       all(isfield(S,{'done','rows','Y','OK','WALL','STAT','ERRM'}))
        done=S.done; Y=S.Y; OK=S.OK; WALL=S.WALL; STAT=S.STAT; rows=S.rows;
        ERRM=S.ERRM;
        fprintf('  Resuming %s: %d of %d solves already complete.\n', ...
                opt.file, sum(done), N);
    else
        fprintf('  %s holds a different design or layout -- starting fresh.\n', ...
                opt.file);
    end
end

% THE FILE IS CREATED BEFORE THE FIRST SOLVE. A study that dies at solve 3 of
% 156 still leaves a file that says what the study was and that three runs are
% in it.
if ~isempty(opt.file)
    save(opt.file,'runs','rows','done','sig','fac','X','V','Xapp', ...
                  'Y','OK','WALL','STAT','ERRM','Delta','lo','hi','isint','-v7');
    fprintf('  Data file: %s  (%d runs planned)\n', opt.file, N);
end

%% ------------------------------------------------------------------ solve
FC0 = Baseline_Input_Variable();
FC0.Nx = opt.Nx;
t0 = tic;
warned = false;

for t = 1:opt.r
    for j = 1:k+1
        id = ID(t,j);
        if done(id), continue, end

        FC = FC0; FC.verbose = false; FC.quiet = true;
        v  = V(id,:);
        for i = 1:k, FC = F(i).setter(FC, v(i)); end

        tic_run = tic;
        [res, status, errmsg] = solve_one(FC, usePar, opt.timeout);
        WALL(id) = toc(tic_run);
        STAT(id) = status;
        ERRM(id) = errmsg;

        if status == "ok"
            % SCORED FROM THE EXTRACT, NOT THE RAW STRUCT. Sensitivity('extract')
            % is where the reported metrics are named and where the derived ones
            % are computed, so a response that exists only after extraction --
            % SEC_external among them -- reads as absent from res, and gv()
            % answers a missing field with NaN rather than an error. Scoring the
            % raw struct therefore screens out every solve in the study and
            % reports it as a validity failure. The raw struct is kept only as
            % the fallback for a failed extraction.
            src = res;
            try
                R = Sensitivity('extract', res, FC);
                R.run_id = id;  R.traj = t;  R.point = j;
                R.factor = string(runs(id).factor);
                R.fac_lbl = string(runs(id).label);
                R.block  = string(runs(id).block);
                R.is_start = runs(id).is_start;
                R.wall_s = WALL(id);
                for i = 1:k, R.(['x_' F(i).name]) = v(i); end
                rows{id} = R;
                src = R;
            catch
                % Extraction is a convenience, not the study. If it fails the
                % elementary effects are still computed from Y.
            end
            m = local_score(src, opt.resp);
            Y(id,:) = m.y(:).';
            OK(id)  = m.ok;

            % A RESPONSE THAT IS NEVER FOUND IS A TYPO OR A RENAME, NOT A
            % VALIDITY FAILURE, AND THE TWO LOOK IDENTICAL IN THE LOG. gv()
            % answers a missing field with NaN, so an unscoreable response
            % screens out all N solves and the study reports itself as a
            % table of zero effects hours later. Said once, on the first
            % solve that succeeds, while the run can still be stopped.
            if ~warned && isstruct(src)
                gone = opt.resp(~isfield(src, opt.resp));
                if ~isempty(gone)
                    warning('Morris_Global:respMissing', ...
                        ['Response(s) %s are not fields of the extract. They ' ...
                         'will score NaN on every solve and every elementary ' ...
                         'effect will be discarded. Stop the run and check ' ...
                         'the names against fieldnames(Sensitivity(''extract'',...)).'], ...
                        strjoin(gone, ', '));
                end
                warned = true;
            end
        end
        done(id) = true;

        if status == "ok"
            tag = tern(OK(id),'[ok]','[screened]');
        else
            tag = char("[" + status + "]");
        end
        fprintf('   %3d/%3d  traj %2d pt %2d  %-14s %s = %10.3f  %-11s %5.0f s\n', ...
                id, N, t, j, runs(id).factor, opt.resp{1}, Y(id,1), tag, WALL(id));

        % THE REASON IS PRINTED WHERE THE FAILURE IS, not left for the summary.
        % Three "errors" in a row mean something different from three scattered
        % through the study, and the difference is only visible in the message.
        if strlength(ERRM(id)) > 0
            first = extractBefore(ERRM(id) + newline, newline);
            fprintf('            %s\n', first);
        end

        % REWRITTEN AFTER EVERY SOLVE, whole file, same variable list as the
        % write before the loop. One rule governs the file's contents at every
        % instant of the run, so a copy taken at any moment reads the same way.
        if ~isempty(opt.file)
            save(opt.file,'runs','rows','done','sig','fac','X','V','Xapp', ...
                          'Y','OK','WALL','STAT','ERRM','Delta','lo','hi','isint','-v7');
        end
    end
end

nTO = sum(STAT == "timeout");
nER = sum(STAT == "error");

%% ------------------------------------------------- elementary effects
% Consecutive points on a trajectory differ in exactly one coordinate by
% +/-Delta. Dividing by that SIGNED step preserves the sign, so mu keeps its
% meaning as a diagnostic alongside mu*.
EE = nan(opt.r, k, nr);
for t = 1:opt.r
    for j = 1:k
        a = ID(t,j);  b = ID(t,j+1);
        if ~(done(a) && done(b)), continue, end
        d = Xapp(b,:) - Xapp(a,:);
        i = find(abs(d) > 1e-9, 1);
        if isempty(i), continue, end
        if opt.screen && ~(OK(a) && OK(b)), continue, end
        for q = 1:nr
            EE(t,i,q) = (Y(b,q) - Y(a,q)) / d(i);
        end
    end
end

% FULL PER-SOLVE TABLE, assembled the way the sweep assembles its own: one
% fixed schema per row, ordered against a prototype, so a plain vertcat is
% safe and stays safe on releases where vertcat is field-order sensitive.
% Solves that timed out or errored contribute no row.
RT  = table();
got = ~cellfun(@isempty, rows);
if any(got)
    parts = rows(got);
    proto = parts{1};
    keep  = true(numel(parts),1);
    for i = 2:numel(parts)
        if isequal(sort(fieldnames(parts{i})), sort(fieldnames(proto)))
            parts{i} = orderfields(parts{i}, proto);
        else
            keep(i) = false;
        end
    end
    RT = struct2table(vertcat(parts{keep}));
end

M = struct('factors',{F},'Delta',Delta,'lo',lo,'hi',hi,'isint',isint, ...
           'X',{X},'Xapp',Xapp,'V',V,'Y',Y,'OK',OK,'done',done,'WALL',WALL, ...
           'STAT',STAT,'ERRM',ERRM,'EE',EE,'plan',runs,'runs',RT,'resp',{opt.resp}, ...
           'opt',opt,'sig',sig);

for q = 1:nr
    e   = EE(:,:,q);
    nEE = sum(~isnan(e),1).';
    mus = mean(abs(e),1,'omitnan').';
    mu  = mean(e,1,'omitnan').';
    sg  = std(e,0,1,'omitnan').';
    Tq  = table(string({F.name}).', string({F.label}).', string({F.unit}).', ...
                mus, mu, sg, sg./sqrt(max(nEE,1)), nEE, ...
          'VariableNames',{'factor','label','unit','mu_star','mu','sigma','sem','n_EE'});
    Tq  = sortrows(Tq,'mu_star','descend');
    Tq.rank = (1:height(Tq)).';
    M.(matlab.lang.makeValidName(opt.resp{q})) = Tq;

    fprintf('\n-----------------------------------------------------------------\n');
    fprintf('  %s -- Morris screening (%d trajectories)\n', opt.resp{q}, opt.r);
    disp(Tq);
    lost = opt.r*k - sum(nEE);
    if lost > 0
        fprintf('  %d of %d effects discarded on validity screens or failed solves.\n', ...
                lost, opt.r*k);
    end
end

fprintf('\n-----------------------------------------------------------------\n');
fprintf('  %d solves in %.1f min.  median %.0f s, slowest %.0f s.\n', ...
        sum(done), toc(t0)/60, median(WALL,'omitnan'), max(WALL));
if nTO > 0
    fprintf('  %d solve(s) hit the %g s cap and were abandoned.\n', nTO, opt.timeout);
end
if nER > 0, fprintf('  %d solve(s) errored.\n', nER); end

%% --------------------------------------------------------------- persist
% M IS APPENDED TO THE SAME FILE. The per-solve record and the statistics
% computed from it belong together: statistics separable from the runs behind
% them are statistics that can be quoted without anyone being able to check
% where the numbers came from.
if ~isempty(opt.file)
    save(opt.file,'M','-append');
    fprintf('  Written: %s  (%d of %d runs complete, M appended)\n', ...
            opt.file, sum(done), N);
end
fprintf('=================================================================\n');

if opt.plot
    h = Morris_Global_figs(M);
    if ~isempty(opt.export)
        Sensitivity('export', opt.export, h);   % these figures only
    end
elseif ~isempty(opt.export)
    warning('Morris_Global:exportNoPlot', ...
            '''export'' has no effect without ''plot'', true.');
end
end


%% ===================================================================== util
function [res, status, errmsg] = solve_one(FC, usePar, timeout)
% One solve, capped in wall time when the toolbox allows it.
%
% THE MESSAGE IS RETURNED, NOT DISCARDED. A status of "error" alone cannot
% distinguish a model that failed on a hard corner of the space from a pool
% that has stopped dispatching, and the two call for opposite responses:
% the first is a point to drop, the second is a run to stop.
res = []; status = "ok"; errmsg = "";
if usePar
    % A CANCELLED FUTURE CAN TAKE THE POOL WITH IT. Cancelling a worker
    % mid-ODE sometimes leaves the pool unusable, and every later parfeval
    % then fails immediately -- a whole study of "errors" that are nothing to
    % do with the model. The pool is checked before dispatch and restarted
    % once; if it will not come back the solve runs in-process, losing the
    % timeout rather than the study.
    p = gcp('nocreate');
    if isempty(p)
        try
            p = parpool('local',1);
        catch
            p = [];
        end
    end
    if isempty(p)
        try
            res = ZLDD_complete_modeling(FC);
            errmsg = "pool lost; solved in-process, uncapped";
        catch ME
            status = "error"; errmsg = string(ME.message);
        end
        return
    end
    fut = parfeval(p, @ZLDD_complete_modeling, 1, FC);
    if ~wait(fut, 'finished', timeout)
        cancel(fut); status = "timeout"; return
    end
    if ~isempty(fut.Error)
        status = "error"; errmsg = string(fut.Error.message); return
    end
    res = fetchOutputs(fut);
else
    try
        res = ZLDD_complete_modeling(FC);
    catch ME
        status = "error"; errmsg = string(ME.message);
    end
end
end

function m = local_score(results, resp)
% Scored responses and the validity flag from one solve. A missing field
% scores NaN rather than throwing: one odd trajectory point must not destroy
% a study whose solves have already been paid for.
m.y = nan(1,numel(resp));
for q = 1:numel(resp), m.y(q) = gv(results, resp{q}); end

% Convention matches ok_all in the OAT study: the recovery ceiling and the
% glazing dew-point margin are model failure and disqualify a point; the
% salinity correlation limit marks a closure used outside its fitted range,
% which is a caveat on the number rather than grounds to discard it.
ok  = isfinite(m.y(1));
dew = gv(results,'dew_margin_K');  if ~isnan(dew), ok = ok && dew > 0; end
rec = gv(results,'R_plant_pct');   if ~isnan(rec), ok = ok && rec <= 100.0001; end
m.ok = ok;
end

function v = gv(S,name)
if isstruct(S) && isfield(S,name)
    x = S.(name);
    if isnumeric(x) && ~isempty(x), v = double(x(end)); return, end
end
v = NaN;
end

function v = tern(c,a,b)
if c, v = a; else, v = b; end
end

%% ==================================================================== figs
function h = Morris_Global_figs(M)
%MORRIS_GLOBAL_FIGS  The (mu*, sigma) plane and the mu* ranking, per response.
%
%   h = Morris_Global_figs(M)
%
%   Draws from a finished study, whether M was just returned or loaded:
%
%       M = Morris_Global('r', 10);          % solve now, no figures
%       load('Global.mat');                  % or come back to it later
%       h = Morris_Global('figs', M);        % draw
%       Sensitivity('export', 'svg', h);     % save
%
%   M is appended to Global.mat only when the study completes, so a part-run
%   file carries its runs but not the statistics these figures plot.
%
%   ONE FIGURE SET WITH THE SWEEP. Canvas geometry, type sizes and the title
%   convention come from morris_style(), which reads the sweep's declaration
%   when Sensitivity.m exposes it and mirrors it otherwise, so these panels
%   land on the page beside the sweep's figures as one set. Sizes are authored
%   at final printed width in centimeters: place a panel at its native width
%   or the type reaches the page below the size it was set in.
%
%   NO MODEL IDENTIFIERS ON THE ARTWORK. mu* carries the unit of its response,
%   which differs between them -- kWh m^-3 for energy, percentage points for
%   recovery -- so the axis label is declared per response below rather than
%   built from the column name.
%
%   Returns the figure handles so the caller exports exactly these and not
%   whatever else happens to be open from another study.
resp = M.resp;  h = gobjects(0);

ST = morris_style();
% Canvas widths, in centimeters at final printed size. The plane labels its
% markers in symbols, so it needs no extra width for text and takes the
% sweep's single-axes frame; the ranking sets full factor names in a left
% margin and takes the tornado's width.
FIG_PL = 12.0;                 % [cm] (mu*, sigma) plane
FIG_PH = 9.5;                  % [cm] its height
FIG_RK = 17.0;                 % [cm] ranking panels, matching the tornado

C_pt  = [0.20 0.45 0.70];      % factor marker and bar
C_ann = [0.45 0.45 0.45];      % secondary annotation
C_ref = [0.60 0.60 0.60];      % reference lines

old_defaults = get(0, {'DefaultAxesFontName','DefaultAxesFontSize', ...
                       'DefaultTextFontName','DefaultLineLineWidth', ...
                       'DefaultAxesLabelFontSizeMultiplier', ...
                       'DefaultAxesTitleFontSizeMultiplier'});
set(0,'DefaultAxesFontName',ST.font,'DefaultAxesFontSize',ST.tick, ...
      'DefaultTextFontName',ST.font,'DefaultLineLineWidth',1.0, ...
      'DefaultAxesLabelFontSizeMultiplier',1.0, ...
      'DefaultAxesTitleFontSizeMultiplier',1.0);
cleanupDefaults = onCleanup(@() set(0, ...
    {'DefaultAxesFontName','DefaultAxesFontSize','DefaultTextFontName', ...
     'DefaultLineLineWidth','DefaultAxesLabelFontSizeMultiplier', ...
     'DefaultAxesTitleFontSizeMultiplier'}, old_defaults)); %#ok<NASGU>

for q = 1:numel(resp)
    Tq = M.(matlab.lang.makeValidName(resp{q}));
    xm = max(Tq.mu_star);
    if ~isfinite(xm) || xm <= 0, continue, end
    [sym, unit, tag] = morris_resp_label(resp{q});

    %% ---- M1  (mu*, sigma) plane ---------------------------------------
    % The diagonal is what makes the plane readable: below it a factor's
    % effect is roughly the same wherever it is measured, above it the effect
    % depends on where the other factors sit. A factor near the origin is
    % inert over the sampled space -- the claim the OAT design cannot support,
    % and the reason this panel exists.
    h(end+1) = newfig_m(sprintf('M1_morris_plane_%s', tag), FIG_PL, FIG_PH); %#ok<AGROW>
    ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'off'); set(ax,'Layer','top')

    xr = xm*1.32;
    plot([0 xr],[0 xr],'-','Color',C_ref,'LineWidth',1.0, ...
         'DisplayName','\sigma = \mu^*');
    plot([0 xr],0.5*[0 xr],':','Color',C_ref,'LineWidth',1.0, ...
         'DisplayName','\sigma = 0.5 \mu^*');
    % MARKER AREA FROM THE DECLARED SIZE. scatter takes an AREA in points
    % squared where every other MATLAB call takes a diameter, so a literal
    % here silently sets a marker unrelated to the one the sweep draws; ms is
    % a diameter in points and is squared to match.
    scatter(Tq.mu_star, Tq.sigma, ST.ms^2, C_pt,'filled', ...
            'MarkerEdgeColor','w','LineWidth',0.4,'HandleVisibility','off');

    % LABELS ARE PLACED AWAY FROM THE CROWD, not always to the right. The
    % weak factors pile into the origin, where a uniform offset overprints
    % them; each label is pushed to whichever side of its marker is emptier.
    yr = max([Tq.sigma; eps])*1.18;

    % LABEL PLACEMENT BY EXTENT, NOT BY OFFSET. A left-or-right rule separates
    % neighbours but not coincidences, and the weak factors are coincident:
    % on the recovery response the recuperator sits at the origin and feed
    % supply temperature a hundredth of the range away, so any purely
    % horizontal offset prints one name through the other. Each label is
    % treated as a rectangle whose size follows from the type size and the
    % canvas, and candidate positions are tried in order until one is clear.
    % Placing the largest mu* first gives the factors the reader looks for
    % their natural position and pushes the crowd at the origin outward.
    ax_w = 0.78*FIG_PL;  ax_h = 0.75*FIG_PH;        % [cm] drawn axes box
    ch_w = 0.50*ST.ann/28.35;                       % [cm] mean glyph width
    lh   = (ST.ann/28.35)*1.45/ax_h*yr;             % [axis] line height
    placed = zeros(0,4);                            % [x0 x1 y0 y1] rectangles

    % THE LEGEND IS RESERVED BEFORE ANY LABEL IS SET. It is drawn last but
    % occupies the upper left from the moment the axes exist, and a label
    % placed there is overprinted by it rather than colliding with anything
    % the collision test can see.
    placed(end+1,:) = [0 0.34*xr 0.80*yr yr];

    [~, ord] = sort(Tq.mu_star, 'descend');
    for ii = 1:numel(ord)
        i  = ord(ii);
        s  = morris_symbol(char(Tq.factor(i)), char(Tq.label(i)));
        % Width is measured on the typeset result, not the markup: the braces
        % and backslash of a TeX subscript are instructions, not glyphs, and
        % counting them would reserve roughly twice the space a symbol needs.
        nglyph = numel(regexprep(s, '[\\{}]|(?<=\\)[a-zA-Z]+', 'x'));
        lw = nglyph*ch_w/ax_w*xr;                   % [axis] label width
        px = Tq.mu_star(i);  py = Tq.sigma(i);
        gap = 0.012*xr;   % symbols sit close; a wider gap detaches them

        best = [];
        % COINCIDENT POINTS NEED VERTICAL TRAVEL, AND NOTHING ELSE DOES. On
        % the recovery response the recuperator returns exactly zero and feed
        % supply temperature a five-hundredth of the range, so the two share a
        % coordinate at plotting resolution and no choice of side parts them.
        % Travel is capped at three lines: enough to separate a coincidence,
        % short enough that the label still reads as belonging to its marker.
        for step = 0:3
            for dy = unique([0 step -step])
                for side = [1 -1]
                    if side > 0
                        x0 = px + gap;              x1 = x0 + lw;  ha = 'left';
                    else
                        x1 = px - gap;              x0 = x1 - lw;  ha = 'right';
                    end
                    y0 = py + dy*lh - 0.5*lh;       y1 = y0 + lh;
                    if x0 < 0 || x1 > xr || y0 < 0 || y1 > yr, continue, end
                    if isempty(placed)
                        hit = false;
                    else
                        hit = any(~(placed(:,1) > x1 | placed(:,2) < x0 | ...
                                    placed(:,3) > y1 | placed(:,4) < y0));
                    end
                    if ~hit
                        best = struct('x', tern(side>0, px+gap, px-gap), ...
                                      'y', py + dy*lh, 'ha', ha, ...
                                      'r', [x0 x1 y0 y1]);
                        break
                    end
                end
                if ~isempty(best), break, end
            end
            if ~isempty(best), break, end
        end
        if isempty(best)     % every candidate blocked: set it and accept
            best = struct('x', px+gap, 'y', py, 'ha', 'left', ...
                          'r', [px+gap px+gap+lw py-0.5*lh py+0.5*lh]);
        end
        placed(end+1,:) = best.r; %#ok<AGROW>

        % A LEADER IS DRAWN ONLY FOR A LABEL THAT MOVED, and almost none do.
        % Where two markers coincide the displaced symbol is otherwise
        % unattributable: the reader sees two names beside one dot and cannot
        % tell which is which. Where the label sits on its own marker no line
        % is drawn, so the plane stays clean.
        if abs(best.y - py) > 0.6*lh
            plot([px best.x], [py best.y], '-', 'Color', C_ref, ...
                 'LineWidth', 0.4, 'HandleVisibility','off');
        end

        % TEX IS ON HERE and off on the ranking chart: the plane sets symbols
        % that need subscripts and Greek, the ranking sets factor names that
        % contain no markup and would only risk being reinterpreted.
        text(best.x, best.y, s, 'FontSize',ST.ann, ...
             'HorizontalAlignment',best.ha, 'VerticalAlignment','middle', ...
             'Interpreter','tex');
    end

    xlabel(sprintf('\\mu^*   mean |elementary effect| in %s  [%s]', sym, unit), ...
           'FontSize',ST.lbl)
    ylabel(sprintf('\\sigma   spread of elementary effects  [%s]', unit), ...
           'FontSize',ST.lbl)
    legend('Location','northwest','Box','off','FontSize',ST.leg)
    xlim([0 xr]); ylim([0 yr]);
    % Caption sentence, not a title: above the diagonal a factor's influence
    % is contingent on the rest of the design, which is the finding.
    local_title_m(sprintf(['Morris screening of %s: factors above the diagonal ' ...
                           'act through interaction'], sym), ...
                  'FontSize',ST.ttl,'FontWeight','normal')

    %% ---- M2  mu* ranking with the sampling error -----------------------
    % Error bars are 2 x s.e.m. A bar whose interval reaches the axis is not
    % separable from zero at this number of trajectories, and two whose
    % intervals overlap are not separated from each other; the honest response
    % is more trajectories rather than a confident ordering.
    n = height(Tq);
    h(end+1) = newfig_m(sprintf('M2_morris_ranking_%s', tag), ...
                        FIG_RK, 0.72*n + 3.4); %#ok<AGROW>
    ax = axes; hold(ax,'on'); box(ax,'off'); grid(ax,'on'); set(ax,'Layer','top')
    yy = (n:-1:1).';

    % Alternating row bands, as on the sweep's tornado: with eleven unlabelled
    % rows the eye loses which label belongs to which bar.
    for i = 1:2:n
        patch([-9e9 9e9 9e9 -9e9], yy(i)+[-0.5 -0.5 0.5 0.5], ...
              [0.955 0.955 0.955],'EdgeColor','none','HandleVisibility','off');
    end

    barh(yy, Tq.mu_star, 0.58,'FaceColor',C_pt,'EdgeColor','none');
    errorbar(Tq.mu_star, yy, 2*Tq.sem,'horizontal','LineStyle','none', ...
             'Color',[.25 .25 .25],'CapSize',3,'LineWidth',0.8);

    % LEFT MARGIN SIZED FROM THE DRAWN CONTENT. The labels are long and set in
    % points, so a margin fixed as a fraction of the data range is either
    % wasteful or, for the factor whose interval reaches furthest left, too
    % small: gap height offset spans below zero and would otherwise print its
    % lower cap through its own name. The margin starts at the leftmost cap.
    xlo = min([0; Tq.mu_star - 2*Tq.sem]);
    xlab = xlo - 0.03*xm;
    for i = 1:n
        text(xlab, yy(i), char(Tq.label(i)),'HorizontalAlignment','right', ...
             'FontSize',ST.leg,'Interpreter','none');
        text(Tq.mu_star(i)+2*Tq.sem(i)+0.03*xm, yy(i), ...
             sprintf('\\sigma/\\mu^* = %.2f', Tq.sigma(i)/max(Tq.mu_star(i),eps)), ...
             'FontSize',ST.ann,'Color',C_ann);
    end

    % THE AXIS IS TICKED ONLY WHERE THE QUANTITY EXISTS. mu* is a mean of
    % absolute effects and cannot be negative, so a negative tick invites the
    % reader to interpret a region the statistic does not occupy. An interval
    % reaching past zero is still drawn, because that a factor is
    % indistinguishable from no effect at this many trajectories is the point.
    set(ax,'YTick',[],'YColor','none');
    xlim([xlab - 0.40*xm, 1.42*xm]); ylim([0.4 n+0.6]);
    tk = get(ax,'XTick');  set(ax,'XTick',tk(tk >= 0));
    xlabel(sprintf('\\mu^*   per unit of the swept range of each factor  [%s]', unit), ...
           'FontSize',ST.lbl)
    local_title_m(['Importance ranking; bars 2\times s.e.m., so overlapping ' ...
                   'intervals are not separated'], ...
                  'FontSize',ST.ttl,'FontWeight','normal')
end
end


function [sym, unit, tag] = morris_resp_label(resp)
%MORRIS_RESP_LABEL  Typeset symbol, unit and filename stem for a response.
%
% NO MATLAB IDENTIFIERS ON THE ARTWORK, for the reason resp_axislabel gives in
% Sensitivity.m: a label built from the column name prints a variable the
% reader has never seen. mu* and sigma both carry the unit of their response,
% and that unit is not the same for the two responses screened here, so it is
% declared beside the identifier it belongs to.
switch resp
    case 'SEC_external'
        sym = 'SEC_{external}';   unit = 'kWh m^{-3}';       tag = 'SEC';
    case 'R_still_pct'
        sym = '\itR\rm_{still}';  unit = 'percentage points'; tag = 'recovery';
    case 'prod_L_m2_day'
        sym = 'areal productivity'; unit = 'L m^{-2} d^{-1}'; tag = 'productivity';
    otherwise
        sym = strrep(resp,'_','\_');  unit = '';  tag = resp;
end
end


function h = newfig_m(name, w_cm, h_cm)
% One creation point, authored at final printed size in centimeters, matching
% newfig() in Sensitivity.m. Position in PIXELS would make the printed size
% depend on the display, which is what breaks a figure set across machines.
%
% THE WINDOW IS FORCED UNDOCKED. A docked figure takes the size of the dock,
% and the export reads Position to set PaperPosition, so with figures docked
% every panel in a set leaves at the same screen-shaped canvas whatever it was
% authored at. Fitting that canvas to a text column then scales 9 pt labels
% down by the width ratio, which is the usual reason a vector figure reaches
% the page looking soft. Position is re-asserted after the style change
% because undocking moves the window.
h = figure('Name',name,'Color','w','WindowStyle','normal');
set(h,'Units','centimeters','Position',[2 2 w_cm h_cm]);
end


function S = morris_style()
%MORRIS_STYLE  Type sizes and figure conventions, shared with the sweep.
%
% ONE DECLARATION WHEN THERE IS ONE, A COPY WHEN THERE IS NOT. The sweep holds
% these values in local_style() inside Sensitivity.m, where they govern every
% figure in that file. Reading them across means an edit there moves the Morris
% panels too, so the two sets cannot drift apart on the page; but the accessor
% that exposes them is a recent addition and a copy of Sensitivity.m without it
% would otherwise take the figures down. The values below mirror that
% declaration and are used only when it cannot be reached.
%
% Authored at final printed size: type is fixed in points while the canvas is
% fixed in centimeters, so a panel placed at other than its native width
% carries its labels onto the page at the wrong size.
S = struct('tick',8, 'lbl',9, 'ttl',9, 'leg',8, 'ann',7.5, 'ms',3.5, ...
           'titles',false, 'font','Arial');
try
    S = Sensitivity('style');
catch
    % Sensitivity.m predates the 'style' accessor; the mirror above stands in.
end
end


function s = morris_symbol(name, label)
%MORRIS_SYMBOL  The manuscript's symbol for a factor, for use on the artwork.
%
% THE PLANE IS LABELLED IN SYMBOLS, THE RANKING IN WORDS. Both carry all
% eleven factors, but the plane must set its labels inside the data area,
% where a four-word name is wider than the gap between neighbouring points and
% the weak factors are close enough to coincide. The symbols are those already
% declared in the nomenclature, so they cost the reader no lookup and need no
% legend; an arbitrary key (F1, F2, ...) would free the same space at the
% price of a cross-reference on every point.
%
% A name not listed here falls back to its full label rather than to a code,
% because an unlabelled or arbitrarily labelled point is worse than a wide one.
switch name
    case 'Qfan',         s = 'Q_{fan}';
    case 'A_plate',      s = 'A_{p}';
    case 'hcomp',        s = 'h_{comp}';
    case 'Vfeed',        s = 'V_{feed}';
    case 'Np',           s = 'N_{p}';
    case 'theta',        s = '\theta';
    case 'T_air_in',     s = 'T_{air,in}';
    case 'T_coil',       s = 'T_{coil}';
    case 'Tfeed_target', s = 'T_{feed,tgt}';
    case 'eps_recup',    s = '\epsilon_{recup}';
    case 'Tfeed',        s = 'T_{feed}';
    otherwise,           s = label;
end
end


function local_title_m(varargin)
%LOCAL_TITLE_M  Title only when the shared declaration asks for one.
%
% Journal figures carry a caption, not a title. The sentence stays in the code
% for whoever writes the caption, and the switch is the one that governs the
% sweep's figures, so the two sets cannot disagree about whether titles are
% drawn.
ST = morris_style();
if ST.titles, title(varargin{:}); end
end