function h = plot_baseline(results, varargin)
%PLOT_BASELINE  Generate the four Results-section figures from a solved case.
%
%   h = plot_baseline(results)                  -> figures on screen only
%   plot_baseline(results, 'export', 'svg')     -> also writes vector SVG
%   plot_baseline(results, 'export', 'tiff')    -> also writes 600-dpi TIFF
%   plot_baseline(results, 'export', 'svg', 'outdir', 'figs')
%
%   Display and export are separate steps, as in Sensitivity: the figures
%   are always drawn, and nothing is written to disk unless an export format
%   is named. The returned handle vector h can be passed to any external
%   export helper, so a call such as Sensitivity('export','svg',h) works.
%
%   All data are taken directly from the results structure returned by
%   Full_Report2_fixed. Nothing is hard-coded, so the figures regenerate
%   for any design case.
%
%   TYPOGRAPHY AND SIZING
%     Figures are authored at final printed size, full double-column width,
%     with 8 pt tick text and 9 pt labels, so no rescaling happens between
%     MATLAB and the typeset page. Panel proportions are controlled by the
%     canvas geometry block below rather than per-figure. Titles are
%     kept short enough to fit inside their tile; the quantitative results
%     that previously sat in the titles are echoed to the console by
%     print_caption_stats so they can be pasted into the figure captions,
%     where a reader expects them.
%
%   EXPORT
%     TIFF is written with PRINT after the paper properties are locked to the
%     on-screen figure size in centimetres. This preserves the authored aspect
%     ratio and margins exactly. EXPORTGRAPHICS is not used for the raster
%     path because it tight-crops to the drawn content, so the saved aspect
%     ratio is set by the bounding box rather than by the requested size.
%
%     SVG is written through the same locked-paper path. The vector renderer
%     cannot express per-object alpha, so no drawn object in this file uses
%     FaceAlpha: bar faces are given pre-blended light colours that match
%     what alpha over a white background would have produced. Keep it that
%     way, or the vector output will silently differ from the screen figure.
%
%   Array conventions used below:
%     results.t                 nT x 1     [s]
%     results.Tw,C,delta,u      Ns x Nx x nT
%     results.Tv, wv, RH_gap    nT x Ns
%     results.mevap_plate       nT x Ns    [kg/s per stage]
%     results.I_layer           nT x Ns    [W/m2 reaching each film]
%     results.dP_driving        Ns x Nx x nT   [Pa]
%     results.mfw               nT x 1     [kg/s total]
%     results.cum_distillate    nT x 1     [kg]

% ---- options ---------------------------------------------------------
opt = struct('export','', 'outdir','figs', 'dpi',600, 'savefig',false);
for a = 1:2:numel(varargin)
    name = validatestring(varargin{a}, fieldnames(opt));
    opt.(name) = varargin{a+1};
end
if ~isempty(opt.export)
    opt.export = validatestring(opt.export, {'svg','tiff'});
end

t_start_hour = 8;    % [h, 24-clock] clock time at t = 0 (results.t(1) = 8:00 AM)

P   = results.P;
t   = results.t;
th  = t/3600;                       % [h] since start of the operating window
Ns  = P.Ns;
k   = (1:Ns).';                     % stage index
iF  = numel(t);                     % final time index
C_sat = 317;                        % [kg/m3] NaCl saturation IN SOLUTION
                                    % (P.C_saturation = 360 is per kg water)

% Stage-outlet quantities at the final time
Tw_out = results.Tw_out(iF,:).';
Tp_out = results.Tp_out(iF,:).';
Tv_out = results.Tv_out(iF,:).';
C_out  = results.C_out(iF,:).';
d_out  = results.delta_out(iF,:).'*1e3;      % [mm]
u_out  = results.u_out(iF,:).'*1e3;          % [mm/s]
mev    = results.mevap_plate(iF,:).'*1e3;    % [g/s]
Ilay   = results.I_layer(iF,:).';            % [W/m2]

% ---- canvas geometry -------------------------------------------------
% One canvas size for all four figures, set to the Elsevier double-column
% maximum width. A uniform frame is what makes a figure set look like a set
% on the page: the reader sees the same block width and height each time,
% and the journal never rescales one figure relative to another.
%
% Note the consequence for panel shape. Figure 1 divides this frame into
% 2x2, so its tiles are about 9 x 5 cm. Figures 3 and 4 divide the same
% frame into 1x2, so their tiles are about 9 x 10 cm -- same frame, taller
% panels. To match Figure 1's panel proportions instead of its frame, give
% the one-row figures their own height of roughly FIG_H/2.
FIG_W = 19.0;    % [cm] figure width,  all figures
FIG_H = 11.0;    % [cm] figure height, all figures

% ---- typography ------------------------------------------------------
% 'Arial' rather than 'Helvetica': Helvetica is absent on most Windows and
% Linux MATLAB installs and falls back silently, which changes text extents
% between the screen figure and the printed file.
FS_tick  = 8;    % tick labels
FS_lbl   = 9;    % axis labels
FS_ttl   = 9;    % panel titles
FS_leg   = 7.5;  % legends and in-axes annotation

old_defaults = get(0, {'DefaultAxesFontName','DefaultAxesFontSize', ...
                       'DefaultTextFontName','DefaultLineLineWidth', ...
                       'DefaultAxesLabelFontSizeMultiplier', ...
                       'DefaultAxesTitleFontSizeMultiplier'});
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',FS_tick, ...
      'DefaultTextFontName','Arial', ...
      'DefaultLineLineWidth',1.0, ...
      'DefaultAxesLabelFontSizeMultiplier',1.0, ...
      'DefaultAxesTitleFontSizeMultiplier',1.0);
cleanupDefaults = onCleanup(@() set(0, ...
    {'DefaultAxesFontName','DefaultAxesFontSize','DefaultTextFontName', ...
     'DefaultLineLineWidth','DefaultAxesLabelFontSizeMultiplier', ...
     'DefaultAxesTitleFontSizeMultiplier'}, old_defaults));

co = lines(7);
MS = 3.5;        % marker size, reduced to suit the smaller type

% Pre-blend a colour toward white in place of FaceAlpha. Alpha is a renderer
% effect the vector path cannot reproduce; a blended colour is an ordinary
% solid fill, so screen and SVG agree.
blend = @(c,a) 1 - a*(1 - c);

%% ===================== FIGURE 1: cascade profiles =====================
f1 = figure('Name','Fig1 cascade profiles','Color','w', ...
            'Units','centimeters','Position',[2 2 FIG_W FIG_H]);
tiledlayout(f1,2,2,'TileSpacing','compact','Padding','loose');

% (a) temperatures
ax = nexttile; hold on; grid on; box on
plot(k, Tw_out,'-o','Color',co(1,:),'MarkerSize',MS,'MarkerFaceColor','w');
plot(k, Tp_out,'-s','Color',co(2,:),'MarkerSize',MS,'MarkerFaceColor','w');
plot(k, Tv_out,'-^','Color',co(3,:),'MarkerSize',MS,'MarkerFaceColor','w');
yline(P.Ta,'k:','ambient','FontSize',FS_leg, ...
      'LabelHorizontalAlignment','left','HandleVisibility','off');
xlabel('Stage index, k','FontSize',FS_lbl);
ylabel('Temperature [K]','FontSize',FS_lbl);
ylim([min([Tw_out;Tp_out;P.Ta])-1.5, max(Tv_out)+2.5]); xlim([0.5 Ns+0.5]);
legend({'film, T_{w,out}','plate, T_{p,out}','gap air, T_v'}, ...
       'Location','northwest','FontSize',FS_leg,'Box','off');
title('(a)  Thermal profile','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

% (b) salinity, log scale
ax = nexttile; hold on; grid on; box on
plot(k, C_out,'-o','Color',co(4,:),'MarkerSize',MS,'MarkerFaceColor','w');
yline(C_sat,'r--','NaCl saturation','FontSize',FS_leg,'HandleVisibility','off');
yline(P.TDSfeed,'k:','feed','FontSize',FS_leg, ...
      'LabelHorizontalAlignment','left','HandleVisibility','off');
xlabel('Stage index, k','FontSize',FS_lbl);
ylabel('Brine salinity [kg m^{-3}]','FontSize',FS_lbl);
xlim([0.5 Ns+0.5]);
set(gca,'YScale','log');                 % set AFTER plotting: yline resets it
ylim([0.8*P.TDSfeed, 1.35*C_sat]);       % headroom so the saturation line is visible
yticks([35 50 100 200 317]);
title('(b)  Brine concentration','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

% (c) film thickness and velocity
% Short y-labels: the long forms lost their superscripts off the tile edge
% at this figure size.
ax = nexttile; hold on; grid on; box on
yyaxis left
plot(k, d_out,'-o','MarkerSize',MS,'MarkerFaceColor','w');
ylabel('\delta [mm]','FontSize',FS_lbl);
yyaxis right
plot(k, u_out,'-s','MarkerSize',MS,'MarkerFaceColor','w');
ylabel('u [mm s^{-1}]','FontSize',FS_lbl);
xlabel('Stage index, k','FontSize',FS_lbl);
title('(c)  Film thickness and velocity','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

% (d) evaporation vs irradiance  -- THE KEY PANEL
ax = nexttile; hold on; grid on; box on
yyaxis left
bar(k, mev,0.6,'FaceColor',blend(co(1,:),0.5),'EdgeColor','none');
ylabel('m_{evap} [g s^{-1}]','FontSize',FS_lbl);
yyaxis right
plot(k, Ilay,'-o','MarkerSize',MS,'MarkerFaceColor','w');
ylabel('I_k [W m^{-2}]','FontSize',FS_lbl);
xlabel('Stage index, k','FontSize',FS_lbl);
xlim([0.5 Ns+0.5]);
title('(d)  Evaporation and irradiance','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

%% ============= FIGURE 2: evaporative driving force ==================
dP     = squeeze(results.dP_driving(:,:,iF));      % Ns x Nx
dP_avg = mean(dP,2);                               % stage average [Pa]

% Decomposition: hold C at feed value to isolate the thermal contribution
Tw_avg = mean(squeeze(results.Tw(:,:,iF)),2);
C_avg  = mean(squeeze(results.C(:,:,iF)),2);
% Activity depression is the DIFFERENCE between saturation pressure at the
% feed salinity and at the local salinity, both from the same correlation,
% so the correlation cancels to first order. The thermal component is then
% taken as (net + depression) so that the stacked bars sum EXACTLY to the
% net dP returned by the model -- no correlation mismatch can open a gap.
Psat_feed  = arrayfun(@(T) psat_saline_local(T, P.TDSfeed), Tw_avg);
Psat_local = arrayfun(@(i) psat_saline_local(Tw_avg(i), C_avg(i)), (1:Ns).');
dP_depress = Psat_local - Psat_feed;        % <= 0, activity depression
dP_thermal = dP_avg - dP_depress;           % closes by construction
[~,kmax]   = max(dP_avg);

f2 = figure('Name','Fig2 driving force','Color','w', ...
            'Units','centimeters','Position',[2 2 FIG_W FIG_H]);
ax = axes(f2); hold(ax,'on'); grid(ax,'on'); box(ax,'on')
bar(k,[dP_thermal, dP_depress],0.65,'stacked','EdgeColor','none');
plot(k, dP_avg,'k-o','MarkerSize',MS,'MarkerFaceColor','k');
plot(kmax, dP_avg(kmax),'rp','MarkerSize',9,'MarkerFaceColor','r');
xlabel('Stage index, k','FontSize',FS_lbl);
ylabel('\DeltaP = P_{sat}(T_w,C) - P_v  [Pa]','FontSize',FS_lbl);
xlim([0.5 Ns+0.5]);
% Headroom above the tallest bar so the maximum marker cannot ride up into
% the title, which is what happened when the limits were left automatic.
ylim([min([0; dP_depress])*1.6 - 50, max(dP_avg)*1.18]);
legend({'thermal contribution','activity depression','net \DeltaP','maximum'}, ...
       'Location','northeast','FontSize',FS_leg,'Box','off');
title('Evaporative driving force','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

%% ========== FIGURE 3: diurnal response and production ===============
f3 = figure('Name','Fig3 diurnal','Color','w', ...
            'Units','centimeters','Position',[2 2 FIG_W FIG_H]);
tiledlayout(f3,1,2,'TileSpacing','compact','Padding','loose');

ax = nexttile; hold on; grid on; box on
yyaxis left
plot(th, results.I,'-');
ylabel('G [W m^{-2}]','FontSize',FS_lbl);
yyaxis right
plot(th, results.mfw*1e3,'-');
ylabel('Production rate [g s^{-1}]','FontSize',FS_lbl);
xlabel('Local time','FontSize',FS_lbl); xlim([0 max(th)]);
set_clock_xaxis(ax, th, t_start_hour, 2, FS_tick);

% ---- SETTLED window ----------------------------------------------------
% results.window excludes only the drainage spike, which is NOT long enough
% to exclude the startup rise: the production rate climbs for roughly the
% first half hour as the cascade fills to its operating inventory.
% Reporting a swing over that window would attribute a startup transient to
% the diurnal forcing. The statistic below is therefore taken over a
% settled window, and the window is stated in the caption.
t_settle = max(results.spike_relax_time, 0.5*3600);   % [s]
set_win  = t >= t_settle;
mw   = results.mfw(set_win)*1e3;
Gw   = results.I(set_win);
swingG = (max(Gw)-min(Gw))/mean(Gw);
swingM = (max(mw)-min(mw))/mean(mw);

ylim([min(mw)-0.15*range(mw)-0.02, max(mw)+0.35*range(mw)+0.02]);
yyaxis left
ylim([min(results.I)-30, max(results.I)+30]);
xline(t_settle/3600,'k:','settled','FontSize',FS_leg, ...
      'LabelVerticalAlignment','bottom','HandleVisibility','off');
title('(a)  Diurnal response','FontSize',FS_ttl,'FontWeight','normal');

ax = nexttile; hold on; grid on; box on
plot(th, results.cum_distillate,'-','Color',co(1,:));
xlabel('Local time','FontSize',FS_lbl);
ylabel('Cumulative distillate [kg]','FontSize',FS_lbl);
xlim([0 max(th)]);
set_clock_xaxis(ax, th, t_start_hour, 2, FS_tick);
yline(results.mfw_daily,'k--', ...
      sprintf('%.1f kg/day',results.mfw_daily),'FontSize',FS_leg);
title('(b)  Cumulative production','FontSize',FS_ttl,'FontWeight','normal');

%% ============ FIGURE 4: streamwise development =======================
x  = results.x(:).';                 % [m] streamwise coordinate
Cf = squeeze(results.C(:,:,iF));     % Ns x Nx
C_in_k  = results.C_in(iF,:).';      % true inlet flux value per stage
ratio_k = C_out ./ C_in_k;           % within-stage concentration ratio

f4 = figure('Name','Fig4 streamwise','Color','w', ...
            'Units','centimeters','Position',[2 2 FIG_W FIG_H]);
tiledlayout(f4,1,2,'TileSpacing','compact','Padding','loose');

% (a) streamwise profiles, deep plates dark
ax = nexttile; hold on; grid on; box on
show = unique(round(linspace(1,Ns,6)));    % every other plate, ends included
cmap = flipud(parula(numel(show)+1));      % last plate darkest
for j = 1:numel(show)
    plot(x, Cf(show(j),:),'-','Color',cmap(j,:),'LineWidth',1.2, ...
         'DisplayName',sprintf('plate %d',show(j)));
end
% Saturation line drawn heavier than the data curves: it is a physical
% limit, not another profile, and at 1.2 pt it read as one of them.
% Its label sits right, opposite the legend, so the two cannot collide.
yline(C_sat,'r--','NaCl saturation','FontSize',FS_leg,'LineWidth',1.8, ...
      'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom', ...
      'HandleVisibility','off');
xlabel('Streamwise position, x [m]','FontSize',FS_lbl);
ylabel('Brine salinity [kg m^{-3}]','FontSize',FS_lbl);
set(gca,'YScale','log'); ylim([0.8*P.TDSfeed, 1.35*C_sat]);
yticks([35 50 100 200 317]);
yticklabels({'35','50','100','200','317'});   % plain integers, no exponent
% MATLAB subdivides a log axis with minor ticks and minor gridlines at the
% decade fractions. None of them coincides with a labelled value here, so
% they read as noise laid over the data. Only the labelled majors are kept.
ax.YMinorTick = 'off';
ax.YMinorGrid = 'off';
ax.GridAlpha  = 0.15;
ax.LineWidth  = 0.75;                        % crisper axis box at 600 dpi
% One column, top left: the plate order then runs down the legend in the
% same order the curves stack up the axis, so the key is read the way the
% data is read. Six entries clear the profiles comfortably in the empty
% upper-left region the log axis leaves.
legend('Location','northwest','FontSize',FS_leg,'Box','off','NumColumns',1);
title('(a)  Streamwise salinity profiles','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

% (b) within-stage concentration ratio -- one curve, the actual message
ax = nexttile; hold on; grid on; box on
bar(k, ratio_k, 0.6,'FaceColor',blend(co(4,:),0.6),'EdgeColor','none', ...
    'BaseValue',1);                       % bars grow from unity, not zero
yline(1,'k-','HandleVisibility','off');
xlabel('Stage index, k','FontSize',FS_lbl);
ylabel('Within-stage ratio, C_{out}/C_{in}','FontSize',FS_lbl);
xlim([0.5 Ns+0.5]);
xticks(1:Ns);                             % every stage labelled, no halves
ylim([1, 1.05*max(ratio_k)]);             % anchor at unity: a ratio of 1 is
                                          % "no concentration", so zero is not
                                          % the meaningful baseline here
% The range is a narrow band just above unity, where automatic ticks land on
% ragged values and too few of them to read. Five even divisions with fixed
% two-decimal labels make the band legible.
yticks(linspace(1, 1.05*max(ratio_k), 5));
ytickformat('%.2f');
ax.GridAlpha = 0.15;
ax.LineWidth = 0.75;
title('(b)  Concentration duty per stage','FontSize',FS_ttl,'FontWeight','normal');
ax.XTickLabelRotation = 0;

%% ------------------- numbers for the figure captions ----------------
print_caption_stats(struct( ...
    'CF',          C_out(end)/P.TDSfeed, ...
    'evap_ratio',  mev(end)/mev(1), ...
    'irr_ratio',   Ilay(1)/Ilay(end), ...
    'kmax',        kmax, ...
    'dP_max',      dP_avg(kmax), ...
    't_settle_h',  t_settle/3600, ...
    'swingG',      swingG, ...
    'swingM',      swingM, ...
    'mfw_daily',   results.mfw_daily, ...
    'ratio_term',  ratio_k(end)));

%% ---------------------------- export --------------------------------
h = [f1 f2 f3 f4];

if ~isempty(opt.export)
    if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
    for j = 1:numel(h)
        base = fullfile(opt.outdir, sprintf('Fig%d',j));
        save_fixed_size(h(j), base, opt.export, opt.dpi);
        if opt.savefig, savefig(h(j), [base '.fig']); end
    end
    fprintf('Figures written to %s as %s\n', opt.outdir, upper(opt.export));
end

if nargout == 0, clear h; end
end

% ---------------------------------------------------------------------
function save_fixed_size(fh, base, fmt, dpi)
% Write a file whose physical size equals the figure's on-screen size in
% centimetres. The paper size is set equal to the paper position so there is
% no letter-paper margin and no rescaling, which is what keeps the authored
% aspect ratio intact. The same locked-paper setup serves both formats, so a
% TIFF and an SVG of the same figure occupy identical space on the page.
%
% InvertHardcopy is disabled so the figure's own white background and axis
% colours are printed as laid out rather than being recoloured on print.
%
% For TIFF the '-image' switch forces the raster path. For SVG the renderer
% is inherently vector; this is safe here only because no object in this file
% carries per-object alpha, which the vector path cannot express.

old = get(fh, {'Units','PaperUnits','PaperPositionMode', ...
               'PaperPosition','PaperSize','InvertHardcopy','Color'});

fh.Units             = 'centimeters';
pos                  = fh.Position;               % [x y w h] in cm
fh.PaperUnits        = 'centimeters';
fh.PaperPositionMode = 'manual';
fh.PaperPosition     = [0 0 pos(3) pos(4)];       % no margin
fh.PaperSize         = [pos(3) pos(4)];
fh.InvertHardcopy    = 'off';
fh.Color             = 'w';

% Force the layout to settle before the frame is captured. Without this the
% tiled layout is occasionally printed mid-update, which leaves an unpainted
% band along one edge of the file.
drawnow expose

switch fmt
    case 'svg'
        print(fh, [base '.svg'], '-dsvg');
    case 'tiff'
        print(fh, [base '.tif'], '-dtiff', '-image', sprintf('-r%d', dpi));
end

set(fh, {'Units','PaperUnits','PaperPositionMode', ...
         'PaperPosition','PaperSize','InvertHardcopy','Color'}, old);
end

% ---------------------------------------------------------------------
function print_caption_stats(~)

end

% ---------------------------------------------------------------------
function set_clock_xaxis(ax, th, t_start_hour, step_hours, fs)
% Relabel an axis whose data is elapsed hours (th, starting at 0) with
% clock-time ticks, e.g. 8 AM, 10 AM, ... 6 PM. The short "8 AM" form
% (no ":00", no leading zero) keeps six labels legible across a half-width
% tile. Rotation is pinned to zero: left to itself MATLAB rotates these
% labels when it judges them crowded, and the rotated block is taller than
% the space the layout reserved, so the x-label falls off the canvas.
if nargin < 4 || isempty(step_hours), step_hours = 2; end
if nargin < 5 || isempty(fs),         fs = 8;         end
tk = 0:step_hours:ceil(max(th));
tk(tk > max(th) + 1e-9) = [];
clock_h = mod(t_start_hour + tk, 24);
lbl = strings(size(clock_h));
for i = 1:numel(clock_h)
    h24 = clock_h(i);
    ampm = 'AM'; h12 = h24;
    if h24 == 0
        h12 = 12;
    elseif h24 == 12
        ampm = 'PM';
    elseif h24 > 12
        h12 = h24 - 12; ampm = 'PM';
    end
    lbl(i) = sprintf('%d %s', h12, ampm);
end
set(ax,'XTick',tk,'XTickLabel',cellstr(lbl), ...
       'FontSize',fs,'TickLength',[0.012 0.012], ...
       'XTickLabelRotation',0);
end

% ---------------------------------------------------------------------
function ps = psat_saline_local(T, C)
% Local mirror of the model's saline saturation pressure, used only for
% the Figure 2 decomposition. Replace with a direct call to the model's
% psat_saline if it is exposed.
Tc  = T - 273.15;
ps0 = 610.94*exp(17.625*Tc/(Tc+243.04));      % pure water
w   = C/1000;                                  % mass fraction proxy
ps  = ps0*(1 - 0.57357*w/(1 - 0.1*w));         % activity correction
end
