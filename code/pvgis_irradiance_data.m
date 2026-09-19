function irr = pvgis_irradiance_data()
% PVGIS_IRRADIANCE_DATA  Monthly-average hourly POA irradiance for
% multiple locations, fixed plane, slope = 21 deg, azimuth = 0 deg
% (south-facing), local time -- as pulled from PVGIS "Monthly data".
%
% Locations:
%   Dhaka  : 23.755 N, 90.393 E   (PVGIS-ERA5)
%   Saudi  : 25.624 N, 42.353 E   (PVGIS-SARAH3, Jan uses PVGIS-ERA5)
%
% USAGE:
%   irr = pvgis_irradiance_data();
%   [t, G] = irr.get('Dhaka', 'July');   % or irr.get('Saudi', 7)
%
% irr.t              : 1 x 24 vector, time-of-day [s since midnight]
% irr.data.<Loc>.G    : 12 x 24 matrix, irradiance [W/m^2], one row per month
% irr.month_names     : 1 x 12 cellstr, row labels matching G's rows
% irr.location_names  : 1 x N cellstr, valid location names
% irr.get(loc, month) : returns (t_irr_data, Gi_irr_data) for one
%                        location + month. loc is a name (case-insensitive);
%                        month is either a name (case-insensitive, full or
%                        3-letter abbreviation) or a number 1-12.

irr.t = [2700 6300 9900 13500 17100 20700 24300 27900 31500 35100 38700 42300 ...
         45900 49500 53100 56700 60300 63900 67500 71100 74700 78300 81900 85500];

irr.month_names = {'January','February','March','April','May','June', ...
                    'July','August','September','October','November','December'};

irr.location_names = {'Dhaka','Saudi'};

% ---- Dhaka (23.755 N, 90.393 E) ----
irr.data.Dhaka.G = [ ...
    0 0 0 0 0   0  0 167 397 611 774 862 840 764 625 438 217   5 0 0 0 0 0 0;  % January
    0 0 0 0 0   0  2 183 418 637 808 903 903 836 699 506 275  67 0 0 0 0 0 0;  % February
    0 0 0 0 0   0 41 229 452 657 816 904 904 835 697 506 282  79 0 0 0 0 0 0;  % March
    0 0 0 0 0   1 76 240 420 587 718 793 807 748 623 446 246  74 0 0 0 0 0 0;  % April
    0 0 0 0 0  10 82 205 340 468 570 636 669 619 516 373 216  74 4 0 0 0 0 0;  % May
    0 0 0 0 0  11 67 164 280 393 483 532 529 488 415 304 182  74 11 0 0 0 0 0; % June
    0 0 0 0 0   6 57 150 264 378 468 515 497 471 398 297 182  76 12 0 0 0 0 0; % July
    0 0 0 0 0   0 51 154 277 406 500 550 542 501 422 310 177  64 3 0 0 0 0 0;  % August
    0 0 0 0 0   0 52 175 320 461 559 605 594 551 444 305 153  35 0 0 0 0 0 0;  % September
    0 0 0 0 0   0 61 219 404 565 668 708 678 606 478 310 131   5 0 0 0 0 0 0;  % October
    0 0 0 0 0   0 39 235 455 647 775 824 781 683 533 334 125   0 0 0 0 0 0 0;  % November
    0 0 0 0 0   0  1 184 399 600 744 811 764 680 540 353 145   0 0 0 0 0 0 0]; % December

% ---- Saudi (25.624 N, 42.353 E) ----
irr.data.Saudi.G = [ ...
    0 0 0 0 0   0   0 119 351  571  747  852  879  825  691 498 264 31 0 0 0 0 0 0;  % January
    0 0 0 0 0   0   0  31 248  518  738  888  954  932  828 640 420 174 0 0 0 0 0 0;  % February
    0 0 0 0 0   0   0  87 320  575  783  937  996  963  862 665 451 213 6 0 0 0 0 0;  % March
    0 0 0 0 0   0   7 145 380  618  814  945  980  922  788 618 407 192 23 0 0 0 0 0; % April
    0 0 0 0 0  27 185 415 645  811  938  977  912  790  604 395 196 42 0 0 0 0 0 0;   % May
    0 0 0 0 0   0  39 211 438  655  827  948  985  955  843 675 469 249 67 0 0 0 0 0; % June
    0 0 0 0 0   0  28 191 417  639  812  945  988  958  847 681 476 259 74 0 0 0 0 0; % July
    0 0 0 0 0   0  11 163 398  638  828  962 1011  964  836 662 442 230 47 0 0 0 0 0; % August
    0 0 0 0 0   0   3 169 437  687  873 1002 1042  986  852 641 418 174 5 0 0 0 0 0;  % September
    0 0 0 0 0   0   0 161 433  677  860  969  996  930  776 573 338 91 0 0 0 0 0 0;   % October
    0 0 0 0 0   0   0  81 319  558  729  843  877  808  681 481 257 19 0 0 0 0 0 0;   % November
    0 0 0 0 0   0   0  25 261  493  679  804  842  798  685 507 271 19 0 0 0 0 0 0];  % December

irr.get = @(location, month) get_month(irr, location, month);

end % pvgis_irradiance_data

function [t_irr_data, Gi_irr_data] = get_month(irr, location, month)
% Resolve a location name + month name/number to a row of
% irr.data.<Location>.G, and return it alongside the shared time
% vector irr.t.

% --- resolve location ---
loc_idx = find(strcmpi(location, irr.location_names), 1);
if isempty(loc_idx)
    error('pvgis_irradiance_data:get:unknownLocation', ...
          'Unrecognized location: %s (valid: %s)', location, ...
          strjoin(irr.location_names, ', '));
end
loc_name = irr.location_names{loc_idx};
G = irr.data.(loc_name).G;

% --- resolve month ---
if ischar(month) || isstring(month)
    month = char(month);
    idx = find(strncmpi(month, irr.month_names, max(3,length(month))), 1);
    if isempty(idx)
        error('pvgis_irradiance_data:get:unknownMonth', ...
              'Unrecognized month name: %s', month);
    end
else
    idx = month;
    if idx < 1 || idx > size(G,1) || idx ~= round(idx)
        error('pvgis_irradiance_data:get:badIndex', ...
              'Month number must be an integer 1-12.');
    end
end

t_irr_data  = irr.t;
Gi_irr_data = G(idx,:);
end % get_month