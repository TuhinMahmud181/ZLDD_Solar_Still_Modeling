function FC = CFD_Comparison_Input_Variable()
%CFD_COMPARISON_INPUT_VARIABLE  Input set of the CFD comparison case.
%
%   FC = CFD_Comparison_Input_Variable() returns the configuration used for the
%   comparison with the ANSYS Fluent model in Section 6 of the article. It is
%   the baseline design case of Baseline_Input_Variable() with two changes:
%
%     1. Feed flow of 115 L per operating window (baseline: 100 L).
%     2. Entrance and turn factors on the gap heat- and mass-transfer
%        coefficients, FC.gap_factor = [top gap, plate gaps 2..Np, floor gap]:
%          top gap, under the glass : 5.6  correlation evaluated on the 50 mm
%                                          inlet slot rather than on the full
%                                          height of the space under the glass
%          plate gaps 2..Np         : 1.7  each entered through a sharp 180-deg
%                                          turn (Chyu, 1991)
%          floor gap                : 1.2  straight entry, entrance-region
%                                          enhancement only (Mills, 1962)
%
%   All results in Sections 4 and 5 of the article use Baseline_Input_Variable(),
%   which carries no gap factors, so the fully developed correlation is used
%   unchanged there. The baseline file is not edited for this case: the
%   returned struct is copied and changed here, as Baseline_Input_Variable
%   prescribes for any study.
%
%   USAGE
%       FC      = CFD_Comparison_Input_Variable();
%       results = ZLDD_complete_modeling(FC);
%
%   The exported results of the Fluent model for the same case are in the
%   cfd/ folder of the repository.
%
%   See also BASELINE_INPUT_VARIABLE, ZLDD_COMPLETE_MODELING.

FC = Baseline_Input_Variable();

FC.Vfeed      = 115;                % [L/day] feed per operating window (baseline 100)
FC.gap_factor = [5.6, 1.7, 1.2];    % [-] [top gap, plate gaps 2..Np, floor gap]
end
