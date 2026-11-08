% Load the extracted outputs
JG = load('JG_border_toa.mat', 'Tcap', 'burned24');

% Check the output matrix dimensions
size(JG.Tcap)  % Expected: 1792 × 33758

% Select one border cell
c = 1;

% Sensitivity of capped arrival time at this cell
[Si, STi, Sij] = sobol(JG.Tcap(:,c));

% Sensitivity of burning within 24 hours at the same cell
[Si_burn, STi_burn, Sij_burn] = sobol(JG.burned24(:,c));