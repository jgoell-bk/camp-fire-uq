function extract_border_toa(station)
% EXTRACT_BORDER_TOA  Arrival times at the Paradise wildland-border cells for every
% run of one station, in design_row order (the Sobol layout), saved to
% <station>_border_toa.mat in the current folder.
%
%   module load matlab/R2023b
%   cd /scratch/users/jgoell/elmfire-class/camp_sobol_d6    % or wherever the .mat files should go
%   matlab -nodisplay -batch "addpath('/oak/stanford/groups/gorle/jgoell/Wildfire/experiment/matlab'); extract_border_toa('JG')"
%
% Saved variables:
%   T        single, 1792 x 33758: arrival time (s after 06:00 PST 8 Nov 2018) at each
%            border cell in each run; NaN = not burned within 24 h
%   X        1792 x 6 inputs, columns a_U b_theta_deg b_1_pp b_10_pp b_100_pp q_spot
%   P_EPS    1792 x 1, = 0.001*50.^q_spot
%   block    1792 x 1 Sobol block of each run (A, AB_<factor>, BA_<factor>, B)
%   run_id   1792 x 1 run IDs (JG_0001 ...)
%   row,col  1-based raster row/column of each cell, in the column order of T
%            (MATLAB's column-major order: not the same order as Python's)
%
% Rows of T follow design_row, so reshape(T, 14, 128, []) gives (block, base row, cell)
% with blocks in the order A, AB_1..AB_6, BA_1..BA_6, B.

% Inputs. Groupmates can't read Oak: point these at the copies in the shared scratch folder.
mask_file = '/oak/stanford/groups/gorle/jgoell/Wildfire/experiment/qgis_setup/paradise_border_mask.tif';
runs_file = ['/oak/stanford/groups/gorle/jgoell/Wildfire/experiment/results/' station '_runs.csv'];

mask = imread(mask_file) == 1;                        % 1476 x 1457, true at the 33,758 border cells
[row, col] = find(mask);                              % same order as a(mask) below
runs = readtable(runs_file, 'TextType', 'string');
assert(isequal(runs.design_row, (1:height(runs))'), 'runs are not in design_row order')

nRuns = height(runs);
T = zeros(nRuns, nnz(mask), 'single');
for k = 1:nRuns
    a = imread(runs.toa_path(k));                     % arrival time raster of run k
    assert(isequal(size(a), size(mask)), 'raster %s is not on the mask grid', runs.toa_path(k))
    T(k, :) = a(mask);
end
T(T <= 0) = NaN;      % -9999 = not burned by 24 h (-1 marks the ignition cell, not in this region)

% New: arrival times capped at the simulation horizon
horizon_s = 24 * 3600;
Tcap = T;
Tcap(isnan(Tcap)) = horizon_s;
Tcap = min(Tcap, horizon_s);

% New: 1 = burned within 24 hours; 0 = did not burn
burned24 = double(~isnan(T) & T <= horizon_s);

X = runs{:, {'a_U', 'b_theta_deg', 'b_1_pp', 'b_10_pp', 'b_100_pp', 'q_spot'}};
P_EPS = runs.P_EPS;
block = runs.sobol_block;
run_id = runs.run_id;
out = [station '_border_toa.mat'];
save(out, 'T', 'Tcap', 'burned24', 'horizon_s', ...
    'X', 'P_EPS', 'block', 'run_id', 'row', 'col', '-v7.3');
fprintf('%s: %d runs x %d cells, %.0f%% of entries burned -> %s\n', ...
        station, nRuns, nnz(mask), 100 * mean(~isnan(T(:))), fullfile(pwd, out));
end
