%% Sobol indices from already-extracted ELMFIRE outputs
% Put this script and sobol.m together. Set data_folder below.
code_folder = fileparts(mfilename('fullpath'));
addpath(code_folder);
data_folder = pwd;  % EDIT if the .mat file is in another folder.

% Extraction only needs to be done once after updating the extractor:
% extract_border_toa('JG');
% extract_border_toa('OP');
% extract_border_toa('CO');

station = 'JG';
% Load selected variables to avoid also loading the large original T matrix.
JG = load(fullfile(data_folder, [station '_border_toa.mat']), ...
    'Tcap', 'burned24', 'block', 'row', 'col');
assert(isfield(JG, 'Tcap') && isfield(JG, 'burned24'), ...
    'Run the updated extractor first to save Tcap and burned24.');
Tcap = JG.Tcap;
burned24 = JG.burned24;
disp(size(Tcap));   % Expected: 1792-by-number_of_border_cells.

% Verify repeated block ordering and matching AB/BA factor names.
labels = reshape(string(JG.block), 14, []);
assert(size(labels, 2) == 128 && size(Tcap, 1) == 1792, ...
    'This main script expects the N=128, D=6 design.');
assert(all(labels == labels(:, 1), 'all'), ...
    'Block ordering differs between batches.');
assert(labels(1,1) == "A" && labels(14,1) == "B" && ...
    all(startsWith(labels(2:7,1), "AB_")) && ...
    all(startsWith(labels(8:13,1), "BA_")) && ...
    all(extractAfter(labels(2:7,1), "AB_") == ...
        extractAfter(labels(8:13,1), "BA_")), ...
    'Expected A, six matching AB and BA factors, then B.');
factor_names = extractAfter(labels(2:7,1), "AB_");

%% Start with ONE border cell
c = 1;  % Change this to inspect another column of Tcap.
[Si, STi, Sij] = sobol(Tcap(:, c));
disp(table(factor_names, Si(:), STi(:), ...
    'VariableNames', {'Factor', 'FirstOrder', 'TotalEffect'}));
disp('Second-order matrix (factor order matches the table):');
disp(Sij);
% Alternative QoI at the same location:
[Si_burn, STi_burn, Sij_burn] = sobol(burned24(:, c));

%% Entire border: uncomment this section when ready
nCells = size(Tcap, 2);
Si_all = nan(nCells, 6);
STi_all = nan(nCells, 6);
Sij_all = nan(6, 6, nCells);
Si_burn_all = nan(nCells, 6);
STi_burn_all = nan(nCells, 6);
Sij_burn_all = nan(6, 6, nCells);
for c = 1:nCells
    [Si_all(c,:), STi_all(c,:), Sij_all(:,:,c)] = sobol(Tcap(:,c));
    [Si_burn_all(c,:), STi_burn_all(c,:), Sij_burn_all(:,:,c)] = ...
        sobol(burned24(:,c));
end
row = JG.row;
col = JG.col;
save(fullfile(data_folder, [station '_sobol_results.mat']), ...
    'Si_all', 'STi_all', 'Sij_all', ...
    'Si_burn_all', 'STi_burn_all', 'Sij_burn_all', ...
    'factor_names', 'row', 'col', '-v7.3');
