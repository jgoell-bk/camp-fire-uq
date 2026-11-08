function [Si, STi, Sij] = sobol(y)
% Input:
%   y: outputs at ONE border cell, in design_row order.
%      Each batch contains:
%      A, AB_1..AB_6, BA_1..BA_6, B
%
% Outputs:
%   Si:  1 × 6 first-order indices
%   STi: 1 × 6 total-effect indices
%   Sij: 6 × 6 pairwise interaction indices

validateattributes(y, {'numeric', 'logical'}, ...
    {'vector', 'nonempty', 'real', 'finite'});

p = 6;
runs_per_batch = 2*p + 2;  % 14

assert(mod(numel(y), runs_per_batch) == 0, ...
    'Expected 14 outputs per complete batch.');

% Use a column vector in double precision.
% Centering improves numerical stability.
y = double(y(:));
y = y - mean(y);

M = numel(y) / runs_per_batch;  % 128 for your design

% Each ROW now contains one batch's 14 outputs
Y = reshape(y, runs_per_batch, M).';

% Separate the reference and hybrid outputs
ya  = Y(:,1);       % 128 × 1
yAB = Y(:,2:7);     % 128 × 6
yBA = Y(:,8:13);    % 128 × 6
yb  = Y(:,14);      % 128 × 1

% Estimate output variance using A and B together
V = var([ya; yb], 1);

% Allocate results
Si  = nan(1,p);
STi = nan(1,p);
Sij = nan(p,p);

% Normalized indices are undefined for constant outputs
if V == 0
    return
end

% First-order and total-effect indices
for i = 1:p
    Si(i) = mean(yb .* (yAB(:,i) - ya)) / V;

    STi(i) = 0.5 * mean((ya - yAB(:,i)).^2) / V;
end

% Second-order indices
for i = 1:p-1
    for j = i+1:p
        joint = mean(yBA(:,i) .* yAB(:,j) - ya .* yb) / V;

        Sij(i,j) = joint - Si(i) - Si(j);
        Sij(j,i) = Sij(i,j);
    end
end

end