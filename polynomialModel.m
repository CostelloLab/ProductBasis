% polynomialModel() converts a set of logic tables in given variables into
% an equivalent set of polynomials in all the variables

function M = polynomialModel(logicTables, varsOfInterest)

bMax = length(logicTables);

if ~exist('varsOfInterest', 'var')
    varsOfInterest = (eye(bMax) == 1);
end

M.xs = varsOfInterest;
M.fs = zeros(0, size(varsOfInterest, 1));
M.ts = ones(0, 1);

M.baseFs = cell(1, bMax);
for loopB = 1:bMax
    tempPoly = logicTableToPoly(logicTables{loopB}{2});
    fullCoefs = false(size(tempPoly{2}, 1), bMax);
    for loopX = 1:length(logicTables{loopB}{1})             % manually loop to avoid problems when two tempPoly{1} variables map to the same
        fullCoefs(:, logicTables{loopB}{1}(loopX)) = ...    % fullCoefs variable, as can happen with periodic boundary conditions (e.g. fly network)
            fullCoefs(:, logicTables{loopB}{1}(loopX)) | tempPoly{2}(:, loopX);
    end
    M.baseFs{loopB} = reducePoly({ tempPoly{1}, fullCoefs });
end

M.cxs = false(0, bMax);
M.cs = cell(0);
M.cts = ones(0, 1);
