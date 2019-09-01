% polynomialModel() converts a set of logic tables in given variables into
% an equivalent set of polynomials in all the variables

function M = polynomialModelT(logicTables, varsOfInterest)

bMax = length(logicTables);

if ~exist('varsOfInterest', 'var')
    varsOfInterest = (eye(bMax) == 1);
end

if isempty(logicTables)
    numGs = 0;
else
    numGs = size(logicTables{1}{3}, 2);
end

M.xs = varsOfInterest;
M.fs = cell(1, 0);
M.ferr = M.fs;
M.gs = zeros(0, numGs);

M.baseFs = cell(1, bMax);
for loopB = 1:bMax
    if isempty(logicTables{loopB}{3})
        logicTables{loopB}{3} = zeros(size(logicTables{loopB}{2}, 1), 0);
    end
    
    for loopG = 1:size(logicTables{loopB}{2}, 1)
        tempPoly = logicTableToTPoly(logicTables{loopB}{2}(loopG, :));  % get the polynomial in the 'local' variables
        xFactors = false(size(tempPoly{3}, 1), bMax);                   % now move to global variables..
        for loopX = 1:length(logicTables{loopB}{1})                         % (manually loop to avoid problems when two tempPoly{1} variables map
            xFactors(:, logicTables{loopB}{1}(loopX)) = ...                 % to the same fullCoefs variable, as with periodic boundaries (e.g. fly))
                xFactors(:, logicTables{loopB}{1}(loopX)) | tempPoly{3}(:, loopX);
        end
        oneGPoly = reducePolyT({ tempPoly{1}, tempPoly{2}, xFactors, repmat(logicTables{loopB}{3}(loopG, :), length(tempPoly{1}), 1) });
        if loopG == 1
            M.baseFs{loopB} = oneGPoly;
        else
            M.baseFs{loopB} = { [ M.baseFs{loopB}{1}; oneGPoly{1} ], [ M.baseFs{loopB}{2}; oneGPoly{2} ], ...
                [ M.baseFs{loopB}{3}; oneGPoly{3} ], [ M.baseFs{loopB}{4}; oneGPoly{4} ] };
        end
    end
end
