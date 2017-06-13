% evolveTransitionMatrix() returns the time-evolution of a variable of the
% user's choice (bReadout), given an initial state b0.
% 
% xs, fs are the direct output of approxTransitionMatrix();
% x0 and xReadout are row vectors of the starting state and state to
% measure, respectively.

function R = evolveTransitionMatrix(M, initXCoefs, xReadout, numTimePoints, isContinuous)

if ~exist('isContinuous', 'var')
    isContinuous = false;
end

numXs = size(M.fs, 1);
R = zeros(1, numTimePoints+1);


    % do the time evolution manually..  (not by eigendecomposition -- that
    % will miss the generalized null eigenvectors)

if isempty(M.xs)
    return;
end

varPoly = { 1, xReadout };
while true
    whichXs = zeros(1, size(varPoly{1}, 1));
    for loopTerm = 1:size(varPoly{1}, 1)
        theX = find(sum(M.xs(1:numXs, :) ~= repmat(varPoly{2}(loopTerm, :), numXs, 1), 2) == 0, 1, 'first');
        if ~isempty(theX)
            whichXs(loopTerm) = theX;
        end
    end
    
    if sum(whichXs == 0) == 0
        break;
    end
    
    newPoly = { varPoly{1}(whichXs ~= 0, :), varPoly{2}(whichXs ~= 0, :) };
    
    for loopTerm = find(whichXs == 0)
        constrainedTerm = { varPoly{1}(loopTerm, :), varPoly{2}(loopTerm, :) };
        factoringConstraints = find(sum(M.cxs & ~repmat(varPoly{2}(loopTerm, :), size(M.cxs, 1), 1), 2) == 0)';
        if isempty(factoringConstraints)
            error('stuck!')
        end
        for loopConstraint = factoringConstraints
            constrainedTerm = multiplyPoly(constrainedTerm, M.cs{loopConstraint});
        end
        newPoly = { [ newPoly{1}; constrainedTerm{1} ], [ newPoly{2}; constrainedTerm{2} ] };
    end
    
    varPoly = reducePoly(newPoly);
end

if isempty(varPoly{1})
    return
end

xCoefs = initXCoefs;
for t = 1:numTimePoints+1
    R(t) = varPoly{1}' * xCoefs(whichXs);
    if isContinuous
        [ ~, xPts ] = ode45(@doIntegration, 0:1, xCoefs);
        xCoefs = xPts(end, :)';
    else
        xCoefs = M.fs * xCoefs;
    end
end


    function dx = doIntegration(~, xi)
        dx = M.fs*xi;
    end

end