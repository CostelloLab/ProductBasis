% evolveTransitionMatrix() returns the time-evolution of a variable of the
% user's choice (bReadout), given an initial state b0.
% 
% xs, fs are the direct output of approxTransitionMatrix();
% x0 and xReadout are row vectors of the starting state and state to
% measure, respectively.

function [ R, extraVar ] = evolveF(M, initXCoefs, xReadout, numTimePoints, isContinuous)

if ~exist('isContinuous', 'var')
    isContinuous = false;
end

numBs = length(M.baseFs);
isPBN = false;
for loopB = 1:numBs
    if sum(M.baseFs{loopB}{1} ~= round(M.baseFs{loopB}{1} )) > 0
        isPBN = true;
        break;
    end
end

numXs = size(M.fs{1}, 1);       % change!!
if ~isempty(xReadout)
    R = zeros(1, numTimePoints+1);
else
    R = zeros(size(M.xs, 1), numTimePoints+1);
end
extraVar = [];


    % do the time evolution manually..  (not by eigendecomposition -- that
    % will miss the generalized null eigenvectors)

if isempty(M.xs)
    return;
end

if ~isempty(xReadout)
    varPoly = { 1, 0, xReadout };
    while true
        whichXs = zeros(1, size(varPoly{1}, 1));
        for loopTerm = 1:size(varPoly{1}, 1)
            theX = find(sum(M.xs(1:numXs, :) ~= repmat(varPoly{3}(loopTerm, :), numXs, 1), 2) == 0, 1, 'first');
            if ~isempty(theX)
                whichXs(loopTerm) = theX;
            end
        end
        
        if sum(whichXs == 0) == 0
            break;
        end
        
        newPoly = { varPoly{1}(whichXs ~= 0, :), varPoly{2}(whichXs ~= 0, :), varPoly{3}(whichXs ~= 0, :) };
        
        for loopTerm = find(whichXs == 0)
            extraVar = varPoly{3}(loopTerm, :);
            return
        end
        
        varPoly = reducePoly(newPoly);
    end
    
    if isempty(varPoly{1})
        return
    end
end

xCoefs = initXCoefs;
for t = 1:numTimePoints+1
    if ~isempty(xReadout)
        R(t) = varPoly{1}' * xCoefs(whichXs);
    else
        R(:, t) = xCoefs';
    end
    if isContinuous
        [ ~, xPts ] = ode45(@doIntegration, 0:1, xCoefs);
        xCoefs = xPts(end, :)';
    else
        xCoefs = M.fs{1} * xCoefs;
    end
end


    function dx = doIntegration(~, xi)
        dx = M.fs{1}*xi;
    end

end