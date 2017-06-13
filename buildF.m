% buildF() approximates the transition rules of a Boolean
% network.  Each variable 'x' can be one of the Boolean variables 'b' or a
% product of them:  for example x13 is b1*b3.  A transition function is
% denoted 'f'.

function [ M, simStartTime, stats ] = ...
    buildF(M, maxVars, maxFinalVars, eigThreshold, isContinuous)


    % our initial features are just the fundamental Boolean state variables

backspaceString = '';
simStartTime = 1;
finished = false;

numBs = size(M.xs, 2);
stats = zeros(0, 4);

if ~exist('maxVars', 'var')
    maxVars = min(2^numBs, 1e4);
end

if ~exist('maxFinalVars', 'var')
    maxFinalVars = maxVars;
end

if ~exist('eigThreshold', 'var')
    eigThreshold = 0.01;
end

if ~exist('isContinuous', 'var')
    isContinuous = false;
end

xsToKeep = true(1, size(M.xs, 1));
isQueryVar = xsToKeep;

    % first iteration:  find evolution equations for the fundamental variables
    % 
    % 2nd-Nth iterations:  find evolution equations for the higher-order
    % correlation variables that were introduced in earlier iterations
    
while diff(size(M.fs)) ~= 0 || size(M.xs, 1) > maxFinalVars
    
    stats(end+1, :) = [ 0 0 0 0 ];
    
    
        % PART 1
        % weed out the linearly-dependent variables if we are ignoring the
        % transient
    
    removeMoreEigs = true;
    while size(M.fs, 1) > maxFinalVars && removeMoreEigs
        
        printProgress('finding dependencies..');
        
        idx = find(M.ts < simStartTime);
        for loopFindRedundancies = 1:2
            littleF = M.fs(idx, :);
            little = max(1.e-6, max(size(littleF)) * eps(norm(littleF)) );
            
            squareF = M.fs(idx, idx);
            [ eigenVecs, eigenVals ] = eig(squareF' + 1.e-12*eye(size(squareF)), 'nobalance');
            eigenVals = diag(eigenVals)';
            
            littleTs = repmat(M.ts(idx), 1, length(idx));
            littleTs(abs(eigenVecs) < little) = 0;
            if loopFindRedundancies == 1
                deltaTs = simStartTime - max(littleTs, [], 1);
                eigenVals = eigenVals(abs(eigenVals).^deltaTs < eigThreshold);
            else
                eigenVals = eigenVals(abs(eigenVals) < 1 - 1e-6);
            end
            
            [ ~, eigenIdx ] = sort(abs(eigenVals));
            eigenVals = eigenVals(eigenIdx);
            
            eigenGroups = [ 1 find(abs(diff(eigenVals)) > little)+1 length(eigenVals)+1 ];
            isIndependent = true(1, length(idx));
            for loopEig = 1:length(eigenGroups)-1
                sortedIndices = eigenGroups(loopEig):eigenGroups(loopEig+1)-1;
                oneEig = mean(eigenVals(sortedIndices));
                
                adjustedF = littleF - oneEig*eye(size(littleF));
                
                [~, r, e] = qr(adjustedF', 'vector');
                rSquareIdx = 1:min(size(r));
                diagR = [ diag(r(rSquareIdx, rSquareIdx))' zeros(1, max(0, diff(size(r)))) ];
                
                isIndependent(e(abs(diagR) <= little)) = false;
                if sum(~isIndependent) > 0
                    break
                end
            end
            removeMoreEigs = (loopEig < length(eigenGroups)-1);
            
            independentVars = sort(idx(isIndependent));
            dependentVars = sort(idx(~isIndependent));
            if ~isempty(dependentVars)
                break
            end
            
            if loopFindRedundancies == 1
                idx = 1:length(M.ts);
            end
        end
        
        toSolve = M.fs - oneEig*eye(size(M.fs));
        
        allCoefs = zeros(length(dependentVars), size(M.fs, 1));
        allCoefs(:, dependentVars) = -eye(length(dependentVars));
        if ~isempty(dependentVars)
            allCoefs(:, independentVars) = ((toSolve(independentVars, :)') \ (toSolve(dependentVars, :)' ))';
        end
if ~isempty(allCoefs)
    %disp(oneEig)
    if imag(oneEig) ~= 0
    end
end
        
        areComplex = (sum(imag(allCoefs) ~= 0, 2) > 0);
        allCoefs = [ allCoefs(~areComplex, :); real(allCoefs(areComplex, :)); imag(allCoefs(areComplex, :)) ];
        
        if loopFindRedundancies == 2 && ~isempty(dependentVars)
            eigenVecTimes = ones(length(dependentVars)+sum(areComplex), 1) * (M.ts(idx)');
            eigenVecTimes(abs(allCoefs) < 1.e-6) = 0;
            minStartTime = max(eigenVecTimes, [], 2) + max(1, ceil(log(eigThreshold)/log(abs(max(oneEig, 1.e-16)))));
            allCoefsToKeep = (minStartTime == min(minStartTime)) | (minStartTime <= simStartTime);
            allCoefs = allCoefs(allCoefsToKeep, :);
            if max(minStartTime(allCoefsToKeep)) > simStartTime
                simStartTime = max(minStartTime(allCoefsToKeep));
            end
        end
        
        
            % postpone constraints that are going to take us forever --
            % maybe they'll get simpler next round
        
        cSizes = sum(abs(allCoefs) > 1.e-6, 2);
        [ cSizes, idx ] = sort(cSizes);
        allCoefs = allCoefs(idx, :);
        
        if ~isempty(cSizes)
            allCoefs = allCoefs(cSizes <= 10*cSizes(1), :);
        end
        
        
            % add the constraints here
        
        printProgress('solving constraints');
        xsToKeep = true(1, size(M.fs, 2));
        
        startingVarTs = M.ts;       % use the original 'ts' so that constraint #1 doesn't update ts(1) and prevent constraint #2 from working 
        for loopDependency = 1:size(allCoefs, 1)
            printProgress([ num2str(loopDependency) ' / ' num2str(size(allCoefs, 1)) ' new constraints']);
            
            nonzeroCoefs = find(abs(allCoefs(loopDependency, :)) > 1.e-2);
            maxVarTs = max(startingVarTs(nonzeroCoefs));
            if maxVarTs < simStartTime
dbgWhichXs = ~M.xs(nonzeroCoefs, 3);
dbg12 = allCoefs(loopDependency, nonzeroCoefs)';
if abs(sum(dbg12(dbgWhichXs))) > 1.e-3
end
                addConstraints({ allCoefs(loopDependency, nonzeroCoefs)', M.xs(nonzeroCoefs, :) }, maxVarTs+1);
dbgWhichXs = ~M.xs(:, 3);
if ~isempty(M.fs)
if sum(sum(abs(M.fs*dbgWhichXs - dbgWhichXs(1:size(M.fs, 1))) > 1.e-2)) > 0
end
end
            end
        end
        
        stats(end, 2) = stats(end, 2) + size(allCoefs, 1);
        
        M.fs(:, size(M.fs, 2)+1:size(M.xs, 1)) = 0;
        
        
            % put the easiest constraints at the beginning
        
        cSizes = zeros(1, length(M.cs));
        for loopC = 1:length(M.cs)
            cSizes(loopC) = size(M.cs{loopC}{1}, 1);
        end
        
        [~, idx] = sort(cSizes);
        M.cxs = M.cxs(idx, :);
        M.cs = M.cs(:, idx);
        M.cts = M.cts(idx, :);
        
        
            % remove unused variables
        
        varsToDiscard = find(~xsToKeep);
        if (diff(size(M.fs)) == 0) && isempty(varsToDiscard)
            finished = true;
            break
        end
        
        while ~isempty(varsToDiscard)
            fsToDiscard = varsToDiscard(varsToDiscard <= size(M.fs, 1));
            M.xs(varsToDiscard, :) = [];
            M.fs(fsToDiscard, :) = [];
            M.fs(:, varsToDiscard) = [];
            M.ts(fsToDiscard, :) = [];
            xsToKeep(:, varsToDiscard) = [];
            isQueryVar(:, varsToDiscard) = [];
            varsToDiscard = find(sum(abs(M.fs) > 1.e-6, 1) == 0 & ~isQueryVar);
        end
    end
    
    stats(end, 3) = size(M.xs, 1);
    stats(end, 4) = length(M.cs);
    
    if finished, break, end
    
    
        % PART 2
        % add new variables if we've found all of the degeneracies
    
        % add/prioritize low-index variables (if not factored by
        % constraints)
    
    numFs = size(M.fs, 1);
    if numFs > maxFinalVars && diff(size(M.fs)) ~= 0
        maxFactorIndices = min(sum(M.xs(numFs+1:end, :), 2));
        for numFactorIndices = 2:maxFactorIndices
            factorXs = false(0, numBs);
            factorsToKeep = false(0);
            for loopX = find(sum(M.xs(numFs+1:end, :), 2)' >= numFactorIndices) + numFs
                newFactorPositions = nchoosek(find(M.xs(loopX, :)), numFactorIndices-1);
                factorXs = false(size(newFactorPositions, 1), numBs);
                factorsToKeep = false(size(newFactorPositions, 1), 1);
                for loopFactor = 1:size(newFactorPositions, 1)
                    factorXs(loopFactor, newFactorPositions(loopFactor, :)) = true;
                    factorsToKeep(loopFactor) = (sum(sum(~repmat(factorXs(loopFactor, :), size(M.cxs, 1), 1) & M.cxs, 2) == 0) == 0);
                end
            end
            
            uniqueXs = unique([ M.xs; factorXs(factorsToKeep, :) ], 'rows', 'stable');   % faster if we pre-allocate
            if size(uniqueXs, 1) ~= size(M.xs, 1)
                M.xs = uniqueXs;
                break;
            end
        end
        
        M.fs = [ M.fs  zeros(numFs, size(M.xs, 1)-size(M.fs, 2)) ];
        
        numXIndices = sum(M.xs, 2)';
        numXIndices(1:numFs) = -1;
        [~, newXsOrder] = sort(numXIndices);
        
        M.xs = M.xs(newXsOrder, :);             %reorder the xs
        M.fs = M.fs(:, newXsOrder);
        
        numFsToAdd = sum( sum(M.xs(size(M.fs, 1)+1:end, :), 2) == maxFactorIndices );
    else
        numFsToAdd = min(diff(size(M.fs)), maxVars-numFs);
    end
    
    xsToKeep(end+1:size(M.xs, 1)) = true;
    isQueryVar(end+1:size(M.xs, 1)) = false;
    
    
        % loop over each new variable we need to add an evolution equation
        % for
    
    if diff(size(M.fs)) ~= 0 && size(M.fs, 1) == maxVars
        disp('stats:')
        disp(stats)
        error('variable limit reached')
    end
    
    newFs = zeros(numFsToAdd, size(M.xs, 1));
    newTs = ones(numFsToAdd, 1);
    
    for loopX = 1:numFsToAdd
        
        printProgress([ num2str(loopX) ' / ' num2str(numFsToAdd) ' new features']);
        
        
            % xf identifies the variable as a product of state variables
            % (i.e. xf = [ 2 4 5 ] is the variable x_245)
        
        newX = M.xs(size(M.fs, 1)+loopX, :);
        if ~isContinuous
            xfPoly = { 1, false(1, numBs) };
            for loopBaseF = find(newX)
                xfPoly = multiplyPoly(xfPoly, M.baseFs{loopBaseF});
            end
        else
            xfPoly = { zeros(0, 1), false(0, numBs) };
            for loopBaseF = find(newX)
                otherVars = { 1, newX };
                otherVars{2}(loopBaseF) = false;
                xfPoly = addPoly(xfPoly, multiplyPoly(otherVars, M.baseFs{loopBaseF}));
            end
        end
        
        if ~isempty(M.cs)
            [ xfPolyConstrained, polyStartTime ] = constrainPoly(xfPoly, 1);
        else
            xfPolyConstrained = xfPoly;
            polyStartTime = 1;
        end
        
        theVars = findVars(xfPolyConstrained);
        newFs(loopX, theVars) = xfPolyConstrained{1};
        newTs(loopX) = polyStartTime;
    end
    stats(end, 1) = numFsToAdd;
    
    M.fs = [ M.fs  zeros(size(M.fs, 1), size(M.xs, 1)-size(M.fs, 2)); newFs ];
    M.ts = [ M.ts; newTs ];
end

printProgress(sprintf('\b'));



        % addConstraint() adds a constraint to the global list,
        % and calls itself recursively to enforce the further constraints
        % that any variable times the first constraint must give 0 or 1.
    
    function addConstraints(polyToConstrain, constraintStartTime)
        
        for loopConstraintVar = 1:size(polyToConstrain{1}, 1)
            
            idxs = ((1:size(polyToConstrain{1}, 1)) ~= loopConstraintVar);
            constraintPoly = constrainPoly( {    -polyToConstrain{1}(idxs, :) / polyToConstrain{1}(loopConstraintVar, :), ...
                                polyToConstrain{2}(idxs, :) | repmat(polyToConstrain{2}(loopConstraintVar, :), sum(idxs), 1) }, 1 );
            constraintPolyCompare = constrainPoly({ 1, polyToConstrain{2}(loopConstraintVar, :) }, 1);
            
            if ~(polysAreEqual(constraintPoly, { [ 1 1 ], polyToConstrain{2}(loopConstraintVar, :) }) ...
                    || polysAreEqual(constraintPoly, constraintPolyCompare))
                
                
                    % get rid of any pre-existing constraints that are now
                    % redundant
                
                factoredConstraints = find(sum(~M.cxs(:, polyToConstrain{2}(loopConstraintVar, :)), 2)' == 0);
                toDelete = false(1, length(factoredConstraints));
                for factorCounter = 1:length(factoredConstraints)
                    loopFactoredConstraint = factoredConstraints(factorCounter);
                    constraintCompare = multiplyPoly({ 1, M.cxs(loopFactoredConstraint, :) }, constraintPoly);
                    if size(constraintCompare{1}, 1) == size(M.cs{loopFactoredConstraint}{1}, 1)
                        if sum(abs(constraintCompare{1} - M.cs{loopFactoredConstraint}{1}) > 1.e-6) == 0 ...
                                    && sum(sum(constraintCompare{2} ~= M.cs{loopFactoredConstraint}{2})) == 0
                            toDelete(factorCounter) = true;
                        end
                    end
                end
                M.cxs(factoredConstraints(toDelete), :) = [];
                M.cs(factoredConstraints(toDelete)) = [];
                M.cts(factoredConstraints(toDelete), :) = [];
                
                
                    % actually add the new constraint here
                
                priorInstance = find(sum(M.cxs ~= repmat(polyToConstrain{2}(loopConstraintVar, :), size(M.cxs, 1), 1), 2) == 0, 1);
                if isempty(priorInstance)
                    M.cxs(end+1, :) = polyToConstrain{2}(loopConstraintVar, :);
                    M.cs{end+1} = constraintPoly;
                    M.cts(end+1, :) = constraintStartTime;
                else
                    M.cs{priorInstance} = multiplyPoly(M.cs{priorInstance}, constraintPoly);
                    M.cts(priorInstance) = constraintStartTime;
                end
                
                
                    % Find and re-constrain all xs factored by the new constraint,
                    % adding any new variables that appear.  For example,
                    % if our constraint is that x_{12} = x_3 + x_4, then
                    % multiply x_{123} by x_3 + x_4 and by x_{12} - x_4.
                
                factoredVars = find(sum(~M.xs(:, polyToConstrain{2}(loopConstraintVar, :)), 2)' == 0 & xsToKeep);
                for loopFactoredVar = factoredVars
                    if sum(sum(constraintPoly{2} & ~repmat(M.xs(loopFactoredVar, :), size(constraintPoly{2}, 1), 1), 2) == 0) == 0
                        newPoly = constrainPoly({ 1, M.xs(loopFactoredVar, :) }, constraintStartTime);
                        newPolyXs = findVars(newPoly);
                        xsToKeep(loopFactoredVar) = false;
                        isQueryVar(newPolyXs) = true;
                        
                        try     % in case the new fs are being built and M.fs is not updated with zeros yet
                            changedFs = (abs(M.fs(:, loopFactoredVar)) > 1.e-6);
                            M.ts(changedFs) = max(M.ts(changedFs), constraintStartTime);
                            toAdd = M.fs(:, loopFactoredVar)*(newPoly{1}(:, 1)');
                            M.fs(:, loopFactoredVar) = 0;
                        end
                        if max(newPolyXs) > size(M.fs, 2)
                            M.fs(:, max(newPolyXs)) = 0;
                        end
                        M.fs(:, newPolyXs) = M.fs(:, newPolyXs) + toAdd;
                    end
                end
            end
        end
    end



        % constrainPoly() applies the existing constraints to a given
        % polynomial by multiplying each term of the polynomial by its
        % 'factors'
    
    function [ constrainedPoly, startTime ] = constrainPoly(unconstrainedTerms, initialStartTime)
        startTime = initialStartTime;
        constrainedPoly = { zeros(0, 1), false(0, numBs) };
        
        while ~isempty(unconstrainedTerms{1})
            
            madeAChange = false;
            for loopConstraint = 1:length(M.cs)
                newTerms = { zeros(0, 1), false(0, numBs) };
                
                factoredTerms = find(sum(~unconstrainedTerms{2}(:, M.cxs(loopConstraint, :)), 2) == 0)';
                changedTerm = false(1, length(unconstrainedTerms{1}));
                for loopFactoredTerm = factoredTerms
                    if sum(sum(M.cs{loopConstraint}{2}(:, ~unconstrainedTerms{2}(loopFactoredTerm, :)), 2) == 0) == 0
                        originalTerm = { unconstrainedTerms{1}(loopFactoredTerm, :), unconstrainedTerms{2}(loopFactoredTerm, :) };
                        oneConstrainedTerm = multiplyPoly(originalTerm, M.cs{loopConstraint});
                        changedTerm(loopFactoredTerm) = true;
                        madeAChange = true;
                        startTime = max([ startTime; M.cts(loopConstraint) ]);
                        newTerms = { [ newTerms{1}; oneConstrainedTerm{1} ], ...
                                     [ newTerms{2}; oneConstrainedTerm{2} ] };
                    end
                end
                
                [ unconstrainedTerms, addedAConstraint ] = checkOversizeConstraints(reducePoly( ...
                    { [ unconstrainedTerms{1}(~changedTerm); newTerms{1} ], [ unconstrainedTerms{2}(~changedTerm, :); newTerms{2} ] } ));
                if addedAConstraint
                    madeAChange = true;
                    break;
                end
            end
            if ~madeAChange
                break
            end
        end
        
        constrainedPoly = reducePoly(unconstrainedTerms);
    end


    function [ polyCoefs, addedAConstraint ] = checkOversizeConstraints(polyCoefs)
        addedAConstraint = false;
        lookAgain = true;
        while lookAgain
            lookAgain = false;
            
            realCoefs = polyCoefs{1}';
            absCoefs = abs(realCoefs);
            sortedCoefs = sort(absCoefs);
            coefJumps = [ 0 find(abs(diff(sortedCoefs))-1 > 1.e-6) length(sortedCoefs) ] + 1;
            
            for loopCoef = 1:length(coefJumps)-1
                meanCoef = mean(sortedCoefs(coefJumps(loopCoef):coefJumps(loopCoef+1)-1));
                coefMultipliers = round(abs(realCoefs/meanCoef)) .* sign(realCoefs);
                residualCoefs = (realCoefs - meanCoef*coefMultipliers)';
                if meanCoef - sum(abs(residualCoefs)) - 1. > 1.e-2
                    coefsToPullOut = (coefMultipliers ~= 0);
                    addedAConstraint = true;
                    addConstraints({ coefMultipliers(coefsToPullOut)', polyCoefs{2}(coefsToPullOut, :) }, simStartTime);
                    polyCoefs{1}(:, 1) = residualCoefs;
                    polyCoefs = reducePoly(polyCoefs);
                    lookAgain = ~isempty(polyCoefs{1});
                    break;
                end
            end
        end
    end
    
    
        % findVars() finds all correlation variables in a polynomial and
        % returns their position in the full list [ xs, xfs, newXs ].
        % Any new variables get added to newXs[].
    
    function xIDs = findVars(thePoly)
        
        xIDs = zeros(size(thePoly{2}, 1), 1);
        for loopPolyTerm = 1:size(thePoly{2}, 1);
            oneID = find((sum(M.xs == repmat(thePoly{2}(loopPolyTerm, :), size(M.xs, 1), 1), 2) == numBs)' ...
                        & xsToKeep, 1, 'first');
            if isempty(oneID)
                oneID = size(M.xs, 1)+1;
                M.xs(size(M.xs, 1)+1, :) = thePoly{2}(loopPolyTerm, :);
                xsToKeep(1, size(xsToKeep, 2)+1) = true;
                isQueryVar(1, size(isQueryVar, 2)+1) = false;
            end
            xIDs(loopPolyTerm) = oneID;
        end
    end


    function ifEq = polysAreEqual(poly1, poly2)
        ifEq = false;
        if size(poly1{1}, 1) == size(poly2{1}, 1)
            if sum(abs(poly1{1} - poly2{1}) > 1.e-6) == 0 ...
                    && sum(sum(  poly1{2} ~= poly2{2}  )) == 0
                ifEq = true;
            end
        end
    end
    
    
    function printProgress(progressString)
        %return
        disp([backspaceString progressString])
        backspaceString = sprintf(repmat('\b', 1, length(progressString)+1));
    end
end