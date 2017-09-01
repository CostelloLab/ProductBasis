% buildF() approximates the transition rules of a Boolean
% network.  Each variable 'x' can be one of the Boolean variables 'b' or a
% product of them:  for example x13 is b1*b3.  A transition function is
% denoted 'f'.

function [ M, simStartTime, stats ] = ...
    buildF(M, maxVars, maxFinalVars, isContinuous, eigThreshold, allowFactorization)


    % our initial features are just the fundamental Boolean state variables

backspaceString = '';
finished = false;
realTol = 1.e-6;
maxAllowedCoef = realTol*flintmax();

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

if ~exist('allowFactorization', 'var')
    allowFactorization = true;
end

xsToKeep = true(1, size(M.xs, 1));

isPBN = false;
oversizeCoefThreshold = 1;
for loopB = 1:numBs
    if sum(M.baseFs{loopB}{1} ~= round(M.baseFs{loopB}{1} )) > 0
        isPBN = true;
        allowFactorization = false;
        oversizeCoefThreshold = 1 / eigThreshold;
        break;
    end
end


    % don't allow equation-reduction for PBNs or continuous-time networks
    % (it doesn't really work..)
    
% if isPBN || isContinuous
%     maxFinalVars = maxVars;
% end


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
        
        printProgress('finding decaying modes..');
        
        little = max(size(M.fs)) * eps(norm(M.fs)) + norm(M.ferr);
        
        if isPBN
            squareF = M.fs(:, 1:size(M.fs, 1));
            [ eigenVecs, eigenValsMat, eigenErrs ] = condeig(squareF');
            eigenVals = diag(eigenValsMat)';
            eigenErrs = eigenErrs*norm(M.ferr(:, 1:size(M.fs, 1)) + eps(squareF));
            
            [ ~, eigenIdx ] = sort(abs(eigenVals));
            eigenVals = eigenVals(eigenIdx);
            eigenVecs = eigenVecs(:, eigenIdx);
            eigenErrs = eigenErrs(eigenIdx);
            eigenVecErr = eigenVecs'*M.fs(:, size(M.fs, 1)+1:end);
            
            eigenGroups = [ 1 find(abs(diff(eigenVals, 1, 2)) > little)+1 length(eigenVals)+1 ];
            removableGroups = ((diff(eigenGroups) == 1) & (sqrt(mean(eigenVecErr(eigenGroups(1:end-1), :).^2, 2)) > little)');
            eigenVals(eigenGroups(removableGroups)) = [];
            eigenVecs(:, eigenGroups(removableGroups)) = [];
            eigenErrs(eigenGroups(removableGroups)) = [];
            
            eigenVecs = eigenVecs(:, abs(eigenVals) < 1 - little);
            eigenErrs = eigenErrs(abs(eigenVals) < 1 - little);
            eigenVals = eigenVals(abs(eigenVals) < 1 - little);
            eigenGroups = [ 1 find(abs(diff(eigenVals)) > little)+1 length(eigenVals)+1 ];
        else
            eigenVals = 0;
            eigenErrs = 0;
            eigenGroups = [ 1 2 ];
        end
        
        isIndependent = true(1, size(M.fs, 1));
        numDependentVars = 0;
        if ~isempty(eigenVals)
            for loopEig = 1:length(eigenGroups)-1
                sortedIndices = eigenGroups(loopEig):eigenGroups(loopEig+1)-1;
                oneEig = mean(eigenVals(sortedIndices));
                eigErr = max(eigenErrs(sortedIndices));
%                 if isPBN
%                     try
%                         [ polishedEigs, ~, polishedEigErrs ] = polishEigs(squareF', oneEig, zeros(0, length(sortedIndices)));
%                         oneEig = abs(mean(diag(polishedEigs)));
%                         eigErr = sum(diag(polishedEigErrs));
%                     catch
%                         eigErr = Inf;
%                     end
%                 end
                
                adjustedF = M.fs - oneEig*eye(size(M.fs));

                [~, r, e] = qr(adjustedF', 'vector');
                rSquareIdx = 1:min(size(r));
                diagR = [ diag(r(rSquareIdx, rSquareIdx))' zeros(1, max(0, diff(size(r)))) ];

                isIndependent(e(abs(diagR) <= little)) = false;
                numDependentVars = sum(~isIndependent);
                if numDependentVars > 0
                    break
                end
            end
            removeMoreEigs = (loopEig < length(eigenGroups)-1);
        else
            removeMoreEigs = false;
        end
        
        addedAConstraint = false;
        if numDependentVars > 0
            
            toSolve = M.fs - oneEig*eye(size(M.fs));
            toSolveErr = M.ferr + eigErr*eye(size(M.ferr));
            
            
                % Calculate the coefficients of the dependency polynomials by
                % matrix division (dependent / independent variables).
            
            allCoefs = zeros(numDependentVars, size(M.fs, 1));
            allCoefErrs = allCoefs;
            allCoefErrs2 = allCoefs;
            
            allCoefs(:, ~isIndependent) = -eye(numDependentVars);
            allCoefs(:, isIndependent) = toSolve(~isIndependent, :) / toSolve(isIndependent, :);
            
%             minSV = min(svd(toSolve(isIndependent, :)));
%             allCoefErrs(:, isIndependent) = repmat( ...
%                 rowVecNorm(toSolve(~isIndependent, :)) * norm(toSolveErr(isIndependent, :) + eps(toSolve(isIndependent, :))) / minSV^2 ...
%                     + rowVecNorm(toSolveErr(~isIndependent, :) + eps(toSolve(~isIndependent, :))) / minSV, ...
%                 1, sum(isIndependent) );
%             indepInv = eye(size(M.fs, 2)) / toSolve(isIndependent, :);
%             allCoefErrs2(:, isIndependent) = toSolveErr(~isIndependent, :) * abs(indepInv) ...
%                 + abs( abs(allCoefs(:, isIndependent))*toSolveErr(isIndependent, :) / abs(toSolve(isIndependent, :)) );
%             allCoefErrs = allCoefErrs2;
            independentSVs = svd(toSolve(isIndependent, :));
            if isempty(independentSVs)
                minSV = 1;
                maxSV = 1;
            else
                minSV = min(independentSVs);
                maxSV = max(independentSVs);
            end
            allCoefErrs = repmat( sqrt(sum((eps(toSolve(~isIndependent, :))).^2, 2))*maxSV/minSV, 1, size(M.fs, 1)) ...
                + norm(toSolveErr(~isIndependent, :))*minSV ...
                + norm(toSolve(isIndependent, :))*norm(toSolveErr(~isIndependent, :))*minSV/maxSV;
            
            areComplex = (sum(abs(imag(allCoefs)) > allCoefErrs, 2) > 0);
            allCoefs = [ allCoefs(~areComplex, :); real(allCoefs(areComplex, :)); imag(allCoefs(areComplex, :)) ];
            allCoefErrs = [ allCoefErrs(~areComplex, :); allCoefErrs(areComplex, :); allCoefErrs(areComplex, :) ];
            
            eigenVecTimes = ones(numDependentVars+sum(areComplex), 1) * (M.ts');
            eigenVecTimes(abs(allCoefs) <= allCoefErrs) = 0;
            constraintTimeConst = log(max(abs(oneEig), 1.e-16));
            minStartTime = max(max(eigenVecTimes, [], 2), 1) + log(eigThreshold)/constraintTimeConst;
            allCoefsToKeep = (minStartTime == min(minStartTime)) | (minStartTime <= max([ M.ts; M.cts ]));
            allCoefs = allCoefs(allCoefsToKeep, :);
            allCoefErrs = allCoefErrs(allCoefsToKeep, :);
            minStartTime = minStartTime(allCoefsToKeep);
            
            
                % postpone constraints that are going to take us forever --
                % maybe they'll get simpler next round
            
            cSizes = sum(abs(allCoefs) > allCoefErrs, 2);
            [ cSizes, idx ] = sort(cSizes);
            allCoefs = allCoefs(idx, :);
            allCoefErrs = allCoefErrs(idx, :);
            
            if ~isempty(cSizes)
                minStartTime = minStartTime(idx);
                allCoefs = allCoefs(cSizes <= 2*cSizes(1), :);
                allCoefErrs = allCoefErrs(cSizes <= 2*cSizes(1), :);
            end
            
            
                % add the constraints here
            
            printProgress('solving constraints');
            
            for loopDependency = 1:size(allCoefs, 1)
                printProgress([ num2str(loopDependency) ' / ' num2str(size(allCoefs, 1)) ' new constraints']);
                
                nonzeroCoefs = find(abs(allCoefs(loopDependency, :)) > allCoefErrs(loopDependency, :));
                addConstraints({ allCoefs(loopDependency, nonzeroCoefs)', ...
                        allCoefErrs(loopDependency, nonzeroCoefs)', M.xs(nonzeroCoefs, :) }, minStartTime(loopDependency), constraintTimeConst);
            end
            
            stats(end, 2) = stats(end, 2) + size(allCoefs, 1);
        end
        
        M.fs(:, size(M.fs, 2)+1:size(M.xs, 1)) = 0;
        M.ferr(:, size(M.ferr, 2)+1:size(M.xs, 1)) = 0;
        
        if removeMoreEigs && ~addedConstraint
            printProgress(sprintf('\b'));
            error('buildF:stuck', 'buildF:  got stuck')
        end
        
        
            % put the easiest constraints at the beginning
        
        cSizes = zeros(1, length(M.cs));
        for loopC = 1:length(M.cs)
            cSizes(loopC) = size(M.cs{loopC}{1}, 1);
        end
        
        [~, idx] = sort(cSizes);
        M.cxs = M.cxs(idx, :);
        M.cs = M.cs(:, idx);
        M.cts = M.cts(idx, :);
        
        if (diff(size(M.fs)) == 0) && sum(~xsToKeep) == 0
            finished = true;
            break
        end
        removeUnusedVars()
    end
    
    stats(end, 3) = size(M.xs, 1);
    stats(end, 4) = length(M.cs);
    
    if finished, break, end
    
    
        % PART 2
        % add new variables if we've found all of the degeneracies
    
        % add/prioritize low-index variables (if not factored by
        % constraints)
    
    numFs = size(M.fs, 1);
    if numFs > maxFinalVars && diff(size(M.fs)) ~= 0 && allowFactorization
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
        M.ferr = [ M.ferr  zeros(numFs, size(M.xs, 1)-size(M.ferr, 2)) ];
        
        numXIndices = sum(M.xs, 2)';
        numXIndices(1:numFs) = -1;
        [~, newXsOrder] = sort(numXIndices);
        
        M.xs = M.xs(newXsOrder, :);             %reorder the xs
        M.fs = M.fs(:, newXsOrder);
        M.ferr = M.ferr(:, newXsOrder);
        
        numFsToAdd = sum( sum(M.xs(size(M.fs, 1)+1:end, :), 2) == maxFactorIndices );
    else
        numFsToAdd = min(diff(size(M.fs)), maxVars-numFs);
    end
    
    xsToKeep(end+1:size(M.xs, 1)) = true;
    
    
        % loop over each new variable we need to add an evolution equation
        % for
    
    if diff(size(M.fs)) ~= 0 && size(M.fs, 1) == maxVars
        disp('stats:')
        disp(stats)
        error('buildF:varLimit', 'variable limit reached')
    end
    
    newFs = zeros(numFsToAdd, size(M.xs, 1));
    newFerr = newFs;
    newTs = ones(numFsToAdd, 1);
    
    for loopX = 1:numFsToAdd
        
        printProgress([ num2str(loopX) ' / ' num2str(numFsToAdd) ' new features']);
        
        
            % xf identifies the variable as a product of state variables
            % (i.e. xf = [ 2 4 5 ] is the variable x_245)
        
        newX = M.xs(size(M.fs, 1)+loopX, :);
        if ~isContinuous
            xfPoly = { 1, 0, false(1, numBs) };
            for loopBaseF = find(newX)
                xfPoly = multiplyPoly(xfPoly, M.baseFs{loopBaseF});
            end
        else
            xfPoly = { zeros(0, 1), zeros(0, 1), false(0, numBs) };
            for loopBaseF = find(newX)
                otherVars = { 1, 0, newX };
                otherVars{3}(loopBaseF) = false;
                xfPoly = addPoly(xfPoly, multiplyPoly(otherVars, M.baseFs{loopBaseF}));
            end
        end
        
        if ~isempty(M.cs)
            [ xfPolyConstrained, polyStartTime ] = constrainPoly(xfPoly, 1, false);
        else
            xfPolyConstrained = xfPoly;
            polyStartTime = 1;
        end
        
        theVars = findVars(xfPolyConstrained);
        newFs(loopX, theVars) = xfPolyConstrained{1};
        newFerr(loopX, theVars) = xfPolyConstrained{2};
        newTs(loopX) = polyStartTime;
    end
    stats(end, 1) = numFsToAdd;
    
    M.fs = [ M.fs  zeros(size(M.fs, 1), size(M.xs, 1)-size(M.fs, 2)); newFs ];
    M.ferr = [ M.ferr  zeros(size(M.ferr, 1), size(M.xs, 1)-size(M.ferr, 2)); newFerr ];
    M.ts = [ M.ts; newTs ];
    
    removeUnusedVars();         % these can be caused by constrainPoly()
end

simStartTime = max([ M.ts; M.cts ]);
printProgress(sprintf('\b'));



        % remove variables marked in xsToKeep for discarding
    
    function removeUnusedVars()
        
        varsToDiscard = find(~xsToKeep);
        fsToDiscard = varsToDiscard(varsToDiscard <= size(M.fs, 1));
        
        M.xs(varsToDiscard, :) = [];
        M.fs(fsToDiscard, :) = [];
        M.fs(:, varsToDiscard) = [];
        M.ferr(fsToDiscard, :) = [];
        M.ferr(:, varsToDiscard) = [];
        M.ts(fsToDiscard, :) = [];
        xsToKeep(:, varsToDiscard) = [];
    end



        % addConstraint() adds a constraint to the global list,
        % and calls itself recursively to enforce the further constraints
        % that any variable times the first constraint must give 0 or 1.
    
    function addConstraints(polyToConstrain, modeEvaporationTime, constraintTimeConstant)
        
        if length(polyToConstrain{1}) > 1
            scaledPoly = polyToConstrain;
            scaledPoly{1} = scaledPoly{1}/min(abs(polyToConstrain{1}));
            [ ~, constraintWasAdded ] = checkOversizeConstraints(scaledPoly, true, false, true);
            if constraintWasAdded
                return
            end
        end
        
        polyToConstrain = constrainPoly(polyToConstrain, ceil(modeEvaporationTime), true);
        if ~allowFactorization && ~isempty(polyToConstrain{1})
            [ ~, varOrder ] = sortrows(polyToConstrain{3});
            termsToConsider = varOrder(1);
        else
            termsToConsider = 1:size(polyToConstrain{1}, 1);
        end
        
        for loopConstraintVar = termsToConsider
            if sum(sum(polyToConstrain{3}(:, ~polyToConstrain{3}(loopConstraintVar, :)), 2) == 0) == 1
                
                idxs = ((1:size(polyToConstrain{1}, 1)) ~= loopConstraintVar);
                if ~allowFactorization
                    newPolyIndices = polyToConstrain{3}(idxs, :);
                else
                    newPolyIndices = polyToConstrain{3}(idxs, :) | repmat(polyToConstrain{3}(loopConstraintVar, :), sum(idxs), 1);
                end
                constraintPoly = constrainPoly( {    -polyToConstrain{1}(idxs, :) / polyToConstrain{1}(loopConstraintVar, :), ...
                                    polyToConstrain{2}(idxs, :) / abs(polyToConstrain{1}(loopConstraintVar, :)) + ...
                                        abs(polyToConstrain{1}(idxs, :)) * polyToConstrain{2}(loopConstraintVar, :) ...
                                                    / polyToConstrain{1}(loopConstraintVar, :).^2, ...
                                    newPolyIndices }, 1, true );
                constraintPolyCompare = constrainPoly({ 1, 0, polyToConstrain{3}(loopConstraintVar, :) }, 1, true);
                
                if ~(polysAreEqual(constraintPoly, { [ 1 1 ], [ 0 0 ], polyToConstrain{3}(loopConstraintVar, :) }) ...
                        || polysAreEqual(constraintPoly, constraintPolyCompare))
                
                
                        % Get rid of any pre-existing constraints that are
                        % now redundant (for deterministic networks only;
                        % for PBNs old constraints can be useful if they
                        % involve fewer indices than the new one, because
                        % they're not allowed to multiply same indices)
                    
                    if allowFactorization
                        factoredConstraints = find(sum(~M.cxs(:, polyToConstrain{3}(loopConstraintVar, :)), 2)' == 0);
                        toDelete = false(1, length(factoredConstraints));
                        for factorCounter = 1:length(factoredConstraints)
                            loopFactoredConstraint = factoredConstraints(factorCounter);
                            constraintCompare = multiplyPoly({ 1, 0, M.cxs(loopFactoredConstraint, :) }, constraintPoly);
                             if polysAreEqual(constraintCompare, M.cs{loopFactoredConstraint})
                                toDelete(factorCounter) = true;
                            end
                        end
                        M.cxs(factoredConstraints(toDelete), :) = [];
                        M.cs(factoredConstraints(toDelete)) = [];
                        M.cts(factoredConstraints(toDelete), :) = [];
                    end
                    
                    
                        % actually add the new constraint here
                    
                    constraintStartTime = ceil(modeEvaporationTime + max(0, log(abs(polyToConstrain{1}(loopConstraintVar, :)))/constraintTimeConstant));
                    priorInstance = find(sum(M.cxs ~= repmat(polyToConstrain{3}(loopConstraintVar, :), size(M.cxs, 1), 1), 2) == 0, 1);
                    if isempty(priorInstance) || ~allowFactorization
                        M.cxs(end+1, :) = polyToConstrain{3}(loopConstraintVar, :);
                        M.cs{end+1} = constraintPoly;
                        M.cts(end+1, :) = constraintStartTime;
                    else
                        constraintPoly = checkOversizeConstraints(multiplyPoly(M.cs{priorInstance}, constraintPoly), true, true, false);
                        M.cs{priorInstance} = constraintPoly;
                        M.cts(priorInstance) = constraintStartTime;
                    end
                    addedConstraint = true;
                    
                    
                        % Find and re-constrain all xs factored by the new constraint,
                        % adding any new variables that appear.  For example,
                        % if our constraint is that x_{12} = x_3 + x_4, then
                        % multiply x_{123} by x_3 + x_4 and by x_{12} - x_4.
                    
                    if ~allowFactorization
                        factoredVars = find(sum(repmat(polyToConstrain{3}(loopConstraintVar, :), size(M.xs, 1), 1) ~= M.xs, 2)' == 0 & xsToKeep);
                    else
                        factoredVars = find(sum(~M.xs(:, polyToConstrain{3}(loopConstraintVar, :)), 2)' == 0 & xsToKeep);
                    end
                    for loopFactoredVar = factoredVars
                        if sum(sum(constraintPoly{3} & ~repmat(M.xs(loopFactoredVar, :), size(constraintPoly{3}, 1), 1), 2) == 0) == 0
                            [ newPoly, eqStartTime ] = constrainPoly({ 1, 0, M.xs(loopFactoredVar, :) }, constraintStartTime, true);
                            newPolyXs = findVars(newPoly);
                            xsToKeep(loopFactoredVar) = false;
                            
                            try     % in case the new fs are being built and M.fs is not updated with zeros yet
                                changedFs = (abs(M.fs(:, loopFactoredVar)) > M.ferr(:, loopFactoredVar));
                                M.ts(changedFs) = max(M.ts(changedFs), eqStartTime);
                                toAdd = M.fs(:, loopFactoredVar)*(newPoly{1}(:, 1)');
                                toAddErr = M.ferr(:, loopFactoredVar)*abs(newPoly{1}(:, 1)') + abs(M.fs(:, loopFactoredVar))*(newPoly{2}(:, 1)');
                                M.fs(:, loopFactoredVar) = 0;
                            catch
                                toAdd = 0;
                                toAddErr = 0;
                            end
                            if max(newPolyXs) > size(M.fs, 2)
                                M.fs(:, max(newPolyXs)) = 0;
                                M.ferr(:, max(newPolyXs)) = 0;
                            end
                            M.fs(:, newPolyXs) = M.fs(:, newPolyXs) + toAdd;
                            M.ferr(:, newPolyXs) = M.ferr(:, newPolyXs) + toAddErr;
                        end
                    end
                end
            end
        end
    end



        % constrainPoly() applies the existing constraints to a given
        % polynomial by multiplying each term of the polynomial by its
        % 'factors'
    
    function [ constrainedPoly, startTime ] = constrainPoly(unconstrainedTerms, initialStartTime, allowedToAddConstraints)
        startTime = initialStartTime;
        
        while ~isempty(unconstrainedTerms{1})
            
            madeAChange = false;
            for loopConstraint = 1:length(M.cs)
                newTerms = { zeros(0, 1), zeros(0, 1), false(0, numBs) };
                
                if ~allowFactorization
                    factoredTerms = find(sum(repmat(M.cxs(loopConstraint, :), length(unconstrainedTerms{1}), 1) ~= unconstrainedTerms{3}, 2)' == 0);
                else
                    factoredTerms = find(sum(~unconstrainedTerms{3}(:, M.cxs(loopConstraint, :)), 2) == 0)';
                end
                changedTerm = false(1, length(unconstrainedTerms{1}));
                for loopFactoredTerm = factoredTerms
                    if sum(sum(M.cs{loopConstraint}{3}(:, ~unconstrainedTerms{3}(loopFactoredTerm, :)), 2) == 0) == 0
                        originalTerm = { unconstrainedTerms{1}(loopFactoredTerm, :), ...
                            unconstrainedTerms{2}(loopFactoredTerm, :), unconstrainedTerms{3}(loopFactoredTerm, :) };
                        constraintIsAllowed = true;
                        if ~allowFactorization
                            originalTerm{3} = (originalTerm{3} & ~M.cxs(loopConstraint, :));
                            if sum(sum(M.cs{loopConstraint}{3}(:, originalTerm{3}))) ~= 0
                                constraintIsAllowed = false;
                            end
                        end
                        
                        if constraintIsAllowed
                            oneConstrainedTerm = multiplyPoly(originalTerm, M.cs{loopConstraint});
                            changedTerm(loopFactoredTerm) = true;
                            madeAChange = true;
                            startTime = max([ startTime; M.cts(loopConstraint) ]);
                            newTerms = { [ newTerms{1}; oneConstrainedTerm{1} ], ...
                                         [ newTerms{2}; oneConstrainedTerm{2} ], ...
                                         [ newTerms{3}; oneConstrainedTerm{3} ] };
                        end
                    end
                end
                
                [ unconstrainedTerms, addedAConstraint ] = checkOversizeConstraints(reducePoly( ...
                    { [ unconstrainedTerms{1}(~changedTerm); newTerms{1} ], [ unconstrainedTerms{2}(~changedTerm); newTerms{2} ], ...
                    [ unconstrainedTerms{3}(~changedTerm, :); newTerms{3} ] } ), allowedToAddConstraints, true, false);
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


    function [ polyCoefs, addedAConstraint ] = checkOversizeConstraints(polyCoefs, allowedToAddConstraints, errorOnOversize, returnOnAdd)
        
        polyStartingTime = max([ M.ts; M.cts ]);
        
        if max(abs(polyCoefs{2})) > realTol
            printProgress(sprintf('\b'));
            error('buildF:overTol', 'buildF:  loss of numerical precision')
        end
        
        addedAConstraint = false;
        lookAgain = true;
        while lookAgain
            lookAgain = false;
            
            realCoefs = polyCoefs{1}';
            absCoefs = abs(realCoefs);
            [ sortedCoefs, sortidx ] = sort(absCoefs);
            sortedErrs = polyCoefs{2}(sortidx)';
            if ~isempty(sortedCoefs) && errorOnOversize
                if abs(sortedCoefs(end)) > maxAllowedCoef
                    printProgress(sprintf('\b'));
                    error('buildF:coeftoobig', 'buildF:  coefficient too big')
                end
            end
            if length(sortedCoefs) > 1
                diffCoefs = abs(diff(sortedCoefs));
            else
                diffCoefs = zeros(1, 0);
            end
            coefJumps = [ 0 find(diffCoefs-1 > sortedErrs(1:end-1)+sortedErrs(2:end)) length(sortedCoefs) ] + 1;
            
            for loopCoef = length(coefJumps)-1:-1:1
                meanCoef = mean(sortedCoefs(coefJumps(loopCoef):coefJumps(loopCoef+1)-1));
                coefMultipliers = round(abs(realCoefs/meanCoef)) .* sign(realCoefs);
                residualCoefs = (realCoefs - meanCoef*coefMultipliers)';
                if meanCoef - sum(abs(residualCoefs)) - oversizeCoefThreshold > sum(polyCoefs{2})
                    coefsToPullOut = (coefMultipliers ~= 0);
                    if allowedToAddConstraints
                        addedAConstraint = true;
                        addConstraints({ coefMultipliers(coefsToPullOut)', polyCoefs{2}(coefsToPullOut), polyCoefs{3}(coefsToPullOut, :) }, polyStartingTime, -inf);
                    end
                    polyCoefs{1}(:, 1) = residualCoefs;
                    polyCoefs = reducePoly(polyCoefs);
                    lookAgain = (~isempty(polyCoefs{1})) && (~returnOnAdd);
                    break;
                end
            end
        end
    end
    
    
        % findVars() finds all correlation variables in a polynomial and
        % returns their position in the full list [ xs, xfs, newXs ].
        % Any new variables get added to newXs[].
    
    function xIDs = findVars(thePoly)
        
        xIDs = zeros(size(thePoly{3}, 1), 1);
        for loopPolyTerm = 1:size(thePoly{3}, 1);
            oneID = find((sum(M.xs == repmat(thePoly{3}(loopPolyTerm, :), size(M.xs, 1), 1), 2) == numBs)' ...
                        & xsToKeep, 1, 'first');
            if isempty(oneID)
                oneID = size(M.xs, 1)+1;
                M.xs(size(M.xs, 1)+1, :) = thePoly{3}(loopPolyTerm, :);
                xsToKeep(1, size(xsToKeep, 2)+1) = true;
            end
            xIDs(loopPolyTerm) = oneID;
        end
    end


    function ifEq = polysAreEqual(poly1, poly2)
        ifEq = false;
        if size(poly1{1}, 1) == size(poly2{1}, 1)
            if sum(abs(poly1{1} - poly2{1}) > (poly1{2} + poly2{2})) == 0 ...
                    && sum(sum(  poly1{3} ~= poly2{3}  )) == 0
                ifEq = true;
            end
        end
    end
    

    function vNorm = rowVecNorm(v0)
        vNorm = sqrt(sum(v0.^2, 2));
    end
    

    function printProgress(progressString)
        disp([backspaceString progressString])
        backspaceString = sprintf(repmat('\b', 1, length(progressString)+1));
    end
end