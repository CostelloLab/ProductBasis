numVars = 10;
maxInputs = 4;
numSims = 100;
numTimeSteps = 100;
doProb = true;
doContinuous = false;
doEqReduction = true;
doSaves = true;

simCounter = 1;
backspaceString = '';
errsString = '';
numVarLimitErrs = 0;
numCoefErrs = 0;
numTolErrs = 0;
numStuckErrs = 0;
numTimeErrs = 0;
varCeiling = 2^numVars;
errHistory = zeros(0, 2);
if doContinuous
    diagConstant = 0;
    subTsteps = 100;
else
    diagConstant = 1;
    subTsteps = 1;
end
%bsAll = zeros(numTimeSteps+1, numVars, numSims);

if doProb || doContinuous
    tol = 1.e-6 + 3/sqrt(numSims);
else
    tol = 1.e-6;
end

dbgCcount = 0;
while true
    
    wasErr = false;
    progressString = [ 'sim #' num2str(simCounter) errsString ];
    disp([backspaceString progressString])
    backspaceString = sprintf(repmat('\b', 1, length(progressString)+1));
    rates = rand(1, ceil(numVars*rand()));
% rates = [ .5 1 ];
    
    logicTables = cell(1, numVars);
    simLogicTables = cell(1, numVars);
    for loopVar = 1:numVars
        numInputs = ceil(maxInputs*rand());
        whichInputs = randperm(numVars, numInputs);
        if doProb || doContinuous
            whichInputs(whichInputs==loopVar) = [];
            whichInputs = [ loopVar whichInputs ];
            oneRate = rates(ceil(length(rates)*rand()));
            simStepRate = oneRate/subTsteps;
            oneTable = oneRate * round(rand(1, 2^(length(whichInputs)-1)));
            simLogicTables{loopVar} = { whichInputs, [ simStepRate*oneTable, simStepRate*oneTable + (1-simStepRate) ] };
            if doContinuous
                logicTables{loopVar} = { whichInputs, [ oneRate*oneTable, oneRate*oneTable - oneRate ] };
            else
                logicTables{loopVar} = simLogicTables{loopVar};
            end
        else
            logicTables{loopVar} = { whichInputs, round(rand(1, 2^numInputs)) };
            simLogicTables{loopVar} = logicTables{loopVar};
        end
    end
    
    varToLookAt = ceil(numVars*rand());
    readoutVector = (1:numVars)==varToLookAt;
    if doEqReduction
        varCeiling = 2^floor(numVars*rand());
    end
% logicTables = holdLogicTables;
% simLogicTables = holdSimLogicTables;
% readoutVector = holdReadoutVector;
% varCeiling = holdVarCeiling;
    M = polynomialModel(logicTables, readoutVector);
    
% holdLogicTables = logicTables; holdSimLogicTables = simLogicTables; holdReadoutVector = readoutVector;  holdVarCeiling = varCeiling;
    if doSaves
        save('~/Desktop/tmp/holdVars', 'logicTables', 'simLogicTables', 'readoutVector', 'varCeiling')
    end
    try
        [ M, tStart ] = buildF(M, 2^numVars, varCeiling, doContinuous);
    catch boolException
        if strcmp(boolException.identifier, 'buildF:varLimit')
            numVarLimitErrs = numVarLimitErrs+1;
            errHistory = [ errHistory; simCounter 1 ];
            wasErr = true;
        elseif strcmp(boolException.identifier, 'buildF:coeftoobig')
            numCoefErrs = numCoefErrs+1;
            errHistory = [ errHistory; simCounter 2 ];
            wasErr = true;
        elseif strcmp(boolException.identifier, 'buildF:overTol')
            numTolErrs = numTolErrs+1;
            errHistory = [ errHistory; simCounter 3 ];
            wasErr = true;
        elseif strcmp(boolException.identifier, 'buildF:stuck')
            numStuckErrs = numStuckErrs+1;
            errHistory = [ errHistory; simCounter 4 ];
            wasErr = true;
        else
            error('unknown error')
        end
        errsString = [ '  (' num2str(numVarLimitErrs) ' + ' num2str(numCoefErrs) ' + ' num2str(numTolErrs) ...
            ' + ' num2str(numStuckErrs) ' + ' num2str(numTimeErrs) ' errs)' ];
    end
    
    if ~isempty(M.cs)
        dbgCcount = dbgCcount + 1;
    end
    while true
        simHistory = zeros(numSims, numTimeSteps+1);
        bs = zeros(numSims, numVars);
        for loopSim = 1:numSims
            [ bsTemp, xs ] = simBoolModel(simLogicTables, round(rand(1, numVars)), readoutVector, numTimeSteps, subTsteps);
            bs(loopSim, :) = bsTemp(min(tStart, size(bsTemp, 1)), :);
            simHistory(loopSim, :) = xs;
        end
        
        allXs = zeros(numSims, size(M.xs, 1));
        for loopSim = 1:numSims
            for loopX = 1:size(M.xs, 1)
                allXs(loopSim, loopX) = prod(bs(loopSim, M.xs(loopX, :)), 2);
            end
        end
        
        [ validHistory, extraVar ] = evolveF(M, mean(allXs, 1)', readoutVector, numTimeSteps-tStart+1, doContinuous);
        if isempty(extraVar) || wasErr
            break
        end
        
        M.xs = [ M.xs; extraVar ];
        M.fs = [ M.fs, zeros(size(M.fs, 1), 1) ];
        [ M, tStart ] = buildF(M, 2^numVars, varCeiling);
    end
    if tStart > numTimeSteps+1
        predictedHistory = mean(simHistory, 1);
        if tStart > 1e5
            numTimeErrs = numTimeErrs+1;
            errsString = [ '  (' num2str(numCoefErrs) ' + ' num2str(numTolErrs) ' + ' num2str(numTimeErrs) ' errs)' ];
        end
    else
        predictedHistory = [ mean(simHistory(:, 1:tStart), 1), validHistory(2:end) ];
    end
    
    if ~wasErr
        mismatches = find(abs(mean(simHistory, 1) - predictedHistory) > tol, 1);
        if ~isempty(mismatches)
            error([ 'oops!!  ' num2str(mismatches) ' mismatches' ])
        end
    end
    
    simCounter = simCounter+1;
end