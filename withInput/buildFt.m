% buildF() approximates the transition rules of a Boolean
% network.  Each variable 'x' can be one of the Boolean variables 'b' or a
% product of them:  for example x13 is b1*b3.  A transition function is
% denoted 'f'.

function [ M, stats ] = buildFt(M, maxVars, isContinuous)


    % our initial features are just the fundamental Boolean state variables

backspaceString = '';

numBs = size(M.xs, 2);
numGs = size(M.gs, 2);
stats = zeros(0, 2);

if ~exist('maxVars', 'var')
    maxVars = min(2^numBs, 1e4);
end

if ~exist('isContinuous', 'var')
    isContinuous = false;
end

if isempty(M.fs)
    numFs = 0;
else
    numFs = size(M.fs{1}, 1);
end


    % first iteration:  find evolution equations for the fundamental variables
    % 
    % 2nd-Nth iterations:  find evolution equations for the higher-order
    % correlation variables that were introduced in earlier iterations
    
while numFs < size(M.xs, 1)
    
    stats(end+1, :) = [ 0 0 ];
    stats(end, 1) = size(M.xs, 1);
    
    
        % PART 2
        % add new variables if we've found all of the degeneracies

    
        % loop over each new variable we need to add an evolution equation
        % for
    
    if size(M.xs, 1) > maxVars
        disp('stats:')
        disp(stats)
        error('buildF:varLimit', 'variable limit reached')
    end
    
    numFsToAdd = size(M.xs, 1)-numFs;
    
    newFs = cell(1, length(M.fs));
    for loopG = 1:length(M.fs)
        newFs{loopG} = zeros(numFsToAdd, size(M.xs, 1));
    end
    newFerr = newFs;
    
    for loopNewF = 1:numFsToAdd
        
        printProgress([ num2str(loopNewF) ' / ' num2str(numFsToAdd) ' new features']);
        
        
            % xf identifies the variable as a product of state variables
            % (i.e. xf = [ 2 4 5 ] is the variable x_245)
        
        newX = M.xs(numFs+loopNewF, :);
        if ~isContinuous
            xfPoly = { 1, 0, false(1, numBs), zeros(1, numGs) };
            for loopBaseF = find(newX)
                xfPoly = multiplyTPoly(xfPoly, M.baseFs{loopBaseF});
            end
        else
            xfPoly = { zeros(0, 1), zeros(0, 1), false(0, numBs), zeros(0, 1) };
            for loopBaseF = find(newX)
                otherVars = { 1, 0, newX, zeros(1, numGs) };
                otherVars{3}(loopBaseF) = false;
                xfPoly = addTPoly(xfPoly, multiplyTPoly(otherVars, M.baseFs{loopBaseF}));
            end
        end
        
        for loopPolyTerm = 1:size(xfPoly{3}, 1);
            varID = find((sum(M.xs == repmat(xfPoly{3}(loopPolyTerm, :), size(M.xs, 1), 1), 2) == numBs)', 1, 'first');
            gID = find((sum(M.gs == repmat(xfPoly{4}(loopPolyTerm, :), size(M.gs, 1), 1), 2) == numGs)', 1, 'first');
            if isempty(varID)
                varID = size(M.xs, 1)+1;
                M.xs(varID, :) = xfPoly{3}(loopPolyTerm, :);
            end
            if isempty(gID)
                gID = size(M.gs, 1)+1;
                M.gs(gID, :) = xfPoly{4}(loopPolyTerm, :);
                M.fs{1, gID} = zeros(numFs, size(M.xs, 1));
                M.ferr{1, gID} = M.fs{1, gID};
                newFs{1, gID} = zeros(numFsToAdd, size(M.xs, 1));
                newFerr{1, gID} = newFs{1, gID};
            end
            
            newFs{gID}(loopNewF, varID) = xfPoly{1}(loopPolyTerm);
            newFerr{gID}(loopNewF, varID) = xfPoly{2}(loopPolyTerm);
        end
    end
    stats(end, 2) = numFsToAdd;
    
    numXs = size(M.xs, 1);
    for loopG = 1:length(M.fs)
        M.fs{loopG} = [ M.fs{loopG}  zeros(numFs, numXs-size(M.fs{loopG}, 2)); newFs{loopG}  zeros(size(newFs{loopG}).*[1 -1]+[0 numXs]) ];
        M.ferr{loopG} = [ M.ferr{loopG}  zeros(numFs, numXs-size(M.ferr{loopG}, 2)); newFerr{loopG}  zeros(size(newFerr{loopG}).*[1 -1]+[0 numXs]) ];
    end
    numFs = numFs + size(newFs{1}, 1);
end

printProgress(sprintf('\b'));
    

    function printProgress(progressString)
        disp([backspaceString progressString])
        backspaceString = sprintf(repmat('\b', 1, length(progressString)+1));
    end
end