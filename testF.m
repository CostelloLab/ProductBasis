% ------------------- EQUATION-BUILDING TESTS ---------------------
% Here we test that buildF() correctly builds a ladder of
% equations that can be used to simulate a network, or a mixed ensemble of
% networks.

%% test the Repressilator network (from Elowitz and Leibler, Nature 2000)

logicTables = { ...
    { [3], [ 1 0 ] }, ...       % variable 1 -- evolves to NOT variable 2
    { [1], [ 1 0 ] }, ...       % var 2 <-- NOT 3
    { [2], [ 1 0 ] }    };      % var 3 <-- NOT 1

M = polynomialModel(logicTables);

[ M, tStart ] = buildF(M);

if sum(sum(M.xs ~= [ 1 0 0; 0 1 0; 0 0 1; 0 0 0 ])) == 0 ...
        && sum(sum(abs(M.fs - [ 0 0 -1 1; -1 0 0 1; 0 -1 0 1; 0 0 0 1 ]) > 1e-6)) == 0 ...
        && tStart == 1
    disp('passed repressilator test')
else
    error('*** ERROR!!! on repressilator test ***')
end



%% test a 2nd-order logic gates (i.e. the only extra variables are x_{34}
% and x_{12}, so the series only involves up to 2nd-order variables)

logicTables = { ...
    { [3 4], [ 0 1 1 0 ] }, ...     % 1 <-- 3 XOR 4
    { [4], [ 1 0 ] }, ...           % 2 <-- NOT 4
    { [1 2], [ 0 0 0 1 ] }, ...     % 3 <-- 1 AND 2
    { [1 2], [ 0 1 1 1 ] }   };     % 4 <-- 1 OR 2

M = polynomialModel(logicTables, [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1 ] == 1);

[ M, tStart ] = buildF(M);

if sum(sum(M.xs ~= [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 1 1; 0 0 0 0; 1 1 0 0 ])) == 0 ...
        && sum(sum(abs(M.fs - [ 0 0 1 1 -2 0 0; 0 0 0 -1 0 1 0; 0 0 0 0 0 0 1; 1 1 0 0 0 0 -1; ...
                                   0 0 0 0 0 0 1; 0 0 0 0 0 1 0; 0 0 1 0 -1 0 0 ]) > 1e-6)) == 0 ...
        && tStart == 1
    disp('passed 2nd-order logic gate test')
else
    error('*** ERROR!!! on 2nd-order logic gate test ***')
end



%% run a number of tests on a 4th-order logic system, involving many
% variables up through x_{1234}

% #1: test our simTransitionMatrix() routine correctly simulates a single
% network

logicTables = { ...
    { [3], [ 1 0 ] }, ...           % 1 <-- NOT 3
    { [1 4], [ 0 1 1 0 ] }, ...     % 2 <-- 1 XOR 4
    { [1 2], [ 0 0 0 1 ] }, ...     % 3 <-- 1 AND 2
    { [3 2], [ 0 1 1 1 ] }   };     % 4 <-- 2 OR 3

M = polynomialModel(logicTables);

x0 = M.xs;
fs0 = M.fs;
[ bs, xs ] = simTransitionMatrix(logicTables, [0 1 1 0], [0 1 0 1], 10);

if sum(sum(bs ~= [ 0 1 1 0; 0 0 0 1; 1 1 0 0; 1 1 1 1; 0 0 1 1; 0 1 0 1;
                  1 1 0 1; 1 0 1 1; 0 0 0 1; 1 1 0 0; 1 1 1 1 ])) == 0 ...
        && sum(xs' ~= [ 0 0 0 1 0 1 1 0 0 0 1 ]) == 0
    disp('passed 4th-order logic gate simulation')
else
    error('*** ERROR!!! on 4th-order logic gate simulation ***')
end


%% #2: do our usual test of buildF(), carrying the
% calculation to completion

[ M, tStart ] = buildF(M);

[fM_vecs, fM_vals] = eig(M.fs, 'vector');
fM_xs = M.xs;

if tStart == 1 ...
        && sum(sum(M.xs ~= [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1;
                         0 0 0 0; 1 0 0 1; 1 1 0 0; 0 1 1 0;
                         0 0 1 1; 1 0 1 0; 1 0 1 1; 1 1 0 1;
                         1 1 1 0; 0 1 0 1; 0 1 1 1; 1 1 1 1 ])) == 0 ...
        && sum(sum(abs(M.fs - [ ...
                        0  0 -1  0  1  0  0  0  0  0  0  0  0  0  0  0;
                        1  0  0  1  0 -2  0  0  0  0  0  0  0  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0;
                        0  1  1  0  0  0  0 -1  0  0  0  0  0  0  0  0;
                        0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0;
                        0  1  0  0  0  0  0 -1  0  0  0  0  0  0  0  0;
                        1  0  0  1  0 -2  0  0 -1 -1  2  0  0  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0 -1  0  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0  0 -1  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0  0 -1  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0 -2 -1  1 -1  2;
                        0  0  0  0  0  0  1  0  0  0  0 -1 -1  0  0  1;
                        0  0  0  0  0  0  1  0  1  1 -2 -2 -1  1 -1  2;
                        0  0  0  0  0  0  1  0  0  0  0 -1  0  0  0  0;
                        0  0  0  0  0  0  1  0  0  0  0 -1 -1  0  0  1   ]) > 1e-6)) == 0
    disp('passed 4th-order logic gate test')
else
    error('*** ERROR!!! on 4th-order logic gate test ***')
end


%% #3: test that evolveTransitionMatrix() using the output of
% buildF() correctly predicts the evolution of (in this
% case) one network

if sum(round(evolveTransitionMatrix(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(4, :), 10)) ...
                    ~= [ 1 1 1 1 0 1 1 1 1 1 1 ]) == 0 && ...
        sum(round(evolveTransitionMatrix(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(1, :), 10)) ...
                    ~= [ 0 1 1 0 1 1 0 0 1 1 0 ]) == 0 && ...
        sum(round(evolveTransitionMatrix(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(6, :), 10)) ...
                    ~= [ 0 1 1 0 0 1 0 0 1 1 0 ]) == 0 && ...
        sum(round(evolveTransitionMatrix(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(15, :), 10)) ...
                    ~= [ 0 0 0 0 0 1 0 0 0 0 0 ]) == 0
    disp('passed 4th-order logic gate evolution')
else
    error('*** ERROR!!! on 4th-order logic gate evolution ***')
end




%% ------------------- EQUATION-PRUNING TESTS ---------------------
% Now we test buildF()'s ability to prune a set of
% equations even as it builds them, in order to analyze the final behavior
% (i.e. attractors) of much larger networks

% #5: test the attractors found by buildF()

M.xs = x0;
M.fs = fs0;
M.cxs = false(0, size(x0, 2));
M.cs = {};
M.ts = zeros(0, 1);
[ M, tStart ] = buildF(M, 20, 0, .1 );

fM_idxs_1 = zeros(1, size(M.xs, 1));
for loopX = 1:size(M.xs, 1)
    fM_idxs_1(loopX) = find(sum(ones(size(fM_xs, 1), 1)*M.xs(loopX, :) ~= fM_xs, 2) == 0);
end
[ ~, fM_idxs_2 ] = sort(pi*real(fM_vals)+imag(fM_vals));
fM_idxs_2 = fM_idxs_2(abs(fM_vals(fM_idxs_2)) > 1.e-4);
[ntM_vecs, ntM_vals] = eig(M.fs, 'vector');
[ ~, ntM_idxs_2 ] = sort(pi*real(ntM_vals)+imag(ntM_vals));
scaled_fM_vecs = fM_vecs(fM_idxs_1, fM_idxs_2);
scaled_fM_vecs = scaled_fM_vecs .* (ones(size(scaled_fM_vecs, 1), 1)*(ntM_vecs(2, ntM_idxs_2)./scaled_fM_vecs(2, :)));

[M.xs,newOrder] = sortrows(M.xs,1);
M.fs = M.fs(newOrder, newOrder);
M.ts = M.ts(newOrder);

Mxs_check = [ 0 0 0 1; 1 1 0 0; 0 0 1 1; 1 0 1 1; 1 1 0 1; 0 1 0 1; 1 1 1 1 ];
[Mxs_check,newOrder2] = sortrows(Mxs_check,1);
Mts_check = [ 3 4 1 2 3 3 2 ];
Mts_check = Mts_check(newOrder2);

if sum(sum(M.xs ~= Mxs_check)) == 0 ...
        && sum(abs(ntM_vals(ntM_idxs_2) - fM_vals(fM_idxs_2))) < 1.e-10 ...
        && sum(sum(abs(ntM_vecs(:, ntM_idxs_2) - scaled_fM_vecs))) < 1.e-10 ...
        && sum(M.ts ~= Mts_check') == 0 ...
        && tStart == 4
    disp('passed 4th-order logic gate test w/o transients')
else
    error('*** ERROR!!! on 4th-order logic gate test w/o transients ***')
end



%% TEST:  2 OR gates and a NOT gate -- test the attractors

logicTables = { ...
    { [2], [ 1 0 ] }, ...           % 1 <-- NOT 2
    { [1 3], [ 0 1 1 1 ] }, ...     % 2 <-- 1 OR 3
    { [1 2], [ 0 1 1 1 ] }   };     % 3 <-- 1 OR 2

M = polynomialModel(logicTables);

[ M, tStart ] = buildF(M, 10, 0, .1);

if sum(sum(M.xs ~= ([ 0 1 1 ] == 1))) == 0 ...
        && sum(sum(abs(M.fs - [ 1 ] > 1.e-12))) == 0 ...
        && M.ts == 5 ...
        && tStart == 5
    disp('passed OR-NOT test')
else
    error('*** ERROR!!! on OR-NOT test ***')
end



%% TEST:  find attractors of a 4-node probabilistic Repressilator-type network

logicTables = { ...
    { [4 1], [ 0 .9 .1 1 ] }, ...
    { [1 2], [ 0 .8 .2 1 ] }, ...
    { [2 3], [ 0 .7 .3 1 ] }, ...
    { [3 4], [ 0 .6 .4 1 ] }   };

M = polynomialModel(logicTables);

[ M, tStart ] = buildF(M, 16, 0, .6);

if sum(sum(M.xs ~= ([ 1 1 1 1 ] == 1))) == 0 ...
        && sum(sum(abs(M.fs - [ 1 ] > 1.e-10))) == 0 ...
        && M.ts == 4 ...
        && tStart == 4
    disp('passed 4-node PBN test')
else
    error('*** ERROR!!! on 4-node PBN test ***')
end




%% FLY DEVELOPMENT TEST -- here we find the attractors of a much larger
% (16-node) network, taken from Albert et al. 2003

FlyDevelopmentTest();
%loadEMT();     % alternative (bigger: ~80 node) network:  Steinway et al. 2014
numBs = size(M.xs, 2);
fullTest = false;

flyStartingStates = round(rand(1,length(logicTables)));
flySim = simTransitionMatrix(logicTables, flyStartingStates, zeros(1,length(logicTables)), 30);
if (sum(diff(flySim(end-1:end, :), 1, 1) ~= 0) ~= 0)
    error('** oops -- no convergence on fly test')
end

if fullTest
    allNums = 0:((2^numBs)-1);
    allBinaries = zeros(numBs, length(allNums));
    for cb = 1:numBs
        allBinaries(16-cb, :) = mod(bitshift(allNums, -cb+1), 2);
    end
end


% iteratively refine our equations

if fullTest
    [ M, tStart, stats, dbgIndices ] = buildF(M, 2000, 0, .01, false, allBinaries);
    allBinaries = allBinaries(:, dbgIndices);
else
    [ M, tStart, stats ] = buildF(M, 2000, 0, .01, false);       % stats = # new vars, # new constraints, total # variables, total # constraints
end

    % check M.fs and M.cs to make sure that they are consistent with
    % the steady state that we observed by simulation

xCoefs_0 = zeros(size(M.fs, 2), 1);
xCoefs_f = zeros(size(M.fs, 1), 1);
cCoefs = zeros(length(M.cs), 1);
for loopX = 1:size(M.fs, 2)
    xCoefs_0(loopX) = prod(flySim(end-1, M.xs(loopX, :) ~= 0));
    if loopX <= size(M.fs, 1)
        xCoefs_f(loopX) = prod(flySim(end, M.xs(loopX, :) ~= 0));
    end
end
for loopX = 1:length(M.cs)
    cCoefs(loopX) = prod(flySim(end-1, M.cxs(loopX, :) ~= 0));
end
if max(abs(M.fs*xCoefs_0 - xCoefs_f)) > 1.e-7
    error('*** ERROR!!! on fly development test ***')
end
polySums = zeros(length(M.cs), 1);
for loopConstraint = 1:length(M.cs)
    for loopVar = 1:size(M.cs{loopConstraint}{1}, 1)
        polySums(loopConstraint) = polySums(loopConstraint) + ...
            M.cs{loopConstraint}{1}(loopVar, 1)*prod(flySim(end-1, M.cs{loopConstraint}{2}(loopVar, :) ~= 0));
    end
    if abs(polySums - cCoefs) > 1.e-10
        error('*** ERROR!!! on fly development test ***')
    end
end

    % in the full test, we check every single possible state that has 
    % not yet been eliminated as not being in steady state, at each 
    % step in the calculation, and make sure that each equation works
    % for each state

if fullTest
    allXBinaries = zeros(size(M.xs, 1), size(allBinaries, 2));
    for loopX = 1:size(M.xs, 1)
        allXBinaries(loopX, :) = prod(allBinaries(M.xs(loopX, :) == 1, :), 1);
    end
    allXOutputs = M.fs*allXBinaries;
    badOutputs = (abs(allXOutputs) > 1.e-6) & (abs(allXOutputs-1) > 1e-6);
    if sum(sum(badOutputs)) ~= 0
        disp('problem here**')
    end

    for loopX = 1:size(M.xs, 1)
        loopXMask = M.xs(loopX, :)'*ones(1, size(allBinaries, 2));
        if sum(sum(loopXMask == 1 & allBinaries == 0, 1) == 0) == 0
            disp('extra variable here..')
        end
    end
end


    % check the eigenvalues of the final state-transition matrix -- should
    % all be of modulus 1

if min( max(abs(abs(eig(M.fs)))), max(abs(abs(eig(M.fs))-1)) ) > 1e-7
    error('*** ERROR!!! bad eigenvalues on fly development test ***')
else
    disp([ 'passed fly development test for starting time ', num2str(tStart) ])
end


    % temp -- make sure all of the other outcomes are forbidden by the
    % constraints

allNums = (0:((2^numBs)-1))';
allBinaries = zeros(length(allNums), numBs);
for cb = 1:numBs
    allBinaries(:, numBs+1-cb) = mod(bitshift(allNums, -cb+1), 2)';
end

for loopConstraint = 1:length(M.cs)
    polySum = 0*allNums;
    for loopTerm = 1:size(M.cs{loopConstraint}{1}, 1)
        variableProduct = 0*allNums+1;
        for loopBool = find(M.cs{loopConstraint}{2}(loopTerm, :))
            variableProduct = variableProduct .* allBinaries(:, loopBool);
        end
        polySum = polySum + M.cs{loopConstraint}{1}(loopTerm, 1)*variableProduct;
    end
    
    variableProduct = 0*allNums+1;
    for loopBool = find(M.cxs(loopConstraint, :))
        variableProduct = variableProduct .* allBinaries(:, loopBool);
    end
    
    toRemove = variableProduct ~= round(polySum);
    allNums(toRemove, :) = [];
    allBinaries(toRemove, :) = [];
end