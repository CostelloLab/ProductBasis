% ------------------- EQUATION-BUILDING TESTS ---------------------
% Here we test that buildFt() correctly builds a ladder of
% equations that can be used to simulate a network, or a mixed ensemble of
% networks.

%% test the Repressilator network (from Elowitz and Leibler, Nature 2000)

logicTables = { ...
    { [3], [ 1 0 ], [ 0 ] }, ...            % variable 1 -- evolves to NOT x3
    { [1], [ 1 0 ], [ 0 ] }, ...            % var 2 <-- NOT x1
    { [2], [ 1 0; 1 1 ], [ 0; 1 ] }    };   % var 3 <-- NOT (2-g(t))

M = polynomialModelT(logicTables);

M = buildFt(M);

if sum(sum(M.xs ~= [ 1 0 0; 0 1 0; 0 0 1; 0 0 0 ])) == 0 ...
        && sum(sum(abs(M.fs{1} - [ 0 0 -1 1; -1 0 0 1; 0 -1 0 1; 0 0 0 1 ]) > 1e-6)) == 0 ...
        && sum(sum(abs(M.fs{2} - [ 0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0 ]) > 1e-6)) == 0 ...
        && (sum(sum(abs(M.gs - [ 0; 1 ])))) == 0
    disp('passed repressilator test')
else
    error('*** ERROR!!! on repressilator test ***')
end



%% test a 2nd-order logic gates (i.e. the only extra variables are x_{34}
% and x_{12}, so the series only involves up to 2nd-order variables)

logicTables = { ...
    { [3 4], [ 0 1 1 0; 0 0 1 1 ], [ 0 0; 2 0 ] }, ...  % 1 <-- x3 XOR x4 + x3*g1^2
    { [4], [ 1 0 ], [ 0 0 ] }, ...                      % 2 <-- NOT x_4
    { [1 2], [ 0 0 0 1 ], [ 0 0 ] }, ...                % 3 <-- x1 AND x2
    { [1 2], [ 0 0 1 1; 0 1 0 0 ], [ 0 0; 1 1 ] }   };  % 4 <-- x1 OR x2*g1*g2

M = polynomialModelT(logicTables, [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1 ] == 1);

M = buildFt(M);

if sum(sum(M.xs ~= [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 1 1; 0 0 0 0; 1 1 0 0 ])) == 0 ...
        && sum(sum(abs(M.fs{1} - [ 0 0 1 1 -2 0 0; 0 0 0 -1 0 1 0; 0 0 0 0 0 0 1; 1 0 0 0 0 0 0; ...
                                   0 0 0 0 0 0 1; 0 0 0 0 0 1 0; 0 0 1 0 -1 0 0 ]) > 1e-6)) == 0 ...
        && sum(sum(abs(M.fs{2} - [ 0 0 1 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; ...
                                   0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 1 0 -1 0 0 ]) > 1e-6)) == 0 ...
        && sum(sum(abs(M.fs{3} - [ 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 1 0 0 0 0 -1; ...
                                   0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 ]) > 1e-6)) == 0 ...
        && (sum(sum(abs(M.gs - [ 0 0; 2 0; 1 1 ])))) == 0
    disp('passed 2nd-order logic gate test')
else
    error('*** ERROR!!! on 2nd-order logic gate test ***')
end



%% run a number of tests on a 4th-order logic system, involving many
% variables up through x_{1234}

% #1: test our simBoolModel() routine correctly simulates a single
% network

logicTables = { ...
    { [3], [ 1 0 ], [ 0 0 ] }, ...                          % 1 <-- NOT x3
    { [1 4], [ 0 0 1 1; 0 1 0 -1 ], [ 0 0; 1 2 ] }, ...     % 2 <-- x1*g1*g2^2 XOR x4
    { [1 2], [ 0 0 0 1; 1 1 1 1 ], [ 0 0; 0 1 ] }, ...      % 3 <-- x1 AND x2 + g2
    { [3 2], [ 0 0 1 1; 0 1 0 0 ], [ 0 0; 0 1 ] }   };      % 4 <-- x2 OR x3*g1

M = polynomialModelT(logicTables);

[ bs, xs ] = simBoolModelT(logicTables, [0 1 1 0], [0 1 0 1], 10);

if sum(sum(bs ~= [ 0 1 1 0; 0 0 0 1; 1 1 0 0; 1 1 1 1; 0 0 1 1; 0 1 0 1;
                  1 1 0 1; 1 0 1 1; 0 0 0 1; 1 1 0 0; 1 1 1 1 ])) == 0 ...
        && sum(xs' ~= [ 0 0 0 1 0 1 1 0 0 0 1 ]) == 0
    disp('passed 4th-order logic gate simulation')
else
    error('*** ERROR!!! on 4th-order logic gate simulation ***')
end


%% #2: do our usual test of buildFt(), carrying the
% calculation to completion

M = buildFt(M);

[fM_vecs, fM_vals] = eig(M.fs{1}, 'vector');
fM_xs = M.xs;

if sum(sum(M.xs ~= [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1;
                         0 0 0 0; 1 0 0 1; 1 1 0 0; 0 1 1 0;
                         0 0 1 1; 1 0 1 0; 1 0 1 1; 1 1 0 1;
                         1 1 1 0; 0 1 0 1; 0 1 1 1; 1 1 1 1 ])) == 0 ...
        && sum(sum(abs(M.fs{1} - [ ...
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
                        0  0  0  0  0  0  1  0  0  0  0 -1 -1  0  0  1   ]) > 1e-6)) == 0 ...
        && (sum(sum(abs(M.gs - zeros(1, 0))))) == 0
    disp('passed 4th-order logic gate test')
else
    error('*** ERROR!!! on 4th-order logic gate test ***')
end


%% #3: test that evolveF() using the output of
% buildFt() correctly predicts the evolution of (in this
% case) one network

if sum(round(evolveFt(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(4, :), 10)) ...
                    ~= [ 1 1 1 1 0 1 1 1 1 1 1 ]) == 0 && ...
        sum(round(evolveFt(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(1, :), 10)) ...
                    ~= [ 0 1 1 0 1 1 0 0 1 1 0 ]) == 0 && ...
        sum(round(evolveFt(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(6, :), 10)) ...
                    ~= [ 0 1 1 0 0 1 0 0 1 1 0 ]) == 0 && ...
        sum(round(evolveFt(M, [ 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 ]', M.xs(15, :), 10)) ...
                    ~= [ 0 0 0 0 0 1 0 0 0 0 0 ]) == 0
    disp('passed 4th-order logic gate evolution')
else
    error('*** ERROR!!! on 4th-order logic gate evolution ***')
end




%% FLY DEVELOPMENT TEST -- here we find the attractors of a much larger
% (16-node) network, taken from Albert et al. 2003

FlyDevelopmentTest();
%loadEMT();     % alternative (bigger: ~80 node) network:  Steinway et al. 2014
numBs = size(M.xs, 2);
fullTest = false;

flyStartingStates = round(rand(1,length(logicTables)));
flySim = simBoolModel(logicTables, flyStartingStates, zeros(1,length(logicTables)), 30);
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
    [ M, stats, dbgIndices ] = buildFt(M, 10000, false, allBinaries);
    allBinaries = allBinaries(:, dbgIndices);
else
    [ M, stats ] = buildFt(M, 10000, false);       % stats = # new vars, # new constraints, total # variables, total # constraints
end

    % check M.fs to make sure that they are consistent with
    % the steady state that we observed by simulation

xCoefs_0 = zeros(size(M.fs{1}, 2), 1);
xCoefs_f = zeros(size(M.fs{1}, 1), 1);
for loopX = 1:size(M.fs{1}, 2)
    xCoefs_0(loopX) = prod(flySim(end-1, M.xs(loopX, :) ~= 0));
    if loopX <= size(M.fs{1}, 1)
        xCoefs_f(loopX) = prod(flySim(end, M.xs(loopX, :) ~= 0));
    end
end
if max(abs(M.fs{1}*xCoefs_0 - xCoefs_f)) > 1.e-7
    error('*** ERROR!!! on fly development test ***')
end
