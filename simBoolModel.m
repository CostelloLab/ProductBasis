% evolveTransitionMatrix() returns the time-evolution of a variable of the
% user's choice (bReadout), given an initial state b0 as a row vector.

function [ bHistory, readout ] = simBoolModel(logicTables, b0, bReadout, tMax, subTsteps)

if ~exist('subTsteps', 'var')
    subTsteps = 1;
end

numBs = length(logicTables);
bHistory = false(tMax+1, numBs);
bHistory(1, :) = b0;
bState = false(2, numBs);
readout = zeros(tMax+1, 1);
readout(1) = prod(b0(bReadout==1));

powersOfTwos = cell(numBs);
for b = 1:numBs
    powerOfTwos{b} = 2.^((length(logicTables{b}{1})-1):-1:0)';
end

for t = 1:tMax
    bState(1, :) = bHistory(t, :);
    brands = rand(subTsteps, numBs);
    for subT = 1:subTsteps
        for b = 1:numBs
            bInput = bState(1, logicTables{b}{1});
            bState(2, b) = (brands(subT, b) < logicTables{b}{2}(bInput*powerOfTwos{b}+1));
        end
        bState(1, :) = bState(2, :);
    end
    bHistory(t+1, :) = bState(2, :);
    readout(t+1) = prod(bHistory(t+1, bReadout==1));
end

end