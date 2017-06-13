% evolveTransitionMatrix() returns the time-evolution of a variable of the
% user's choice (bReadout), given an initial state b0 as a row vector.

function [ bHistory, readout ] = simTransitionMatrix(logicTables, b0, bReadout, tMax) %, isContinuous)

% if ~exist('isContinuous', 'var')
%     isContinuous = false;
% end

numBs = length(logicTables);
bHistory = false(tMax+1, numBs);
bHistory(1, :) = b0;
readout = zeros(tMax+1, 1);
readout(1) = prod(b0(bReadout==1));

% if isContinuous
%     [ ~, bHistory ] = ode45(bEvolutionFunction, 1:tMax, b0');
% end

for t = 1:tMax
%     if ~isContinuous
    for b = 1:numBs
        bInput = bHistory(t, logicTables{b}{1});
        powerOfTwos = 2.^((length(bInput)-1):-1:0)';
        bHistory(t+1, b) = (rand() < logicTables{b}{2}(bInput*powerOfTwos+1));
    end
%     end
    readout(t+1) = prod(bHistory(t+1, bReadout==1));
end
% 
% 
%     function c = bEvolutionFunction(b)
%         powerOfTwos = 2.^((length(b)-1):-1:0)';
%     end

end