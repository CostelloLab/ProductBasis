function [ newVars, newTransitionRules ] = mutateNetwork(allVars, transitionRules)

numVars = length(allVars);

newVars = [ allVars allVars ];
newTransitionRules = [ transitionRules transitionRules ];

for loopVar = 1:numVars
    newVars{numVars+loopVar} = [ 'WT_' allVars{loopVar} ];
    newTransitionRules{numVars+loopVar} = newVars{numVars+loopVar};
    newTransitionRules{loopVar} = [ '(' transitionRules{loopVar} ') and ' newVars{numVars+loopVar} ];
end

end