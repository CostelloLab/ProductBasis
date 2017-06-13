%Input true to run simulation with mutations or false to not
doingMutations = false; %true;
%Input the percent of wild type (non mutated) ex 0.99 for a 1% mutation
wildTypeRate = 0.99;
%input true to run simulation with Monte Carlo or false to not
compareMonteCarlo = false; %true;
% input the amount of runs to generate Monte Carlo
numMonteCarloRuns = 10^4;
% input the number of time steps to simulate
timeSteps= 50;
% input the number of the variable to follow
toLookAt = {[4] [9] [4 9]}; %{ [ 4 9 ] };

% list all variables being simulated
allVars = {'CD45', 'CD4', 'TCRLIG', 'TCRBIND', 'PAGCSK', 'LCK', 'FYN', 'CCBL', 'TCRPHOS', 'RLK', 'ZAP70', ...
    'LATPHOSP', 'GADS', 'ITK', 'IP3', 'CA', 'CALCIN', 'SLP76', 'PLCGB', 'PLCGACT', 'DAG', 'GRB2SOS', 'RASGRP1', ...
    'PKCTH', 'RAS', 'RAF', 'SEK', 'MEK', 'JNK', 'IKKB', 'RSK', 'ERK', 'FOS', 'JUN', 'IKB', 'CREB', 'CRE', 'AP1', ...
    'NFKB', 'NFAT'}; 
% list boolen rules for each variable to turn on, in order with allVars
transitionRules = { ...
  'CD45', 'CD4', 'TCRLIG', 'TCRLIG AND (NOT CCBL)', 'FYN OR (NOT TCRBIND)', 'CD45 AND CD4 AND (NOT PAGCSK)', ...
  '(TCRBIND AND CD45) OR (LCK AND CD45)', 'ZAP70', '(LCK AND TCRBIND) OR FYN', 'FYN', ...
  'FYN AND TCRPHOS AND (NOT CCBL)', 'ZAP70', 'LATPHOSP', 'ZAP70 AND SLP76', 'PLCGACT', 'IP3', 'CA', 'GADS', ...
  'LATPHOSP', '(PLCGB AND SLP76 AND ZAP70 AND ITK) OR (PLCGB AND SLP76 AND ZAP70 AND RLK)', 'PLCGACT', ...
  'LATPHOSP', 'PKCTH', 'DAG', 'GRB2SOS OR RASGRP1', 'RAS', 'PKCTH', 'RAF', 'SEK', 'PKCTH', 'ERK', 'MEK', 'ERK', ...
  'JNK', '(NOT IKKB)', 'RSK', 'CREB', 'FOS AND JUN', '(NOT IKB)', 'CALCIN'}; 
%ready to run
if doingMutations
    [allVars, transitionRules ] = mutateNetwork(allVars, transitionRules);
end
logicTables = cell(1, length(allVars));
for loopVar = 1:length(allVars)
    [ usedVars, theTable ] = wordsToLogicTable(allVars, transitionRules{loopVar});
    logicTables{loopVar} = { usedVars, theTable };
end


if doingMutations
    bHistory=simTransitionMatrix(logicTables,[rand(1,0.5*length(allVars)) < .5, rand(1,0.5*length(allVars)) < wildTypeRate]==1,zeros(1,length(allVars)),30);
else
    bHistory=simTransitionMatrix(logicTables, rand(1, length(allVars)) < .5, zeros(1,length(allVars)),30);
end
plotBoolEvolution();

toLookAtLogical = false(length(toLookAt), length(logicTables));
for loopVar = 1:length(toLookAt)
    toLookAtLogical(loopVar, toLookAt{loopVar}) = true;
end
M = polynomialModel(logicTables, toLookAtLogical);

if compareMonteCarlo
    allSims= zeros(numMonteCarloRuns, (1+timeSteps));
    montC=zeros(length(toLookAt), (timeSteps +1));
    err = zeros(length(toLookAt), (timeSteps +1));
    for montloop = 1:length(toLookAt)
        for i=1:numMonteCarloRuns;
            [ ~, tmp ] = simTransitionMatrix(logicTables,[round(rand(1,0.5*length(allVars))), ones(1,0.5*length(allVars))]==1,toLookAtLogical(montloop,:),timeSteps);
            allSims(i, :) = tmp;
        end
        montC(montloop,:)= mean(allSims,1);
     end  
    err = (sqrt((montC-(montC.^2))/(numMonteCarloRuns-1)));
end

tic;
M = buildF(M, 20000, 20000);
toc

disp(['# variables:  ' num2str(size(M.xs, 1)) ' + ' num2str(diff(size(M.fs))) ])


allVarsToLookAt = find(sum(toLookAtLogical, 1) > 0);
figTitle = allVars{allVarsToLookAt(1)};
for loopName = 2:length(allVarsToLookAt)
    figTitle = [ figTitle ', ' allVars{allVarsToLookAt(loopName)} ];
end

figure, hold on
legendStrings = {};
for loopVar = 1:length(toLookAt)    
    if doingMutations
        plot(0:timeSteps, evolveTransitionMatrix(M, 0.5.^(sum(M.xs(:, 1:0.5*length(allVars)),2)) .* 1.00.^(sum(M.xs(:, 41:end),2)), toLookAtLogical(loopVar, :), timeSteps))
        plot(0:timeSteps, evolveTransitionMatrix(M, 0.5.^(sum(M.xs(:, 1:0.5*length(allVars)),2)) .* wildTypeRate.^(sum(M.xs(:, 41:end),2)), toLookAtLogical(loopVar, :), timeSteps))
    else    
        plot(0:timeSteps, evolveTransitionMatrix(M, 0.5.^(sum(M.xs,2)), toLookAtLogical(loopVar, :), timeSteps))
    end
    if compareMonteCarlo
        errorbar(0:timeSteps, montC(loopVar, :), err(loopVar, :), 'o')
    end  
    legendStrings{end+1} = allVars{toLookAt{loopVar}(1)};
    for loopName = 2:length(toLookAt{loopVar})
        legendStrings{end} = [ legendStrings{end} ' AND ' allVars{toLookAt{loopVar}(loopName)} ];
    end
    if doingMutations
        legendStrings{end+1} = [ legendStrings{end} ' mutated' ];
    end
    if compareMonteCarlo
        legendStrings{end+1} = [ legendStrings{end} ' Monte Carlo' ];
    end
end

if doingMutations
    title([ figTitle ' mutated' ], 'FontSize', 12)
else
    title([ figTitle ], 'FontSize', 12)
end
xlabel('time \rightarrow')
ylabel('population fraction')
xlim([0 timeSteps])
legend(legendStrings)
set(gcf, 'Position', [0, 3000, 350, 200])