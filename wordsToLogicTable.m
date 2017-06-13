% wordsToLogicTable() turns an expression such as 'p53 AND not DEAD' into a
% logic table, given variable names (in this case 'p53' and 'DEAD')

function [ usedVars, logicTable ] = wordsToLogicTable(varNames, transitionRule, globalCode)


    % compute the new name for each variable so that they refer to their
    % respective columns in the inputTable (defined below);
    % also change 'and' to '&', etc.

subs = { '&', '|', '~','&', '|', '~' };
numSubs = length(subs);
varNames(end+1:end+numSubs) = { 'and', 'or', 'not', 'AND', 'OR', 'NOT' };

varLengths = zeros(1, length(varNames));
for loopVar = 1:length(varNames)
    varLengths(loopVar) = length(varNames{loopVar});
end
[ ~, varsInOrder ] = sort(varLengths, 'descend');

usedVars = [];
numUsedVars = 0;

transitionRule = [ ' ' transitionRule ' ' ];        % so we can check the letter after every word

for loopVar = varsInOrder
    textEndpoints = [ strfind(transitionRule, varNames{loopVar}) - 1, length(transitionRule) ];
    textBeginnings = [ 1, textEndpoints(1:end-1) + length(varNames{loopVar}) + 1 ];
    
    firstLetters = transitionRule(textBeginnings(2:end));
    idx = ((firstLetters == ' ') | (firstLetters == ')'));
    if sum(~idx) > 0
    end
    textEndpoints = textEndpoints([ idx true ]);
    textBeginnings = textBeginnings([ true idx ]);
    
    if length(textEndpoints) > 1
        if loopVar <= length(varNames)-numSubs
            numUsedVars = numUsedVars + 1;
            usedVars(1, numUsedVars) = loopVar;
            newVarName = [ 'inputTable(:, ' num2str(numUsedVars) ')' ];
        else
            newVarName = subs{loopVar-length(varNames)+numSubs};
        end
        recodedRule = transitionRule(textBeginnings(1):textEndpoints(1));
        for loopInsertion = 2:length(textEndpoints)
            recodedRule = [ recodedRule newVarName ...
                transitionRule(textBeginnings(loopInsertion):textEndpoints(loopInsertion)) ];
        end
        transitionRule = recodedRule;
    end
end


    % construct a table containing all possible inputs

inputTable = zeros(2^numUsedVars, numUsedVars);
allNums = (0:(2^numUsedVars-1))';
for loopVar = 1:numUsedVars
    inputTable(:, loopVar) = mod(bitshift(allNums, loopVar-numUsedVars), 2);
end

if exist('globalCode', 'var')
    eval(globalCode);
    %transitionRule = [ globalCode transitionRule ];
end


    % evaluate the rule as a MATLAB expression

try
    logicTable = eval(transitionRule);
catch
    error('wordsToLogicTable():  problem with inputs')
end

end