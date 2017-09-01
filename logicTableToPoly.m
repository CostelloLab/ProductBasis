% logicTableToPoly() fits a 2^N-length list of Boolean variables (N is the
% number of variables) with a polynomial in those variables.
% 
% theTable is a column vector of Boolean outcomes; thePoly{1} is the
% coefficients of the polynomial fitting theTable as a column vector, and
% each row of thePoly{2} gives the power of each of the Boolean variables
% (including 1) corresponding to one coefficient

function thePoly = logicTableToPoly(outcomes)


    % calculate the coefficients of the polynomial that describes
    % the transition table
    % 
    % one column of theTable per variable (e.g. 1 or x3 or x24)
    % one row of theTable per input (e.g. x1 = 1, x2 = 0, ...)

if size(outcomes, 1) == 1
    outcomes = outcomes';
end

numVars = round(log2(length(outcomes)));

polyPowers = false(length(outcomes), numVars);
polyBasis = ones(length(outcomes));
count = 0:(length(outcomes)-1);
for loopRedParam = 0:(numVars-1)
    polyPowers(:, loopRedParam+1) = mod(bitshift(count, 1-numVars+loopRedParam), 2);
    singleParamFeature = mod(floor(count/2^loopRedParam), 2);
    polyBasis(:, singleParamFeature == 1) = polyBasis(:, singleParamFeature == 1) ...
        .* (singleParamFeature'*ones(1, length(outcomes)/2));
end

polyCoefficients = polyBasis\outcomes;
polyErrors = (abs(outcomes-polyBasis*polyCoefficients)+eps(outcomes))*cond(polyBasis);
nonzeroEntries = abs(polyCoefficients) > polyErrors;
thePoly = { polyCoefficients(nonzeroEntries, :), polyErrors(nonzeroEntries, :), polyPowers(nonzeroEntries, :) == 1 };


end