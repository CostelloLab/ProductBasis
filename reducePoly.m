% reducePoly() reduces a polynomial by a) combining terms in the same
% variable, and b) eliminating terms with a zero coefficient

function pReduced = reducePoly(pUnreduced)

pReduced = cell(1, 2);


    % sort the terms in order of their variable

[ pReduced{2}, indexOrder ] = sortrows(pUnreduced{2});
pReduced{1} = pUnreduced{1}(indexOrder, :);

if isempty(pReduced{1})
    return
end


    % find terms in the same variable (which are now consecutive in our
    % terms list), and combine them

uniqueCoefs = [ true; sum(abs(diff(pReduced{2}, 1, 1)), 2) ~= 0 ];
nonuniqueCoefs = find(~uniqueCoefs);
for c1 = 1:length(nonuniqueCoefs)
    if c1 == 1
        lastUniqueCoef = nonuniqueCoefs(1)-1;
    elseif nonuniqueCoefs(c1-1) < nonuniqueCoefs(c1)-1
        lastUniqueCoef = nonuniqueCoefs(c1)-1;
    end
    
    pReduced{1}(lastUniqueCoef) = pReduced{1}(lastUniqueCoef) + pReduced{1}(nonuniqueCoefs(c1));
end

pReduced{1} = pReduced{1}(uniqueCoefs, :);
pReduced{2} = pReduced{2}(uniqueCoefs, :);


    % finally, eliminate terms with a zero coefficient

nonzeroElements = (abs(pReduced{1}) > 1.e-6);
pReduced{1} = pReduced{1}(nonzeroElements, :);
pReduced{2} = pReduced{2}(nonzeroElements, :);
