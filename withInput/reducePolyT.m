% reducePoly() reduces a polynomial by a) combining terms in the same
% variable, and b) eliminating terms with a zero coefficient

function pReduced = reducePolyT(pUnreduced)

pReduced = cell(1, 4);


    % sort the terms in order of their variable & power in g_i

[ pReduced{4}, indexOrder ] = sortrows(pUnreduced{4});
pReduced{1} = pUnreduced{1}(indexOrder, :);
pReduced{2} = pUnreduced{2}(indexOrder, :);
pReduced{3} = pUnreduced{3}(indexOrder, :);

[ pReduced{3}, indexOrder ] = sortrows(pReduced{3});
pReduced{1} = pReduced{1}(indexOrder, :);
pReduced{2} = pReduced{2}(indexOrder, :);
pReduced{4} = pReduced{4}(indexOrder, :);

if isempty(pReduced{1})
    return
end


    % find terms in the same variable (which are now consecutive in our
    % terms list), and combine them

uniqueCoefs = [ true; sum(abs(diff(pReduced{3}, 1, 1)), 2) ~= 0 | sum(abs(diff(pReduced{4}, 1, 1)), 2) ~= 0 ];
nonuniqueCoefs = find(~uniqueCoefs);
for c1 = 1:length(nonuniqueCoefs)
    if c1 == 1
        lastUniqueCoef = nonuniqueCoefs(1)-1;
    elseif nonuniqueCoefs(c1-1) < nonuniqueCoefs(c1)-1
        lastUniqueCoef = nonuniqueCoefs(c1)-1;
    end
    
    pReduced{1}(lastUniqueCoef) = pReduced{1}(lastUniqueCoef) + pReduced{1}(nonuniqueCoefs(c1));
    pReduced{2}(lastUniqueCoef) = pReduced{2}(lastUniqueCoef) + pReduced{2}(nonuniqueCoefs(c1));
end

pReduced{1} = pReduced{1}(uniqueCoefs, :);
pReduced{2} = pReduced{2}(uniqueCoefs, :);
pReduced{3} = pReduced{3}(uniqueCoefs, :);
pReduced{4} = pReduced{4}(uniqueCoefs, :);


    % finally, eliminate terms with a zero coefficient

nonzeroElements = (abs(pReduced{1}) > pReduced{2});
pReduced{1} = pReduced{1}(nonzeroElements, :);
pReduced{2} = pReduced{2}(nonzeroElements, :);
pReduced{3} = pReduced{3}(nonzeroElements, :);
pReduced{4} = pReduced{4}(nonzeroElements, :);
