% multiplyPoly() returns the products of two polynomials p1 and p2

function pProd = multiplyTPoly(p1, p2)

numCoefs1 = size(p1{1}, 1);
numCoefs2 = size(p2{1}, 1);

pProd = { reshape(p2{1}*(p1{1}'), [], 1), ...
    reshape(abs(p2{1})*(p1{2}') + p2{2}*abs(p1{1}'), [], 1), ...
    false(numCoefs1*numCoefs2, size(p1{3}, 2)), ...
    zeros(numCoefs1*numCoefs2, size(p1{4}, 2)) };

idxCounter = 0;
for c1 = 1:numCoefs1
    for c2 = 1:numCoefs2
        idxCounter = idxCounter+1;
        pProd{3}(idxCounter, :) = p1{3}(c1, :) | p2{3}(c2, :);
        pProd{4}(idxCounter, :) = p1{4}(c1, :) + p2{4}(c2, :);
    end
end

pProd = reducePolyT(pProd);