% multiplyPoly() returns the products of two polynomials p1 and p2

function pProd = multiplyPoly(p1, p2)

numCoefs1 = size(p1{1}, 1);
numCoefs2 = size(p2{1}, 1);

pProd = { reshape(p2{1}*(p1{1}'), [], 1), false(numCoefs1*numCoefs2, size(p1{2}, 2)) };

idxCounter = 0;
for c1 = 1:numCoefs1
    for c2 = 1:numCoefs2
        idxCounter = idxCounter+1;
        pProd{2}(idxCounter, :) = p1{2}(c1, :) | p2{2}(c2, :);
    end
end

pProd = reducePoly(pProd);