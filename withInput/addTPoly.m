% addPoly() returns the sum of two polynomials p1 and p2

function pSum = addPoly(p1, p2)

pSum = reducePolyT({ [ p1{1}; p2{1} ], [ p1{2}; p2{2} ], [ p1{3}; p2{3} ], [ p1{4}; p2{4} ] });

end