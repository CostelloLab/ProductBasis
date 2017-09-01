% polishEigs() improves an eigenvalue/column eigenvectors estimate, and
% returns error bounds for the eigenvectors using the method of Dongarra
% et al. (1983) -- in particular Section 3.
% 
% All eigenvectors are assumed to have the same eigenvalue.  The set of
% eigenvectors must be complete (span the subspace of the eigenvalue).
% 
% All variable symbols are taken from their paper.  Note:  our 'm' is only
% the off-diagonal part of the m-matrix in the paper; the diagonal part is
% 'lambda'; however, 'mu' is the full correction matrix (diagonal and
% off-diagonal elements both).

function [ lambda, x, lambdaErr, xErr ] = polishEigs(A, lambda, x, maxIterations)

n = length(A);
if ~exist('x', 'var')
    x = [];
end

if isempty(x)
    [ U, T ] = schur(A, 'complex');
    T_lI = T - lambda*eye(n);
    
    if size(x, 2) == 0
        epsT = eps(max(max(A)))^.5;
    else
        srtT = sort(abs(diag(T_lI)));
        epsT = srtT(size(x, 2));
    end
    l = [ 0, find(abs(diag(T_lI)') <= epsT) ];
    k = length(l)-1;
    m = zeros(k);
    
    z = zeros(n, k);
    lowerIdx = zeros(1, 0);
    zSolve = T_lI;
    for loopZ = 1:k
        lowerIdx = [ l(loopZ+1)-1:-1:l(loopZ)+1, lowerIdx ];
        z(l(loopZ+1), loopZ) = 1;
        for loopIdx = lowerIdx
            z(loopIdx, loopZ) = -zSolve(loopIdx, :)*z(:, loopZ) / zSolve(loopIdx, loopIdx);
        end
        m(1:loopZ-1, loopZ) = T_lI(l(2:loopZ), :)*z(:, loopZ);
        
        zSolve = zSolve - z(:, loopZ)*T_lI(l(1+loopZ), :);
    end
    
    x = U*z;
else
    k = size(z, 2);
    m = zeros(k);
end

if ~exist('maxIterations', 'var')
    maxIterations = 100;
end

[ ~, maxIDs ] = max(abs(x), [], 1);
xNorm = x((0:k-1)*n+maxIDs);
x = x ./ repmat(xNorm, n, 1);
m = m .* repmat(xNorm', 1, k) ./ repmat(xNorm, k, 1);


    % find a large element of each 'x' for insertion into the respective
    % column of 'B'

xsToSort = x;
s = zeros(1, k);
remainingXs = 1:k;
%twiddleIdx = false(1, k*n);
for loopX = 1:k
    [ xMaxs, maxIDs ] = max(abs(xsToSort(:, remainingXs)));
    [ ~, whichRemainingX ] = max(xMaxs);
    whichX = remainingXs(whichRemainingX);
    s(whichX) = maxIDs(whichRemainingX);
    remainingXs(whichRemainingX) = [];
    
    xScale = xsToSort(s(whichX), remainingXs)/xsToSort(s(whichX), whichX);
    xsToSort(:, remainingXs) = xsToSort(:, remainingXs) - xsToSort(:, whichX) * xScale;
    
    %twiddleIdx((loopX-1)*n+(1:n)) = ((1:n) ~= s(whichX));
end
%sIdx = (0:k-1)*n + s;
notS = 1:n;
notS(s) = [];


    % run the approximation algorithm

tol = eps(min(abs(lambda), min(min(abs(x))))) / 2;

B = A - lambda*eye(n);
B(:, s) = -x;
X = inv(B);

deltaTwiddle = zeros(n, k);
yTwiddle = deltaTwiddle;

r = lambda*x - (A*x - x*m);
epskap = norm(X, Inf)^2 * max(max(abs(r)));
if epskap >= 1/8
    error('root not convergent')
end
gamma = 4*epskap * sqrt(2 / (1 - 4*epskap + sqrt(1 - 8*epskap)));

y = zeros(n, k);
delta = y;
for iteration = 1:maxIterations
    y0 = y;
    deltaTwiddle = zeros(n, k);
    for loopX = 1:k
        delta(:, loopX) = B \ r(:, loopX) + deltaTwiddle*m(:, loopX);     % ? - sign difference between paper in last term
        deltaTwiddle(notS, loopX) = delta(notS, loopX);
        y(:, loopX) = y(:, loopX) + delta(:, loopX);
    end
    
    yTwiddle(notS, :) = y(notS, :);
    
    maxErr = gamma * sum(max(abs(delta), [], 1)) / (1-gamma);
    if maxErr < tol
        break
    end
    
    r = (r - B*(delta - deltaTwiddle*m)) + deltaTwiddle*y0(s, :) + yTwiddle*delta(s, :);
end

lambda = lambda*eye(k) + m + y(s, :);
lambdaErr = maxErr + eps(lambda);
x = x + yTwiddle;
xErr = maxErr + eps(x);
