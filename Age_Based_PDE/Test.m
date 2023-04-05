S = [ 1 2 3 4 5]'
y = [0 2 4 5 7]'
for k = 1:10
  betaHat = repmat(beta', 1, 5);
end
betaHat = betaHat';
betaHatSUM = sum(betaHat,2)

betaHatSUM(3) = 0
 S.*y.*(betaHatSUM)