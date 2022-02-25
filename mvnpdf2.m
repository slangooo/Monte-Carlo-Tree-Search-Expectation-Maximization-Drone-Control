function y = mvnpdf2(X, Mu, Sigma)
[n,d] = size(X);
X0 = X - Mu;
[R,err] = cholcov(Sigma,0);
if err ~= 0
    error(message('stats:mvnpdf:BadMatrixSigma'));
end
% Create array of standardized data, and compute log(sqrt(det(Sigma)))
xRinv = X0 / R;
logSqrtDetSigma = sum(log(diag(R)));
quadform = sum(xRinv.^2, 2);
y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);
