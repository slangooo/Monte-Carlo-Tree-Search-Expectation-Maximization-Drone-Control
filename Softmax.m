function res = Softmax (X)
A = exp(X)';
res = (A./sum(A,2))';
res(isnan(res))=1;