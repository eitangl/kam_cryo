function Q = randR(n)

tmp = randn(n);
[Q,R] = qr(tmp);
D = diag(sign(diag(R)));
Q = Q*D;