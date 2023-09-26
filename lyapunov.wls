AA = (DiagonalMatrix[{-m[i], - d[i], -d[i], - d[i]}] + DiagonalMatrix[{d[i],d[i], d[i]},-1]);
AA // MatrixForm
A = ArrayFlatten[{{AA /. i -> 1, 0}, {0, AA /. i -> 2}}];
A // MatrixForm
B = (DiagonalMatrix[{1/(2*m[1]), 1/(2*d[1]), 1/(2*d[1]), 1/(2*d[1]),1/(2*m[2]), 1/(2*d[2]), 1/(2*d[2]), 1/(2*d[2])}]);
Eigenvalues[B]
CC = B.A + Transpose[A].B;
CC // MatrixForm // FullSimplify
Do[Print[Det[CC[[1 ;; i, 1 ;; i]]] // FullSimplify], {i, 8}]
x = {{NN[1], x1[1], x1[2], x1[3], NN[2], x2[1], x2[2], x2[3]}}
gx = {{r[1]*e[1]/d[1]*x1[3]*Exp[-\[CapitalDelta][1]*(a[1]*f[1] + b[2]*f[2])]}, {0}, {0}, {0}, {r[2]*e[2]/d[2]*x2[3]*Exp[-\[CapitalDelta][2]*(a[2]*f[2] + b[1]*f[1])]}, {0}, {0}, {0}}
DD = 2*x.B.gx;
DD // MatrixForm // FullSimplify
VV = x.CC.Transpose[x]
VV // MatrixForm // FullSimplify

