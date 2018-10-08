/*
This file contains code that tranlates into Stata/Mata language the original R code presented in the appendix of the following paper:

Rizzi, S., Halekoh, U., Thinggaard, M., Engholm, G., Christensen, N., Johannesen, T. B., & Lindahl-Jacobsen, R. (2018). **How to estimate mortality trends from grouped vital statistics**. _International journal of epidemiology._ [[PubMed]](https://www.ncbi.nlm.nih.gov/pubmed/30256946)
*/

version 14.2

mata:
mata clear 
mata set matastrict on 

real matrix function diff_R(real matrix X, real scalar diff){
	real matrix Y
	real vector B
	real scalar i, d

	Y = J(rows(X)-diff, cols(X), .) 
	
	for (i=1; i<=cols(X); i++) {
		for (d=0; d<=diff; d++) {
			if (d == 0) {
				B = X[, i]
			}
			else {
				B = B[(1+1)..rows(B)] - B[1..(rows(B)-1)]
			}
		}
		Y[., i] = B
	}
	return(Y)
}

real vector function MortSmooth_tpower_R(real vector x, real vector t, real scalar p) {
	return((x :- t):^p :* (x :> t))
}

real matrix function MortSmooth_bbase_R(real vector x, real scalar xl, real scalar xr, real scalar ndx, real scalar deg) {
	real scalar dx, i
	real vector knots
	real matrix P, D, B
	
	dx = (xr - xl)/ndx
	knots = range(xl - deg * dx, xr + deg * dx, dx)
	P = J(length(x), length(knots), .)
	for (i=1; i<=length(knots); i++) {
		P[., i] = MortSmooth_tpower_R(x, knots[i], deg)
	}
	D = diff_R(I(cols(P)), deg + 1):/(gamma(deg + 1) * dx^deg)
    B = (-1)^(deg + 1) :* P * D'
	return(B)
}

function pclm2D(real vector y, real matrix C, real matrix B, real matrix P, b, V) {
//   # Fit a 2D PCLM (estimate A in E(y) = C %*% exp(B %*% A))
//   # y = the matrix of observed counts of length IxN in vector format
//   # C = the composition matrix 
//   # B = the two-dimensional B-spline basis
//   # P = penalty matrix

	real scalar nx, ly, it, da, trace, dev, psi2, aic, bic
	real vector Astart, A, A0, eta, mu, w, gam
	real matrix Q, z, QWQ, H, H0, H1
	
//   # Preparations
	nx = cols(B)
	ly = length(y)
	it = 0
	Astart = log(sum(y) / ly)
	A = J(nx, 1, Astart)
	da = 1

//  # Iterations
	while(abs(da) > 1e-6) {
		it++
		A0 = A
		eta = B * A
		gam = exp(eta)
		mu = C * gam
		w = mu
		
		Q = (C :* ((1 :/ mu) * gam')) * B
		z = y :- mu :+ w :* Q * A
		QWQ = Q' * (w :* Q)
		A = cholsolve(QWQ :+ P, Q' * z)
		da = max(A - A0)
		printf("Iteration %f: diff = %f\n", it, da)
		displayflush()
	}
	
	H = cholsolve(QWQ :+ P, QWQ)
	H0 = invsym(QWQ  :+ P)
	H1 = H0 * QWQ  * H0

	b = A
	V = H1
	
	trace = trace(H)
	dev = 2 :* sum(select(y, y:>0) :* log(select(y, y:>0) :/ select(mu, y:>0)))
	psi2 = dev / (length(y)-trace)
	aic = dev + 2 * trace
	bic = dev + log(length(y)) * trace
	printf("\r Trace = %f, Dev = %f, Psi2 = %f, AIC = %f, BIC = %f\n", trace, dev, psi2, aic, bic)
}

end
