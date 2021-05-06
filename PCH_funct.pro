; NAME:
;
; LOCATION:
;
; PURPOSE:
;	This is the probability function that models the pixel's histogram. This function will be called in CG calculation (PCH method)
;
; INPUTS:
;	X		- x values (2D array)
;	A		- Initial guess for parameter
;	F		- function of
;	pder	- partial derivative with respect to each param
; KEYWORD PARAMETERS:
;	A		- dimension equals to number of parameters. In this case, there are 2: random noise (u_n) and quanta exposure (H)
;	A[0]: u_n , A[1]: H
;
; OUTPUTS:
;	value of y at the intial guess A. CURVEFIT will compare this with actual y values and find optimal parameters
;	partial derivatives of each param. partial derivative should be an array whose dimension equals to X's.
;	A larger array should contain each pder, see example CUVEFIT on harrisgeospatial
;
; EXAMPLE:
; IDL> X = [ array of X locations ]
; IDL> A = [readnoise_estimate,average_estimate,offset_estimate]
; IDL> PCH_funct, X, A, F, pder
; IDL> plot,X,F
;
; MODIFICATION HISTORY:
;
;	Written by: Long Nguyen, CfD, November 19 , 2020
;----------------------------------------------------------------------------------------------------------------------------------
	;A[0]: RN ; A[1]: H;  A[2]: offset

pro PCH_funct, X, A, F, pder
	element_X = X[n_elements(X)-1] - X[0]
	k = double(findgen(element_X + 20)) ; number of repetition, theoretically it should go to inf
	F = double(make_array (n_elements(X)))
	A = double(A)
	;poisson_theory = double(make_array (n_elements(X)))
	gauss_theory = double(make_array (n_elements(X)))
	poisson = exp(-A[1])*A[1]^k/factorial(k)
	for i= 0, n_elements(X)-1 do begin
		norm_const = 1/(A[0]*sqrt(2*!pi))
		in_exp = -(X[i]+A[2]-k)^2/(2*A[0]^2)
		gauss = exp(in_exp)
		gauss_theory[i] = total(gauss)
		;poisson = exp(-A[1])*A[1]^k/factorial(k)
		F[i] =total(norm_const*gauss*poisson) ; need to make all data the same data type!!!
	endfor

;------------------partial derivatives---------------------------------;

	pder_un = double(make_array (n_elements(X))) ; pder[0] is derivative of u_n, pder[1] is derivative of H
	pder_H  = double(make_array (n_elements(X)))
	pder_kos  = double(make_array (n_elements(X)))

	for i=0, n_elements(X)-1 do begin
		norm_const = 1/(A[0]^2*sqrt(2*!pi))
		in_exp = -(X[i]+A[2]-k)^2/(2*A[0]^2)
		gauss = exp(in_exp)
		add_factor = ((X[i]+A[2]-k)/A[0])^2-1
		;poissson as above
		pder_un[i] = total(norm_const*add_factor*gauss*poisson)
	endfor

	for i=0, n_elements(X)-1 do begin
		norm_const = 1/(A[0]*sqrt(2*!pi))
		in_exp = -(X[i]+A[2]-k)^2/(2*A[0]^2)
		gauss = exp(in_exp)
		add_factor = (k/A[1] -1)
		pder_H[i] = total(norm_const*add_factor*gauss*poisson)
	endfor

	for i=0, n_elements(X)-1 do begin
		norm_const = 1/(A[0]*sqrt(2*!pi))
		in_exp = -(X[i]+A[2]-k)^2/(2*A[0]^2)
		gauss = exp(in_exp)
		add_factor =  -(X[i]+A[2]-k)/A[0]^2
		pder_kos[i] = total(norm_const*add_factor*gauss*poisson)
	endfor

	pder = [[pder_un], [pder_H], [pder_kos]]
	;stop
end