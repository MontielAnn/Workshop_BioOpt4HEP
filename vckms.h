/* Rutinas para calcular los elementos de la matriz CKM a partir
 de una matriz de masa

*/
#ifndef VCKMS2_H
#define VCKMS2_H

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define Pi 3.14159265358979323846264338327

/*Sean las masas de los quarks*/
long double A, C, EE, lambda1, lambda2, lambda3;
/*
long double mu = 2.16;
long double md = 4.67;
long double ms = 93;
long double mc = 1270;
long double mb = 4180;
long double mt = 172760;
*/
#define 	mu 	 1.23
#define 	md  	2.67
#define 	ms  	53.16
#define 	mc  	620
#define 	mb  	2839
#define 	mt  	168260

#define minAu mc
#define maxAu mt
#define minEu mu
#define maxEu mc
#define minAd ms
#define maxAd mb
#define minEd md
#define maxEd ms
#define minPhi1 0
#define maxPhi1 2*Pi
#define minPhi2 0
#define maxPhi2 2*Pi

/*Valores teóricos de los elementos de la matriz CKM*/
#define 	vcb_th		0.04053
#define 	vus_th		0.22650
#define 	vub_th		0.00361
#define 	jarslkog_th	0.00003

/*donde la incertidumbre de cada uno es*/
#define 	sigma_vcb	0.00061
#define 	sigma_vus	0.00048
#define 	sigma_vub	0.00009
#define 	sigma_jars	0.0000009

//__BEGIN_DECLS

long double BB(long double, long double, long double, long double, long double);
long double DD(long double, long double, long double, long double, long double);
long double O11(long double, long double, long double, long double, long double);
long double O22(long double, long double, long double, long double, long double);
long double O33(long double, long double, long double, long double, long double);
long double O12(long double, long double, long double, long double, long double);
long double O21(long double, long double, long double, long double, long double);
long double O32(long double, long double, long double, long double, long double);
long double O23(long double, long double, long double, long double, long double);
long double O13(long double, long double, long double, long double, long double);
long double O31(long double, long double, long double, long double, long double);
long double complex CKM11(long double, long double, long double, long double, long double, long double);
long double complex CKM22(long double, long double, long double, long double, long double, long double);
long double complex CKM33(long double, long double, long double, long double, long double, long double);
long double complex CKM12(long double, long double, long double, long double, long double, long double);
long double complex CKM21(long double, long double, long double, long double, long double, long double);
long double complex CKM13(long double, long double, long double, long double, long double, long double);
long double complex CKM31(long double, long double, long double, long double, long double, long double);
long double complex CKM32(long double, long double, long double, long double, long double, long double);
long double complex CKM23(long double, long double, long double, long double, long double, long double);
long double JJ(long double, long double, long double, long double, long double, long double);
long double chi2CKM23(long double, long double, long double, long double, long double, long double);
long double chi2CKM12(long double, long double, long double, long double, long double, long double);
long double chi2CKM13(long double, long double, long double, long double, long double, long double);
long double chi2JJ(long double, long double, long double, long double, long double, long double);
long double chi2(long double, long double, long double, long double, long double, long double);

//__END_DECLS

#endif
