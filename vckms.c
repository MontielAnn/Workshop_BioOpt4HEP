#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vckms.h"


/*Se definen los parámetros de la siguiente manera: */
long double BB(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double bb;
        bb=sqrt(((A-lambda1)*(A-lambda2)*(lambda3-A))/(A-EE));
return bb;
}
long double DD(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double dd;
        dd=sqrt(((EE-lambda1)*(lambda2-EE)*(lambda3-EE))/(A-EE));
return dd;
}

/* los elementos diagonales */
long double O11(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o11;
        o11=1/(sqrt(1+((EE-lambda1)*(A-EE)/((lambda2-EE)*(lambda3-EE)))+((EE-lambda1)*(A-lambda2)*(lambda3-A)/((lambda2-EE)*(lambda3-EE)*(A-lambda1)))));
return o11;
}
long double O22(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o22;
        o22=1/(sqrt(1+((EE-lambda1)*(lambda3-EE)/((lambda2-EE)*(A-EE)))+((A-lambda1)*(lambda3-A)/((A-EE)*(A-lambda2)))));
return o22;
}
long double O33(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o33;
        o33=1/(sqrt(1+((lambda3-A)*(A-EE)/((A-lambda2)*(A-lambda1)))+((EE-lambda1)*(lambda2-EE)*(lambda3-A)/((A-lambda2)*(lambda3-EE)*(A-lambda1)))));
return o33;
}
/*los elementos fuera de la diagonal, donde todos dependen de A, EE, lambda1, lambda2, lambda3*/
long double O12(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o12;
        o12=fabsl(DD(A, EE, lambda1, lambda2, lambda3))/(lambda2-EE)*O22(A, EE, lambda1, lambda2, lambda3);
return o12;
}
long double O21(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o21;
        o21=(-(EE-lambda1)/(fabsl(DD(A, EE, lambda1, lambda2, lambda3))))*O11(A, EE, lambda1, lambda2, lambda3);
return o21;
}
long double O32(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o32;
        o32=(-(fabsl(BB(A, EE, lambda1, lambda2, lambda3)/(A-lambda2))))*O22(A, EE, lambda1, lambda2, lambda3);
return o32;
}
long double O23(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o23;
        o23=((lambda3-A)/fabsl(BB(A, EE, lambda1, lambda2, lambda3)))*O33(A, EE, lambda1, lambda2, lambda3);
return o23;
}
long double O13(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o13;
        o13=fabsl(DD(A, EE, lambda1, lambda2, lambda3)/(lambda3-EE))*((lambda3-A)/(fabsl(BB(A, EE, lambda1, lambda2, lambda3))))*O33(A, EE, lambda1, lambda2, lambda3);
return o13;
}
long double O31(long double A, long double EE, long double lambda1, long double lambda2, long double lambda3){
        long double o31b;
        o31b=(fabsl(BB(A, EE, lambda1, lambda2, lambda3))/(A-lambda1))*((EE-lambda1)/fabsl(DD(A, EE, lambda1, lambda2, lambda3)))*O11(A, EE, lambda1, lambda2, lambda3);
return o31b;
}
/* Los elementos de la VCKM */
long double complex CKM11(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm11;
        ckm11=(O11(Au, EEu, mu, mc, mt)*O11(Ad, EEd, md, ms, mb)) + (O21(Au, EEu, mu, mc, mt)*O21(Ad, EEd, md, ms, mb)*cexp(I*phi1)) + (O31(Au, EEu, mu, mc, mt)*O31(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm11;
}
long double complex CKM22(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm22;
        ckm22=(O12(Au, EEu, mu, mc, mt)*O12(Ad, EEd, md, ms, mb))+(O22(Au, EEu, mu, mc, mt)*O22(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O32(Au, EEu, mu, mc, mt)*O32(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm22;
}
long double complex CKM33(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm33;
        ckm33=(O13(Au, EEu, mu, mc, mt)*O13(Ad, EEd, md, ms, mb))+(O23(Au, EEu, mu, mc, mt)*O23(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O33(Au, EEu, mu, mc, mt)*O33(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm33;
}
long double complex CKM12(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm12;
        ckm12=(O11(Au, EEu, mu, mc, mt)*O12(Ad, EEd, md, ms, mb))+(O21(Au, EEu, mu, mc, mt)*O22(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O31(Au, EEu, mu, mc, mt)*O32(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm12;
}
long double complex CKM21(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm21;
        ckm21=(O12(Au, EEu, mu, mc, mt)*O11(Ad, EEd, md, ms, mb))+(O22(Au, EEu, mu, mc, mt)*O21(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O32(Au, EEu, mu, mc, mt)*O31(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm21;
}
long double complex CKM13(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm13;
        ckm13=(O11(Au, EEu, mu, mc, mt)*O13(Ad, EEd, md, ms, mb))+(O21(Au, EEu, mu, mc, mt)*O23(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O31(Au, EEu, mu, mc, mt)*O33(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm13;
}
long double complex CKM31(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm31;
        ckm31=(O13(Au, EEu, mu, mc, mt)*O11(Ad, EEd, md, ms, mb))+(O23(Au, EEu, mu, mc, mt)*O21(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O33(Au, EEu, mu, mc, mt)*O31(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm31;
}
long double complex CKM32(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm32;
        ckm32=(O13(Au, EEu, mu, mc, mt)*O12(Ad, EEd, md, ms, mb))+(O23(Au, EEu, mu, mc, mt)*O22(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O33(Au, EEu, mu, mc, mt)*O32(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm32;
}
long double complex CKM23(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double complex ckm23;
        ckm23=(O12(Au, EEu, mu, mc, mt)*O13(Ad, EEd, md, ms, mb))+(O22(Au, EEu, mu, mc, mt)*O23(Ad, EEd, md, ms, mb)*cexp(I*phi1))+(O32(Au, EEu, mu, mc, mt)*O33(Ad, EEd, md, ms, mb)*cexp(I*(phi1+phi2)));
return ckm23;
}
/*El invariante de Jarlskog*/
long double JJ(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double jj;
        jj=cimagf(CKM12(Au, Ad, EEu, EEd, phi1, phi2)*CKM23(Au, Ad, EEu, EEd, phi1, phi2)*conjf(CKM13(Au, Ad, EEu, EEd, phi1, phi2))*conjf(CKM22(Au, Ad, EEu, EEd, phi1, phi2)));
return jj;
}

long double chi2CKM23(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double chi23;
        chi23=pow((vcb_th -cabsl(CKM23(Au, Ad, EEu, EEd, phi1, phi2))),2)/pow(sigma_vcb,2);
return chi23;
}
long double chi2CKM12(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double chi12;
        chi12=pow(vus_th -cabsl(CKM12(Au, Ad, EEu, EEd, phi1, phi2)),2)/pow(sigma_vus,2);
return chi12;
}
long double chi2CKM13(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double chi13;
        chi13=(pow(vub_th -cabsl(CKM13(Au, Ad, EEu, EEd, phi1, phi2)),2))/(pow(sigma_vub,2));
return chi13;
}
long double chi2JJ(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double chijj;
        chijj=(pow(jarslkog_th -JJ(Au, Ad, EEu, EEd, phi1, phi2),2))/pow(sigma_jars,2);
return chijj;
}
long double chi2(long double Au, long double Ad, long double EEu, long double EEd, long double phi1, long double phi2){
        long double chi2_;
        chi2_=chi2CKM23(Au, Ad, EEu, EEd, phi1, phi2)+chi2CKM12(Au, Ad, EEu, EEd, phi1, phi2)+chi2CKM13(Au, Ad, EEu, EEd, phi1, phi2)+chi2JJ(Au, Ad, EEu, EEd, phi1, phi2);
return chi2_;
}


