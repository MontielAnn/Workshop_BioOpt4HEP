#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include <unistd.h>
#include "vckms.h"

#define PI 3.14159265358979323846264338327
#define N 10000
#define dim 6
#define kmax 60
#define c1min 0.5
#define c1max 2
#define c2min 0
#define c2max 2.5
#define w 0.9
#define c3 0.7

long double Par[8]; //={N,dim,kmax,c1max,c2max,c1min,c2min,c3};
long double lb[dim]; //={mc, ms, mu, md, 0, 0};
long double ub[dim]; //={mt, mb, mc, ms, 2*PI, 2*PI};
long double x[N][dim];
long double p[N][dim];
long double pi[N];
long double vi[N];
long double v[N][dim];
long double chi[N][dim];
long double xi[dim];
long double fx[N];
long double fp[N];
long double gfit;
long double gfit2;
long double gfitkmas1;
long double g[dim];
long double c1;
long double c2;
long double phi;
long double ch;
int ind;
int k;
long double cad[kmax];


void PSO(long double Par[], long double lb[], long double ub[]){
    int i,j;
//Semilla de números aleatorios
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand48(t1.tv_usec * t1.tv_sec);
    
//Inicialización de las paetículas y su velocidad
    for (i=0;i<N;i++){
        for(j=0;j<dim;j++){
            x[i][j]=((long double) drand48())*(ub[j]-lb[j])+lb[j];
            v[i][j]=0;
        }
    }
    
//Evaluación de las partículas iniciales en la función objetivo
    for(i=0;i<N;i++){
        fx[i]=chi2(x[i][0],x[i][1],x[i][2],x[i][3],x[i][4],x[i][5]);
    }
    
//Registro de la mejor partícula global y las mejores locales
    gfit=fx[0];
    ind=0;
    for(i=1;i<N;i++){
        if(gfit>fx[i]){
            gfit=fx[i];
            ind=i;
        }
    }
    for(j=0;j<dim;j++){
        g[j]=x[ind][j];
    }
    
    for(i=0;i<N;i++){
        fp[i]=fx[i];
        for(j=0;j<dim;j++){
            p[i][j]=x[i][j];
        }
    }

//PROCESO ITERATIVO
    k=0;
    while(k<kmax){
        sleep(((long double) drand48()));
        
    //Calculo de la nueva velocidad para cada partícula
        for(i=0;i<N;i++){
            c1=((c1min+c1max)/2)+((c1max-c1min)/2);
            c2=((c2min+c2max)/2)+((c2max-c2min)/2);
            phi=((c1*(long double) drand48())+(c2*(long double) drand48()));
            ch=2*c3/(fabsl(2-phi-sqrt(phi*(phi))));
            for(j=0;j<dim;j++){
                v[i][j]=ch*((w*v[i][j])+(c1*((long double) drand48())*(p[i][j]-x[i][j]))+(c2*((long double) drand48())*(g[j]-x[i][j])));
            }
        }
        
    //Cálculo de la nueva posición de cada partícula
        for(i=0;i<N;i++){
            for(j=0;j<dim;j++){
                x[i][j]=x[i][j]+v[i][j];
            }
        }
        
    //Verificar que las partículas no se salgan de los límites
        for(i=0;i<N;i++){
            for(j=0;j<dim;j++){
                if(x[i][j]<lb[j]){
                    x[i][j]=lb[j]+(0.25*((long double) drand48())*lb[j]);
                }
                else if(x[i][j]>ub[j]){
                    x[i][j]=ub[j]-(0.25*(((long double) drand48())*ub[j]));
                }
            }
        }
        
    //Evaluación de las nuevas partículas en la función objetivo
        for(i=0;i<N;i++){
            fx[i]=chi2(x[i][0],x[i][1],x[i][2],x[i][3],x[i][4],x[i][5]);
        }
    //Registro de la mejor partícula global y las mejores locales
        gfitkmas1=fx[0];
        for (i=1;i<N;i++){
            if (gfitkmas1>fx[i]){
                gfitkmas1=fx[i];
                ind=i;
            }
            if (gfit>gfitkmas1){
                gfit=gfitkmas1;
                ind=i;
            }
            //for(j=0;j<dim;j++){
                //g[j]=x[ind][j];
            //}
        }
        
        if(x[ind][1]!=x[ind][1]){
            k=kmax;
        }
        else{
            for(i=1;i<N;i++){
                for(j=0;j<dim;j++){
                    g[j]=x[ind][j];
                }
            }
            printf("\n%d %2.15Lf ",k,gfit);
            for(j=0;j<dim;j++){
                printf("%.9Lf, ",g[j]);
            }
            printf("\n");    
        }
        
    
        
    //Actualiza el fitness de las mejores partículas locales
        for (i=0;i<N;i++){
            if(fx[i]<fp[i]){
                fp[i]=fx[i];
                for (j=0;j<dim;j++){
                    p[i][j]=x[i][j];
                }
            }
        }
        k++;
        

    }
    if(gfit<1){
        FILE* dat1;
        dat1=fopen ("EXITOSAS.txt", "a");
        for(i=0;i<dim;i++){
            fprintf(dat1, "%.9Lf ",g[i]);    
        }
        fprintf(dat1, "%2.15Lf\n",gfit);
        fclose(dat1);
    }
    
    return;
}

int main(){
    int i,j;
    long double Par[8]={N,dim,kmax,c1max,c2max,c1min,c2min,c3};
    long double lb[dim]={mc, ms, mu, md, 0, 0};
    long double ub[dim]={mt, mb, mc, ms, 2*PI, 2*PI};
    for(int m=0;m<1;m++){
        printf("\n\nSolución %d\n",m+1);
        PSO(Par,lb,ub);
        printf("\ngfit: %2.15Lf",gfit);
        printf("\n g: ");
        for(j=0;j<dim;j++){
            printf("%Lf, ",g[j]);
        }   
    }
    return 0;
}
