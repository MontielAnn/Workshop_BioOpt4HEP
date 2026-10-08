//
// Created by Sippawit Thammawiset on 29/1/2024 AD.
// https://github.com/sippathamm/BenchmarkFunction
//

// ====================================================================
//  BENCHMARK FUNCTIONS
// ====================================================================
// Contiene 10 funciones benchmark clásicas para probar
// ====================================================================

#ifndef BENCHMARK_FUNCTION_H
#define BENCHMARK_FUNCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


#define SPHERE            0
#define SCHWEFEL_S_2_22   1
#define SCHWEFEL_S_1_20   2
#define ROSENBROCK        3
#define STEP              4
#define QUARTIC_NOISE     5
#define SCHWEFEL_S_2_26   6
#define RASTRIGIN         7
#define ACKLEY            8
#define GRIEWANK          9


// Generar número aleatorio en [LowerBound, UpperBound]
static inline double GenerateRandom(double LowerBound, double UpperBound) {
    return LowerBound + ((double)rand() / RAND_MAX) * (UpperBound - LowerBound);
}


// ====================================================================
// FUNCIONES BENCHMARK
// ====================================================================

// --- Sphere ---
// Mínimo global: f(0,...,0) = 0
static inline double Function_Sphere(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable; i++) {
        Sum += Position[i] * Position[i];
    }
    return Sum;
}

// --- Schwefel 2.22 ---
// Mínimo global: f(0,...,0) = 0
static inline double Function_Schwefel_s_2_22(const double *Position, int NVariable) {
    double Term1 = 0.0;
    double Term2 = 1.0;
    for (int i = 0; i < NVariable; i++) {
        Term1 += fabs(Position[i]);
        Term2 *= fabs(Position[i]);
    }
    return Term1 + Term2;
}

// --- Schwefel 1.20 ---
// Mínimo global: f(0,...,0) = 0
static inline double Function_Schwefel_s_1_20(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable; i++) {
        double InnerSum = 0.0;
        for (int j = 0; j <= i; j++) {
            InnerSum += Position[j];
        }
        Sum += InnerSum * InnerSum;
    }
    return Sum;
}

// --- Rosenbrock ---
// Mínimo global: f(1,...,1) = 0
static inline double Function_Rosenbrock(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable - 1; i++) {
        double Term1 = 100.0 * (Position[i + 1] - Position[i] * Position[i])
                             * (Position[i + 1] - Position[i] * Position[i]);
        double Term2 = (Position[i] - 1.0) * (Position[i] - 1.0);
        Sum += Term1 + Term2;
    }
    return Sum;
}

// --- Step ---
// Mínimo global: f(-0.5,...,-0.5) = 0
static inline double Function_Step(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable; i++) {
        double val = floor(Position[i] + 0.5);
        Sum += val * val;
    }
    return Sum;
}

// --- Quartic con ruido ---
//El mínimo teórico está en(0,...,0), pero el ruido aleatorio lo perturba.
static inline double Function_QuarticNoise(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable; i++) {
        double p2 = Position[i] * Position[i];
        Sum += (double)i * p2 * p2;    // Bug original preservado
    }
    return Sum + GenerateRandom(0.0, 1.0);
}

// --- Schwefel 2.26 ---
// Mínimo global: f(420.9687,...,420.9687) ≈ -418.9829 * N
static inline double Function_Schwefel_s_2_26(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable; i++) {
        Sum += Position[i] * sin(sqrt(fabs(Position[i])));
    }
    return -Sum;
}

// --- Rastrigin ---
// Mínimo global: f(0,...,0) = 0
static inline double Function_Rastrigin(const double *Position, int NVariable) {
    double Sum = 0.0;
    for (int i = 0; i < NVariable; i++) {
        Sum += Position[i] * Position[i]
             - 10.0 * cos(2.0 * M_PI * Position[i])
             + 10.0;
    }
    return Sum;
}

// --- Ackley ---
// Mínimo global: f(0,...,0) = 0
static inline double Function_Ackley(const double *Position, int NVariable) {
    double SumSquare = 0.0;
    double SumCosine = 0.0;
    for (int i = 0; i < NVariable; i++) {
        SumSquare += Position[i] * Position[i];
        SumCosine += cos(2.0 * M_PI * Position[i]);
    }
    double Term1 = -20.0 * exp(-0.2 * sqrt(SumSquare / NVariable));
    double Term2 = -exp(SumCosine / NVariable);
    return Term1 + Term2 + 20.0 + exp(1.0);
}

// --- Griewank ---
// Mínimo global: f(0,...,0) = 0
static inline double Function_Griewank(const double *Position, int NVariable) {
    double SumSquare = 0.0;
    double ProductCosine = 1.0;
    for (int i = 0; i < NVariable; i++) {
        SumSquare += Position[i] * Position[i];
        ProductCosine *= cos(Position[i] / sqrt((double)(i + 1)));
    }
    return SumSquare / 4000.0 - ProductCosine + 1.0;
}


// ====================================================================
// selecciona la función benchmark por nombre
// ====================================================================
static inline double BenchmarkFunction(int FUNCTION, const double *Position, int NVariable) {
    switch (FUNCTION) {
        case SPHERE:           return Function_Sphere(Position, NVariable);
        case SCHWEFEL_S_2_22:  return Function_Schwefel_s_2_22(Position, NVariable);
        case SCHWEFEL_S_1_20:  return Function_Schwefel_s_1_20(Position, NVariable);
        case ROSENBROCK:       return Function_Rosenbrock(Position, NVariable);
        case STEP:             return Function_Step(Position, NVariable);
        case QUARTIC_NOISE:    return Function_QuarticNoise(Position, NVariable);
        case SCHWEFEL_S_2_26:  return Function_Schwefel_s_2_26(Position, NVariable);
        case RASTRIGIN:        return Function_Rastrigin(Position, NVariable);
        case ACKLEY:           return Function_Ackley(Position, NVariable);
        case GRIEWANK:         return Function_Griewank(Position, NVariable);
        default:               return Function_Sphere(Position, NVariable);
    }
}


// ====================================================================
// PROPIEDADES: Dimensión, límites inferior y superior de cada función
// --------------------------------------------------------------------
// ====================================================================

static inline void Property_Sphere(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -100.0;
        (*UpperBound)[i] =  100.0;
    }
}

static inline void Property_Schwefel_s_2_22(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -10.0;
        (*UpperBound)[i] =  10.0;
    }
}

static inline void Property_Schwefel_s_1_20(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -100.0;
        (*UpperBound)[i] =  100.0;
    }
}

static inline void Property_Rosenbrock(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -30.0;
        (*UpperBound)[i] =  30.0;
    }
}

static inline void Property_Step(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -100.0;
        (*UpperBound)[i] =  100.0;
    }
}

static inline void Property_QuarticNoise(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -1.28;
        (*UpperBound)[i] =  1.28;
    }
}

static inline void Property_Schwefel_s_2_26(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -500.0;
        (*UpperBound)[i] =  500.0;
    }
}

static inline void Property_Rastrigin(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -5.12;
        (*UpperBound)[i] =  5.12;
    }
}

static inline void Property_Ackley(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -32.0;
        (*UpperBound)[i] =  32.0;
    }
}

static inline void Property_Griewank(int *NVariable, double **LowerBound, double **UpperBound) {
    *NVariable = 30;
    *LowerBound = (double *)malloc(sizeof(double) * (*NVariable));
    *UpperBound = (double *)malloc(sizeof(double) * (*NVariable));
    for (int i = 0; i < *NVariable; i++) {
        (*LowerBound)[i] = -600.0;
        (*UpperBound)[i] =  600.0;
    }
}


// ====================================================================
// Selecciona las propiedades por nombre de función
// ====================================================================
static inline void BenchmarkProperty(int FUNCTION,
                                     int *NVariable,
                                     double **LowerBound,
                                     double **UpperBound) {
    switch (FUNCTION) {
        case SPHERE:           Property_Sphere(NVariable, LowerBound, UpperBound);           break;
        case SCHWEFEL_S_2_22:  Property_Schwefel_s_2_22(NVariable, LowerBound, UpperBound);  break;
        case SCHWEFEL_S_1_20:  Property_Schwefel_s_1_20(NVariable, LowerBound, UpperBound);  break;
        case ROSENBROCK:       Property_Rosenbrock(NVariable, LowerBound, UpperBound);       break;
        case STEP:             Property_Step(NVariable, LowerBound, UpperBound);             break;
        case QUARTIC_NOISE:    Property_QuarticNoise(NVariable, LowerBound, UpperBound);     break;
        case SCHWEFEL_S_2_26:  Property_Schwefel_s_2_26(NVariable, LowerBound, UpperBound);  break;
        case RASTRIGIN:        Property_Rastrigin(NVariable, LowerBound, UpperBound);        break;
        case ACKLEY:           Property_Ackley(NVariable, LowerBound, UpperBound);           break;
        case GRIEWANK:         Property_Griewank(NVariable, LowerBound, UpperBound);         break;
        default:               Property_Sphere(NVariable, LowerBound, UpperBound);           break;
    }
}

#endif // BENCHMARK_FUNCTION_H