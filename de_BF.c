// =====================================================================================
// EVOLUCIÓN DIFERENCIAL PARA EL MINICURSO ALGORITMOS BIOINSPIRADOS
// =====================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h> 

#include "benchmark_function.h"

// Detección de Sistema Operativo ---
#ifdef _WIN32
    #include <process.h>
    #include <windows.h>
    #define GET_PID _getpid
#else
    #include <unistd.h>
    #define GET_PID getpid
#endif


// ====================================================================
//  SELECCIÓN DE LA FUNCIÓN OBJETIVO
//   Opciones disponibles:
//   SPHERE, SCHWEFEL_S_2_22, SCHWEFEL_S_1_20, ROSENBROCK, STEP,
//   QUARTIC_NOISE, SCHWEFEL_S_2_26, RASTRIGIN, ACKLEY, GRIEWANK
// ====================================================================
#define FUNCION_OBJETIVO  RASTRIGIN


// ====================================================================
//  CONFIGURACIÓN DE PARÁMETROS DEL DE
// ====================================================================
#define NP    50      // Tamaño de la población 
#define KMAX  500     // Máximo número de iteraciones (generaciones)

#define F     0.5     // Factor de mutación [0.4, 1.0]
#define CR    0.9     // Probabilidad de cruce [0.7, 1.0]


// Genera un número aleatorio entre 0 y 1
double rand_double() {
    return (double)rand() / RAND_MAX;
}

// Genera un entero aleatorio en [min, max] 
int rand_int(int min, int max) {
    return min + rand() % (max - min + 1);
}

// --------------------------------------------------------------------
// Genera una semilla única combinando tiempo actual en milisegundos y el PID del proceso
// --------------------------------------------------------------------
unsigned int generarSemillaUnica() {
    uint64_t tiempo_ms = 0;
    uint64_t pid = (uint64_t)GET_PID();

    #ifdef _WIN32
        tiempo_ms = GetTickCount64();
    #else
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        tiempo_ms = ((uint64_t)ts.tv_sec * 1000ULL) + ((uint64_t)ts.tv_nsec / 1000000ULL);
    #endif

    uint64_t semilla_larga = tiempo_ms * pid;
    unsigned int semilla_final = (unsigned int)(semilla_larga ^ (semilla_larga >> 32));

    return semilla_final;
}


// ====================================================================
// MAIN
// ====================================================================

int main() {

    unsigned int semilla = generarSemillaUnica();
    srand(semilla);
    printf(">>> Semilla generada: %u\n", semilla);

    // --- PROPIEDADES DE LA FUNCIÓN ---
    int D;              // Dimensión
    double *lb;         // Límite inferior
    double *ub;         // Límite superior
    
    BenchmarkProperty(FUNCION_OBJETIVO, &D, &lb, &ub);

    printf(">>> Funcion objetivo: %d | Dimension: %d\n", FUNCION_OBJETIVO, D);
    printf(">>> Limites: [%g, %g]\n\n", lb[0], ub[0]);

    // --- VARIABLES ---
    double x[NP][D];      // Población actual
    double trial[NP][D];  // Vectores prueba (mutantes + crossover)
    double fx[NP];        // Fitness de la población actual
    double ftrial[NP];    // Fitness de los vectores prueba
    double g[D];          // Mejor individuo global
    double gfit;          // Mejor fitness global

    // --- INICIALIZACIÓN DE LA POBLACIÓN ---
    for (int i = 0; i < NP; i++) {
        for (int j = 0; j < D; j++) {
            x[i][j] = rand_double() * (ub[j] - lb[j]) + lb[j];
        }
        fx[i] = BenchmarkFunction(FUNCION_OBJETIVO, x[i], D);
    }

    // --- ENCONTRAR EL MEJOR INDIVIDUO INICIAL ---
    gfit = fx[0];
    int ind = 0;
    for (int i = 1; i < NP; i++) {
        if (fx[i] < gfit) {
            gfit = fx[i];
            ind = i;
        }
    }
    for (int j = 0; j < D; j++) {
        g[j] = x[ind][j];
    }

    // ================================================================
    // PROCESO ITERATIVO
    // ================================================================
    for (int k = 1; k <= KMAX; k++) {

        for (int i = 0; i < NP; i++) {

            // --------------------------------------------------------
            // SELECCIÓN DE LOS 3 VECTORES ALEATORIOS (r1, r2, r3)
            // --------------------------------------------------------
            int r1, r2, r3;
            do { r1 = rand_int(0, NP - 1); } while (r1 == i);
            do { r2 = rand_int(0, NP - 1); } while (r2 == i || r2 == r1);
            do { r3 = rand_int(0, NP - 1); } while (r3 == i || r3 == r1 || r3 == r2);

            // --------------------------------------------------------
            // MUTACIÓN: 
            // --------------------------------------------------------
            double mutant[D];
            for (int j = 0; j < D; j++) {
                mutant[j] = x[r1][j] + F * (x[r2][j] - x[r3][j]);
            }

            // --------------------------------------------------------
            // CROSSOVER BINOMIAL
            // --------------------------------------------------------
            int jrand = rand_int(0, D - 1);  
            for (int j = 0; j < D; j++) {
                if (rand_double() < CR || j == jrand) {
                    trial[i][j] = mutant[j];
                } else {
                    trial[i][j] = x[i][j];
                }
            }

            // --------------------------------------------------------
            // VERIFICACIÓN DE LÍMITES
            // --------------------------------------------------------
            for (int j = 0; j < D; j++) {
                if (trial[i][j] < lb[j]) {
                    trial[i][j] = lb[j] + (0.25 * rand_double() * (ub[j] - lb[j]));
                }
                else if (trial[i][j] > ub[j]) {
                    trial[i][j] = ub[j] - (0.25 * rand_double() * (ub[j] - lb[j]));
                }
            }

            // --------------------------------------------------------
            // EVALUACIÓN DEL VECTOR PRUEBA
            // --------------------------------------------------------
            ftrial[i] = BenchmarkFunction(FUNCION_OBJETIVO, trial[i], D);

            // --------------------------------------------------------
            // SELECCIÓN 
            // --------------------------------------------------------
            if (ftrial[i] < fx[i]) {
                fx[i] = ftrial[i];
                for (int j = 0; j < D; j++) {
                    x[i][j] = trial[i][j];
                }
            }
        }

        // ------------------------------------------------------------
        // ACTUALIZACIÓN DEL MEJOR GLOBAL
        // ------------------------------------------------------------
        double mejor_fx_iter = fx[0];
        int mejor_ind_iter = 0;
        for (int i = 1; i < NP; i++) {
            if (fx[i] < mejor_fx_iter) {
                mejor_fx_iter = fx[i];
                mejor_ind_iter = i;
            }
        }

        if (mejor_fx_iter < gfit) {
            gfit = mejor_fx_iter;
            for (int j = 0; j < D; j++) {
                g[j] = x[mejor_ind_iter][j];
            }
        }

        // ------------------------------------------------------------
        // RESULTADOS
        // ------------------------------------------------------------
        printf("Gen %3d | gbest = %12.8f | x = (", k, gfit);
        for (int j = 0; j < D; j++) {
            printf("%9.6f%s", g[j], (j < D - 1) ? ", " : "");
        }
        printf(")\n");
    }

    // ================================================================
    // RESULTADOS FINALES
    // ================================================================
    printf("\n========== FIN DE LA CORRIDA ==========\n");
    printf("Mejor solucion : [");
    for (int j = 0; j < D; j++) {
        printf("%.8f%s", g[j], (j < D - 1) ? ", " : "");
    }
    printf("]\n");
    printf("Mejor fitness  : %.10f\n", gfit);
    printf("========================================\n");

    free(lb);
    free(ub);

    return 0;
}