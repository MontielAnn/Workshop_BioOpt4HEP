// =====================================================================================
// PSO BASE PARA EL MINICURSO ALGORITMOS BIOINSPIRADOS APLICADOS A LA FÍSICA DE PARTÍCULAS
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
//  CONFIGURACIÓN DE PARÁMETROS DEL PSO
// ====================================================================
#define N 20          // Número de partículas
#define KMAX 100      // Máximo número de iteraciones

#define ALPHA 1.0     // Aceleración
#define OMEGA 0.729   // Inercia
#define C1 1.494      // Coeficiente cognitivo 
#define C2 1.494      // Coeficiente social


// ====================================================================
// FUNCIONES AUXILIARES
// ====================================================================

double rand_double() {
    return (double)rand() / RAND_MAX;
}

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
    // --- INICIALIZACIÓN DE LA SEMILLA ALEATORIA ---
    unsigned int semilla = generarSemillaUnica();
    srand(semilla);
    printf(">>> Semilla generada: %u\n", semilla);

    // --- PROPIEDADES DE LA FUNCIÓN OBJETIVO ---
    int D;              // Dimensión
    double *lb;         // Límite inferior
    double *ub;         // Límite superior
    
    BenchmarkProperty(FUNCION_OBJETIVO, &D, &lb, &ub);

    printf(">>> Funcion objetivo: %d | Dimension: %d\n", FUNCION_OBJETIVO, D);
    printf(">>> Limites: [%g, %g]\n\n", lb[0], ub[0]); 

    double x[N][D];     // Posiciones de las partículas
    double v[N][D];     // Velocidades de las partículas
    double p[N][D];     // Mejores posiciones locales (pbest)
    double fp[N];       // Mejores fitness locales
    double g[D];        // Mejor posición global (gbest)
    double gfit;        // Mejor fitness global
    double fx[N];       // Fitness actual de cada partícula

    // --- INICIALIZACIÓN DE LAS PARTÍCULAS ---
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            x[i][j] = rand_double() * (ub[j] - lb[j]) + lb[j];
            v[i][j] = 0.0;
        }
    }

    // --- EVALUACIÓN INICIAL ---
    for (int i = 0; i < N; i++) {
        fx[i] = BenchmarkFunction(FUNCION_OBJETIVO, x[i], D);
        fp[i] = fx[i];
        for (int j = 0; j < D; j++) {
            p[i][j] = x[i][j];
        }
    }

    // --- ENCONTRAR LA MEJOR GLOBAL INICIAL ---
    gfit = fx[0];
    int ind = 0;
    for (int i = 1; i < N; i++) {
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

        // ------------------------------------------------------------
        // DEFINICIÓN DE PARÁMETROS
        // ------------------------------------------------------------
        double alpha = ALPHA;
        double omega = OMEGA;
        double c1    = C1;
        double c2    = C2;

        // ------------------------------------------------------------
        //  ACTUALIZACIÓN DE VELOCIDAD Y POSICIÓN
        // ------------------------------------------------------------
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < D; j++) {

                // Inercia 
                double term_inercia = omega * v[i][j];

                // Cognitivo 
                double r1 = rand_double();
                double term_cognitivo = c1 * r1 * (p[i][j] - x[i][j]);

                // Social 
                double r2 = rand_double();
                double term_social = c2 * r2 * (g[j] - x[i][j]);

                // Velocidad
                v[i][j] = alpha * (term_inercia + term_cognitivo + term_social);

                // Actualización de posición
                x[i][j] = x[i][j] + v[i][j];
            }
        }

        // ------------------------------------------------------------
        // VERIFICACIÓN DE LÍMITES
        // ------------------------------------------------------------
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < D; j++) {
                if (x[i][j] < lb[j]) {
                    x[i][j] = lb[j] + (0.25 * rand_double() * (ub[j] - lb[j]));
                }
                else if (x[i][j] > ub[j]) {
                    x[i][j] = ub[j] - (0.25 * rand_double() * (ub[j] - lb[j]));
                }
            }
        }

        // ------------------------------------------------------------
        // EVALUACIÓN DE LAS NUEVAS POSICIONES
        // ------------------------------------------------------------
        for (int i = 0; i < N; i++) {
            fx[i] = BenchmarkFunction(FUNCION_OBJETIVO, x[i], D);
        }

        // ------------------------------------------------------------
        // ACTUALIZACIÓN DE MEJORES LOCALES Y GLOBAL
        // ------------------------------------------------------------
        double mejor_fx_iter = fx[0];
        int mejor_ind_iter = 0;
        for (int i = 1; i < N; i++) {
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

        for (int i = 0; i < N; i++) {
            if (fx[i] < fp[i]) {
                fp[i] = fx[i];
                for (int j = 0; j < D; j++) {
                    p[i][j] = x[i][j];
                }
            }
        }

        // ------------------------------------------------------------
        // RESULTADOS
        // ------------------------------------------------------------
        printf("%3d | gbest = %12.8f | x = (%9.6f, %9.6f)\n",
               k, gfit, g[0], g[1]);
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

    // --- LIBERAR MEMORIA DE LOS LÍMITES ---
    free(lb);
    free(ub);

    return 0;
}