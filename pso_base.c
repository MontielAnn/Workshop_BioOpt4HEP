// =====================================================================================
// PSO BASE PARA EL MINICURSO ALGORITMOS BIOINSPIRADOS APLICADOS A LA FÍSICA DE PARTÍCULAS
// =====================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h> 

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
//  CONFIGURACIÓN DE PARÁMETROS
// ====================================================================

#define N 20          // Número de partículas
#define D 2           // Dimension del problema
#define KMAX 100      // Máximo número de iteraciones

#define ALPHA 1.0     // Aceleracion
#define OMEGA 0.729   // Inercia
#define C1 1.494      // Coeficiente cognitivo 
#define C2 1.494      // Coeficiente social

// Límites del espacio de búsqueda
#define LB -10.0
#define UB  10.0


// ====================================================================
// FUNCIÓN OBJETIVO (BENCHMARK)
// ====================================================================
// Sphere: Mínimo global en (0,0) con fitness 0. 
//    return x1 * x1 + x2 * x2;
//
// Rastrigin: Mínimo global en (0,0) con fitness 0. 
//    return 20 + (x1*x1 - 10*cos(2*M_PI*x1)) + (x2*x2 - 10*cos(2*M_PI*x2));
//
// Rosenbrock: Mínimo global en (1,1) con fitness 0. 
//    return 100*pow(x2 - x1*x1, 2) + pow(1 - x1, 2);
// ====================================================================

double funObj(double x1, double x2) {
    return x1 * x1 + x2 * x2;
}


// ====================================================================
// FUNCIONES AUXILIARES
// ====================================================================

// Genera un número aleatorio entre 0 y 1
double rand_double() {
    return (double)rand() / RAND_MAX;
}

// --------------------------------------------------------------------
// Genera una semilla única combinando tiempo actual en milisegundos y el PID del proceso
// --------------------------------------------------------------------
unsigned int generarSemillaUnica() {
    uint64_t tiempo_ms = 0;
    uint64_t pid = (uint64_t)GET_PID();

    #ifdef _WIN32
        // Windows: GetTickCount64 da milisegundos desde el arranque
        tiempo_ms = GetTickCount64();
    #else
        // Linux/macOS: clock_gettime
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        tiempo_ms = ((uint64_t)ts.tv_sec * 1000ULL) + ((uint64_t)ts.tv_nsec / 1000000ULL);
    #endif

    // Multiplicación en 64 bits 
    uint64_t semilla_larga = tiempo_ms * pid;

    // Mezcla de bits
    unsigned int semilla_final = (unsigned int)(semilla_larga ^ (semilla_larga >> 32));

    return semilla_final;
}


// ====================================================================
// MAIN
// ====================================================================

int main() {
    // INICIALIZACIÓN DE LA SEMILLA ALEATORIA 
    unsigned int semilla = generarSemillaUnica();
    srand(semilla);

    // VARIABLES DEL ENJAMBRE 
    double x[N][D];     // Posiciones de las partículas
    double v[N][D];     // Velocidades de las partículas
    double p[N][D];     // Mejores posiciones locales (pbest)
    double fp[N];       // Mejores fitness locales
    double g[D];        // Mejor posición global (gbest)
    double gfit;        // Mejor fitness global
    double fx[N];       // Fitness actual de cada partícula

    // --- INICIALIZACIÓN DE LAS PARTÍCULAS ---
    // Posiciones aleatorias dentro de [LB, UB] y velocidad inicial cero.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            x[i][j] = rand_double() * (UB - LB) + LB;
            v[i][j] = 0.0;
        }
    }

    // --- EVALUACIÓN INICIAL ---
    for (int i = 0; i < N; i++) {
        fx[i] = funObj(x[i][0], x[i][1]);
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
        // --------------------------------------------------------------------
        // EJMEPLOS
        //
        //   Inercia dinámica (lineal):
        //       double omega = 0.9 - (0.5 * (double)k / KMAX);
        //
        //   c1 con decaimiento exponencial:
        //       double c1 = 2.5 * exp(-2.0 * (double)k / KMAX);
        //
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
                
                //Verificación de límites
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < D; j++) {
                        if (x[i][j] < LB) {
                            x[i][j] = LB + (0.25 * ((long double)drand48()) * LB);
                        }
                        else if (x[i][j] > UB) {
                            x[i][j] = UB - (0.25 * ((long double)drand48()) * UB);
                        }
                    }
                }
            }
        }

        // ------------------------------------------------------------
        // EVALUACIÓN DE LAS NUEVAS POSICIONES
        // ------------------------------------------------------------
        for (int i = 0; i < N; i++) {
            fx[i] = funObj(x[i][0], x[i][1]);
        }

        // ------------------------------------------------------------
        // ACTUALIZACIÓN DE MEJORES LOCALES Y GLOBAL
        // ------------------------------------------------------------
        // Encontrar el mejor de esta iteración
        double mejor_fx_iter = fx[0];
        int mejor_ind_iter = 0;
        for (int i = 1; i < N; i++) {
            if (fx[i] < mejor_fx_iter) {
                mejor_fx_iter = fx[i];
                mejor_ind_iter = i;
            }
        }

        // Actualizar mejor global si es necesario
        if (mejor_fx_iter < gfit) {
            gfit = mejor_fx_iter;
            for (int j = 0; j < D; j++) {
                g[j] = x[mejor_ind_iter][j];
            }
        }

        // Actualizar mejores locales de cada partícula
        for (int i = 0; i < N; i++) {
            if (fx[i] < fp[i]) {
                fp[i] = fx[i];
                for (int j = 0; j < D; j++) {
                    p[i][j] = x[i][j];
                }
            }
        }

        // ------------------------------------------------------------
        //RESULTADOS
        // ------------------------------------------------------------
        printf("%3d | gbest = %12.8f | x = (%9.6f, %9.6f)\n",
               k, gfit, g[0], g[1]);
    }

    // ================================================================
    // RESULTADOS FINALES
    // ================================================================
    printf("\n========== FIN DE LA CORRIDA ==========\n");
    printf("Mejor solucion : [%.8f, %.8f]\n", g[0], g[1]);
    printf("Mejor fitness  : %.10f\n", gfit);
    printf("========================================\n");

    return 0;
}
