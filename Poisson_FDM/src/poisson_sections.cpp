// --- ARCHIVO: src/poisson_sections.cpp (FDM Refactorizado - MODO DUAL) ---
// Este archivo puede funcionar en dos modos:
// 1. Benchmark (por defecto): Mide el tiempo y lo imprime en una línea.
// 2. Solution: Imprime la malla completa (x,y,V) para graficar.
// Usa '#pragma omp sections' para paralelizar la inicialización y el cálculo de la fuente.

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <omp.h>

// --- VARIABLES GLOBALES DE CONFIGURACIÓN ---
int case_selector = 1;
double x_ini, x_fin, y_ini, y_fin;
static constexpr double EPS = 1e-12;
static constexpr int MAX_ITERATIONS = 500000;
static constexpr double TOLERANCE = 1e-7;

// --- FUNCIONES DEL PROBLEMA FÍSICO (UNIFICADAS) ---
double source_term(double x, double y) {
    switch (case_selector) {
        case 1: return (x*x + y*y) * std::exp(x * y);
        case 2: return 0.0;
        case 3: return 4.0;
        case 4: return x / y + y / x;
        default: return 0.0;
    }
}

double boundary_condition(double x, double y) {
    switch (case_selector) {
        case 1:
            if (std::abs(x - x_ini) < EPS) return 1.0;
            if (std::abs(x - x_fin) < EPS) return std::exp(2.0 * y);
            if (std::abs(y - y_ini) < EPS) return 1.0;
            if (std::abs(y - y_fin) < EPS) return std::exp(x);
            break;
        case 2:
            if (std::abs(x - x_ini) < EPS) return std::log(y*y + 1.0);
            if (std::abs(x - x_fin) < EPS) return std::log(y*y + 4.0);
            if (std::abs(y - y_ini) < EPS) return 2.0 * std::log(x);
            if (std::abs(y - y_fin) < EPS) return std::log(x*x + 1.0);
            break;
        case 3:
            if (std::abs(x - x_ini) < EPS) return (1.0 - y) * (1.0 - y);
            if (std::abs(x - x_fin) < EPS) return (2.0 - y) * (2.0 - y);
            if (std::abs(y - y_ini) < EPS) return x * x;
            if (std::abs(y - y_fin) < EPS) return (x - 2.0) * (x - 2.0);
            break;
        case 4:
            if (std::abs(x - x_ini) < EPS) return y * std::log(y);
            if (std::abs(x - x_fin) < EPS) return 2.0 * y * std::log(2.0 * y);
            if (std::abs(y - y_ini) < EPS) return x * std::log(x);
            if (std::abs(y - y_fin) < EPS) return 2.0 * x * std::log(2.0 * x);
            break;
    }
    return 0.0;
}

// --- FUNCIONES DE CÁLCULO ---
void initialize_grid(int M, int N, std::vector<std::vector<double>> &V, double h, double k) {
    V.assign(M + 1, std::vector<double>(N + 1, 0.0));
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (i == 0 || i == M || j == 0 || j == N) {
                V[i][j] = boundary_condition(x_ini + i * h, y_ini + j * k);
            }
        }
    }
}

void compute_source_matrix(int M, int N, std::vector<std::vector<double>> &F, double h, double k) {
    F.assign(M + 1, std::vector<double>(N + 1, 0.0));
    for (int i = 1; i < M; ++i) {
        for (int j = 1; j < N; ++j) {
            F[i][j] = source_term(x_ini + i * h, y_ini + j * k);
        }
    }
}

void solve_poisson_fdm_sections(int M, int N, std::vector<std::vector<double>> &V, double h, double k) {
    std::vector<std::vector<double>> F;

    // Paraleliza la inicialización y el cálculo de la fuente
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            initialize_grid(M, N, V, h, k);
        }
        #pragma omp section
        {
            compute_source_matrix(M, N, F, h, k);
        }
    }

    // El bucle de solución principal se ejecuta en paralelo después
    std::vector<std::vector<double>> V_old = V;
    double delta = TOLERANCE + 1.0;
    int iterations = 0;
    const double h2 = h * h;
    const double k2 = k * k;
    const double denom = 2.0 * (h2 + k2);

    while (delta > TOLERANCE && iterations < MAX_ITERATIONS) {
        delta = 0.0;
        V_old = V;
        
        #pragma omp parallel for collapse(2) reduction(max:delta)
        for (int i = 1; i < M; ++i) {
            for (int j = 1; j < N; ++j) {
                double numer = (V_old[i+1][j] + V_old[i-1][j]) * k2
                             + (V_old[i][j+1] + V_old[i][j-1]) * h2 - F[i][j] * h2 * k2;
                
                double V_new = numer / denom;
                double diff = std::abs(V_new - V[i][j]);
                
                V[i][j] = V_new;
                if (diff > delta) {
                    delta = diff;
                }
            }
        }
        iterations++;
    }
}

// --- FUNCIÓN MAIN CON LÓGICA DE MODO DUAL ---
int main(int argc, char *argv[]) {
    // 1. Parseo de argumentos
    if (argc < 5) {
        std::cerr << "Uso: " << argv[0] << " <M> <N> <hilos> <caso> [--mode solution]\n";
        return 1;
    }
    const int M = std::atoi(argv[1]);
    const int N = std::atoi(argv[2]);
    const int num_threads = std::atoi(argv[3]);
    case_selector = std::atoi(argv[4]);
    std::string mode = "benchmark";
    if (argc > 5 && std::string(argv[5]) == "--mode" && argc > 6) {
        mode = std::string(argv[6]);
    }

    // 2. Configuración del problema
    omp_set_num_threads(num_threads);
    if      (case_selector == 1) { x_ini=0.0; x_fin=2.0; y_ini=0.0; y_fin=1.0; }
    else if (case_selector == 2) { x_ini=1.0; x_fin=2.0; y_ini=0.0; y_fin=1.0; }
    else if (case_selector == 3) { x_ini=1.0; x_fin=2.0; y_ini=0.0; y_fin=2.0; }
    else if (case_selector == 4) { x_ini=1.0; x_fin=2.0; y_ini=1.0; y_fin=2.0; }
    else { std::cerr << "Caso no válido: " << case_selector << "\n"; return 1; }

    // 3. Preparación y Medición
    double h = (x_fin - x_ini) / M;
    double k = (y_fin - y_ini) / N;
    std::vector<std::vector<double>> V;

    double start_time = omp_get_wtime();
    solve_poisson_fdm_sections(M, N, V, h, k);
    double end_time = omp_get_wtime();

    // 4. Salida según el modo
    if (mode == "solution") {
        std::cout << "x,y,V\n";
        for (int j = 0; j <= N; ++j) {
            for (int i = 0; i <= M; ++i) {
                double x = x_ini + i * h;
                double y = y_ini + j * k;
                std::cout << x << "," << y << "," << V[i][j] << "\n";
            }
        }
    } else { // modo "benchmark" por defecto
        std::cout << M << "," << N << "," << num_threads << "," << (end_time - start_time) << std::endl;
    }
    return 0;
}