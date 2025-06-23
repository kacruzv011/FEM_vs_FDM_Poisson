// --- ARCHIVO: src/poisson_sincronizacion.cpp (MODO DUAL) ---
// Este archivo puede funcionar en dos modos:
// 1. Benchmark (por defecto): Mide el tiempo y lo imprime en una línea.
// 2. Solution: Imprime la malla completa (x,y,V) para graficar.
// ADVERTENCIA: Esta implementación es conceptualmente diferente (usa matriz densa,
// solucionador iterativo) y es extremadamente ineficiente en memoria y tiempo.

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <omp.h>
#include <algorithm>

// --- ESTRUCTURAS Y TIPOS GLOBALES ---
struct Node { double x, y; bool is_boundary; };

// --- VARIABLES GLOBALES DE CONFIGURACIÓN ---
int case_selector = 1;
double x_ini, x_fin, y_ini, y_fin;
static constexpr double EPS = 1e-12;
static constexpr int MAX_ITERATIONS = 20000; // Reducido para evitar timeouts
static constexpr double TOLERANCE = 1e-6;

// --- FUNCIONES DEL PROBLEMA FÍSICO ---
double source_term(double x, double y) {
    switch (case_selector) {
        case 1: return (x * x + y * y) * std::exp(x * y);
        case 2: return 0.0;
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
    }
    return 0.0;
}

// --- FUNCIÓN PRINCIPAL DE CÁLCULO (HÍBRIDA/INEFICIENTE) ---
void solve_fem_sync(int M, int N, const std::vector<Node>& nodes, std::vector<double>& u) {
    int n_nodes = nodes.size();

    // 1. Ensamblaje usando matriz densa
    std::vector<std::vector<double>> K(n_nodes, std::vector<double>(n_nodes, 0.0));
    std::vector<double> F(n_nodes, 0.0);

    #pragma omp parallel for
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            int n0 = j * (M + 1) + i, n1 = n0 + 1, n2 = n0 + (M + 1), n3 = n2 + 1;
            int tri1_nodes[] = {n0, n1, n3};
            int tri2_nodes[] = {n0, n3, n2};
            
            for (const auto* tri_nodes : {tri1_nodes, tri2_nodes}) {
                double x[] = {nodes[tri_nodes[0]].x, nodes[tri_nodes[1]].x, nodes[tri_nodes[2]].x};
                double y[] = {nodes[tri_nodes[0]].y, nodes[tri_nodes[1]].y, nodes[tri_nodes[2]].y};
                double area = 0.5 * std::abs((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
                if (area < 1e-12) continue;

                double b[3] = {y[1]-y[2], y[2]-y[0], y[0]-y[1]};
                double c[3] = {x[2]-x[1], x[0]-x[2], x[1]-x[0]};
                double Ke[3][3];
                for(int r=0; r<3; ++r) for(int s=0; s<3; ++s) Ke[r][s] = (b[r]*b[s]+c[r]*c[s])/(4.0*area);

                double xc = (x[0]+x[1]+x[2])/3.0, yc = (y[0]+y[1]+y[2])/3.0;
                double f_val = source_term(xc, yc) * area / 3.0;
                
                #pragma omp critical
                {
                    for (int r = 0; r < 3; ++r) {
                        if (!nodes[tri_nodes[r]].is_boundary) {
                            F[tri_nodes[r]] += f_val;
                            for (int s = 0; s < 3; ++s) {
                                if (!nodes[tri_nodes[s]].is_boundary) {
                                    K[tri_nodes[r]][tri_nodes[s]] += Ke[r][s];
                                } else {
                                    F[tri_nodes[r]] -= Ke[r][s] * boundary_condition(nodes[tri_nodes[s]].x, nodes[tri_nodes[s]].y);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 2. Aplicar condiciones de frontera
    for (int i = 0; i < n_nodes; ++i) {
        if (nodes[i].is_boundary) {
            for (int j = 0; j < n_nodes; ++j) K[i][j] = 0.0;
            K[i][i] = 1.0;
            F[i] = boundary_condition(nodes[i].x, nodes[i].y);
        }
    }

    // 3. Solucionador iterativo (Jacobi)
    u.assign(n_nodes, 0.0);
    std::vector<double> u_old = u;
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        double max_diff = 0.0;
        #pragma omp parallel for reduction(max:max_diff)
        for (int i = 0; i < n_nodes; ++i) {
            if (nodes[i].is_boundary) {
                u[i] = F[i];
                continue;
            }
            double sigma = 0.0;
            for (int j = 0; j < n_nodes; ++j) {
                if (i != j) {
                    sigma += K[i][j] * u_old[j];
                }
            }
            double new_val = (F[i] - sigma) / K[i][i];
            double diff = std::abs(new_val - u[i]);
            if (diff > max_diff) max_diff = diff;
            u[i] = new_val;
        }
        if (max_diff < TOLERANCE) break;
        u_old = u;
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
    if (case_selector == 1) { x_ini = 0.0; x_fin = 2.0; y_ini = 0.0; y_fin = 1.0; }
    else if (case_selector == 2) { x_ini = 1.0; x_fin = 2.0; y_ini = 0.0; y_fin = 1.0; }
    else { std::cerr << "Caso no válido: " << case_selector << "\n"; return 1; }

    // 3. Generación de la malla
    int n_nodes = (M + 1) * (N + 1);
    std::vector<Node> nodes(n_nodes);
    double hx = (x_fin - x_ini) / M;
    double hy = (y_fin - y_ini) / N;
    for (int j = 0; j <= N; ++j) {
        for (int i = 0; i <= M; ++i) {
            int idx = j * (M + 1) + i;
            nodes[idx].x = x_ini + i * hx;
            nodes[idx].y = y_ini + j * hy;
            nodes[idx].is_boundary = (i == 0 || i == M || j == 0 || j == N);
        }
    }
    
    // 4. Ejecución y medición
    std::vector<double> u;
    double start_time = omp_get_wtime();
    solve_fem_sync(M, N, nodes, u);
    double end_time = omp_get_wtime();

    // 5. Salida según el modo
    if (mode == "solution") {
        std::cout << "x,y,V\n";
        for (int i = 0; i < n_nodes; ++i) {
            std::cout << nodes[i].x << "," << nodes[i].y << "," << u[i] << "\n";
        }
    } else { // modo "benchmark" por defecto
        std::cout << M << "," << N << "," << num_threads << "," << (end_time - start_time) << std::endl;
    }
    return 0;
}