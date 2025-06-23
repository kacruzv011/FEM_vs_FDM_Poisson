// --- ARCHIVO: src/poisson_schedule.cpp (MODO DUAL) ---
// Este archivo puede funcionar en dos modos:
// 1. Benchmark (por defecto): Mide el tiempo y lo imprime en una línea.
// 2. Solution: Imprime la malla completa (x,y,V) para graficar.
// Usa la cláusula 'schedule' de OpenMP para probar políticas static y dynamic.

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <omp.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// --- ESTRUCTURAS Y TIPOS GLOBALES ---
struct Node { double x, y; bool is_boundary; };
using SpMat = Eigen::SparseMatrix<double>;
using T = Eigen::Triplet<double>;

// --- VARIABLES GLOBALES DE CONFIGURACIÓN ---
int case_selector = 1;
double x_ini, x_fin, y_ini, y_fin;
static constexpr double EPS = 1e-12;

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

// --- FUNCIÓN PRINCIPAL DE CÁLCULO (PARALELA CON SCHEDULE) ---
void solve_fem_schedule(int M, int N, int analysis_type, const std::vector<Node>& nodes, Eigen::VectorXd& u) {
    int n_nodes = nodes.size();
    std::vector<std::pair<int, int>> quads;
    quads.reserve(M * N);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            quads.emplace_back(i, j);
        }
    }

    std::vector<T> global_triplets;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(n_nodes);

    #pragma omp parallel
    {
        std::vector<T> local_triplets;
        Eigen::VectorXd local_F = Eigen::VectorXd::Zero(n_nodes);

        if (analysis_type == 1) { // Estático
            #pragma omp for schedule(static)
            for (size_t k = 0; k < quads.size(); ++k) {
                int i = quads[k].first, j = quads[k].second;
                int n0 = j*(M+1)+i, n1=n0+1, n2=n0+(M+1), n3=n2+1;
                int tri1_nodes[] = {n0, n1, n3};
                int tri2_nodes[] = {n0, n3, n2};

                for (const auto* tri_nodes : {tri1_nodes, tri2_nodes}) {
                    double x[]={nodes[tri_nodes[0]].x, nodes[tri_nodes[1]].x, nodes[tri_nodes[2]].x};
                    double y[]={nodes[tri_nodes[0]].y, nodes[tri_nodes[1]].y, nodes[tri_nodes[2]].y};
                    double area=0.5*std::abs((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
                    if(area<1e-12) continue;
                    Eigen::Matrix<double,2,3> grad;
                    grad(0,0)=y[1]-y[2]; grad(0,1)=y[2]-y[0]; grad(0,2)=y[0]-y[1];
                    grad(1,0)=x[2]-x[1]; grad(1,1)=x[0]-x[2]; grad(1,2)=x[1]-x[0];
                    grad/=(2.0*area);
                    Eigen::Matrix3d Ke=area*(grad.transpose()*grad);
                    double xc=(x[0]+x[1]+x[2])/3.0, yc=(y[0]+y[1]+y[2])/3.0;
                    double f_val=source_term(xc,yc)*area/3.0;
                    for(int i_node=0;i_node<3;++i_node){
                        if(!nodes[tri_nodes[i_node]].is_boundary){
                            local_F(tri_nodes[i_node])+=f_val;
                            for(int j_node=0;j_node<3;++j_node){
                                if(!nodes[tri_nodes[j_node]].is_boundary){
                                    local_triplets.emplace_back(tri_nodes[i_node],tri_nodes[j_node],Ke(i_node,j_node));
                                }else{
                                    local_F(tri_nodes[i_node])-=Ke(i_node,j_node)*boundary_condition(nodes[tri_nodes[j_node]].x,nodes[tri_nodes[j_node]].y);
                                }
                            }
                        }
                    }
                }
            }
        } else { // Dinámico
            #pragma omp for schedule(dynamic)
            for (size_t k = 0; k < quads.size(); ++k) {
                int i = quads[k].first, j = quads[k].second;
                int n0 = j*(M+1)+i, n1=n0+1, n2=n0+(M+1), n3=n2+1;
                int tri1_nodes[] = {n0, n1, n3};
                int tri2_nodes[] = {n0, n3, n2};

                for (const auto* tri_nodes : {tri1_nodes, tri2_nodes}) {
                    double x[]={nodes[tri_nodes[0]].x, nodes[tri_nodes[1]].x, nodes[tri_nodes[2]].x};
                    double y[]={nodes[tri_nodes[0]].y, nodes[tri_nodes[1]].y, nodes[tri_nodes[2]].y};
                    double area=0.5*std::abs((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
                    if(area<1e-12) continue;
                    Eigen::Matrix<double,2,3> grad;
                    grad(0,0)=y[1]-y[2]; grad(0,1)=y[2]-y[0]; grad(0,2)=y[0]-y[1];
                    grad(1,0)=x[2]-x[1]; grad(1,1)=x[0]-x[2]; grad(1,2)=x[1]-x[0];
                    grad/=(2.0*area);
                    Eigen::Matrix3d Ke=area*(grad.transpose()*grad);
                    double xc=(x[0]+x[1]+x[2])/3.0, yc=(y[0]+y[1]+y[2])/3.0;
                    double f_val=source_term(xc,yc)*area/3.0;
                    for(int i_node=0;i_node<3;++i_node){
                        if(!nodes[tri_nodes[i_node]].is_boundary){
                            local_F(tri_nodes[i_node])+=f_val;
                            for(int j_node=0;j_node<3;++j_node){
                                if(!nodes[tri_nodes[j_node]].is_boundary){
                                    local_triplets.emplace_back(tri_nodes[i_node],tri_nodes[j_node],Ke(i_node,j_node));
                                }else{
                                    local_F(tri_nodes[i_node])-=Ke(i_node,j_node)*boundary_condition(nodes[tri_nodes[j_node]].x,nodes[tri_nodes[j_node]].y);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            global_triplets.insert(global_triplets.end(), local_triplets.begin(), local_triplets.end());
            F += local_F;
        }
    }

    SpMat K(n_nodes, n_nodes);
    K.setFromTriplets(global_triplets.begin(), global_triplets.end());
    for (int i = 0; i < n_nodes; ++i) {
        if (nodes[i].is_boundary) {
            K.coeffRef(i, i) = 1.0;
            F(i) = boundary_condition(nodes[i].x, nodes[i].y);
        }
    }
    
    Eigen::SparseLU<SpMat> solver(K);
    u = solver.solve(F);
}

// --- FUNCIÓN MAIN CON LÓGICA DE MODO DUAL ---
int main(int argc, char *argv[]) {
    // 1. Parseo de argumentos (este necesita uno más)
    if (argc < 6) {
        std::cerr << "Uso: " << argv[0] << " <M> <N> <hilos> <caso> <schedule_type (1=static, 2=dynamic)> [--mode solution]\n";
        return 1;
    }
    const int M = std::atoi(argv[1]);
    const int N = std::atoi(argv[2]);
    const int num_threads = std::atoi(argv[3]);
    case_selector = std::atoi(argv[4]);
    const int analysis_type = std::atoi(argv[5]);
    std::string mode = "benchmark";
    if (argc > 6 && std::string(argv[6]) == "--mode" && argc > 7) {
        mode = std::string(argv[7]);
    }

    // 2. Configuración del problema
    omp_set_num_threads(num_threads);
    if (case_selector == 1) { x_ini = 0.0; x_fin = 2.0; y_ini = 0.0; y_fin = 1.0; }
    else if (case_selector == 2) { x_ini = 1.0; x_fin = 2.0; y_ini = 0.0; y_fin = 1.0; }
    else { std::cerr << "Caso no válido\n"; return 1; }
    if (analysis_type != 1 && analysis_type != 2) { std::cerr << "Schedule no válido\n"; return 1; }

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
    Eigen::VectorXd u;
    double start_time = omp_get_wtime();
    solve_fem_schedule(M, N, analysis_type, nodes, u);
    double end_time = omp_get_wtime();

    // 5. Salida según el modo
    if (mode == "solution") {
        std::cout << "x,y,V\n";
        for (int i = 0; i < n_nodes; ++i) {
            std::cout << nodes[i].x << "," << nodes[i].y << "," << u(i) << "\n";
        }
    } else {
        // En modo benchmark, la salida incluye el tipo de schedule para el CSV
        std::cout << M << "," << N << "," << num_threads << "," << analysis_type << "," << (end_time - start_time) << std::endl;
    }
    return 0;
}