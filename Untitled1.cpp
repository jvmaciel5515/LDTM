#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm> // Para std::min e std::max
#include <map>       // Para o mapa de elementos para busca
#include <numeric>   // Para std::iota

// Inclusões da biblioteca Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeSolvers>

// Tipos Eigen para facilitar a leitura
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

// --- Definições de Tolerâncias e Constantes Globais ---
// Usando structs para unidades para segurança de tipo e clareza
struct Length {
    double value;
    explicit Length(double v) : value(v) {}
    Length operator+(const Length& other) const { return Length(value + other.value); }
    Length operator-(const Length& other) const { return Length(value - other.value); }
    Length operator*(double factor) const { return Length(value * factor); }
    Length operator/(double factor) const { return Length(value / factor); }
    Length operator/(const Length& other) const { return Length(value / other.value); } // Divisão entre Lengths
};

struct Force {
    double value;
    explicit Force(double v) : value(v) {}
    Force operator+(const Force& other) const { return Force(value + other.value); }
    Force operator-(const Force& other) const { return Force(value - other.value); }
    Force operator*(double factor) const { return Force(value * factor); }
    Force operator/(double factor) const { return Force(value / factor); }
};

struct Pressure { // Pascals
    double value;
    explicit Pressure(double v) : value(v) {}
    Pressure operator*(double factor) const { return Pressure(value * factor); }
    Force operator*(const Length& area) const { return Force(value * area.value); } // Simplified for 2D thickness
    Pressure operator/(double factor) const { return Pressure(value / factor); }
};

struct Strain { // Adimensional
    double value;
    explicit Strain(double v) : value(v) {}
    Strain operator*(double factor) const { return Strain(value * factor); }
    Strain operator/(double factor) const { return Strain(value / factor); }
};

struct Temperature { // Celsius
    double value;
    explicit Temperature(double v) : value(v) {}
};

// Operadores de conversão de double para as unidades
Length m(double v) { return Length(v); }
Force N(double v) { return Force(v); }
Pressure Pa(double v) { return Pressure(v); }
Strain dimensionless(double v) { return Strain(v); }
Temperature C(double v) { return Temperature(v); }


// Constantes
const double PI = 3.14159265358979323846;
const double SMALL_VALUE_DOUBLE = 1e-12; // Para evitar divisão por zero ou comparação com zero
const Pressure TOLERANCE_PRESSURE = Pa(1.0); // Tolerância para a força residual (Newton-Raphson) em Pascals -> Convertido para N depois
const Length TOLERANCE_DISPLACEMENT = m(1e-9); // Tolerância para o deslocamento incremental (Newton-Raphson) em metros
const Length CONTACT_GAP_TOLERANCE = m(1e-6); // Tolerância para detecção de contato (gap) em metros
const int MAX_NR_ITERATIONS = 50; // Número máximo de iterações do Newton-Raphson
const Temperature HOT_PROCESS_TEMPERATURE_THRESHOLD = C(500.0); // Limite para regime "a quente"

// Parâmetros para o passo adaptativo do Newton-Raphson
const double TARGET_NR_ITERATIONS = 7.0; // Número de iterações "ideal" para o NR
const double MIN_STEP_SIZE_FACTOR = 0.5; // Fator mínimo de ajuste do passo (para não reduzir demais)
const double MAX_STEP_SIZE_FACTOR = 1.5; // Fator máximo de ajuste do passo (para não aumentar demais)

const double decrease_factor = 0.5; // Fator de redução de alpha no line search

// Offset para a deformação plástica inicial na lei de Hollomon,
// para evitar log(0) ou potências de zero, e simular o "yield offset"
const double initial_plastic_strain_offset = 1e-6; // Pequeno valor

// Variável global para o número do passo atual (usado para exportação VTU)
int step_number = 0;
// Próximo ID de nó e elemento (para o gerenciamento de IDs únicos na remalhagem/refinamento)
int next_node_id = 0;
int next_element_id = 0;

// --- Estruturas de Dados ---
struct MaterialProperties {
    Pressure youngs_modulus;
    double poisson_ratio;
    Pressure yield_strength;
    Pressure strength_coefficient; // K na lei de hardening de Hollomon (Pa)
    double hardening_exponent;    // n na lei de hardening de Hollomon (adimensional)
};

// Dados por ponto de Gauss
struct GaussPointData {
    Vector sigma_cauchy; // Tensão de Cauchy [sigma_xx, sigma_yy, sigma_xy] (Pa)
    Matrix F_elastic;    // Gradiente de deformação elástico
    Matrix F_plastic;    // Gradiente de deformação plástico
    double equivalent_plastic_strain; // Deformação plástica equivalente acumulada
    Pressure current_yield_strength;  // Resistência ao escoamento atual (Pa)

    GaussPointData() : sigma_cauchy(3), F_elastic(2,2), F_plastic(2,2) {
        sigma_cauchy.setZero();
        F_elastic = Matrix::Identity(2,2);
        F_plastic = Matrix::Identity(2,2);
        equivalent_plastic_strain = 0.0;
        current_yield_strength = Pa(0.0); // Será inicializado com material_props
    }
    // Copiar construtor para facilitar a cópia em interpolação
    GaussPointData(const GaussPointData& other) = default;
};

struct Node {
    int id;
    Length x, y; // Coordenadas originais (referencial)
    Length ux, uy; // Deslocamentos acumulados
    Force fx_reaction, fy_reaction; // Forças de reação (para nós com BC)
    int contact_state; // 0: FREE, 1: STICK, 2: SLIP

    Node(int id, Length x_coord, Length y_coord) :
        id(id), x(x_coord), y(y_coord), ux(m(0.0)), uy(m(0.0)),
        fx_reaction(N(0.0)), fy_reaction(N(0.0)), contact_state(0) {}

    // Retorna a coordenada atual (deformada)
    Length current_x() const { return x + ux; }
    Length current_y() const { return y + uy; }

    void reset_displacements_and_forces() {
        ux = m(0.0);
        uy = m(0.0);
        fx_reaction = N(0.0);
        fy_reaction = N(0.0);
        contact_state = 0;
    }
};

struct Element {
    int id;
    std::vector<int> nodes; // IDs dos 4 nós do elemento (conectividade)
    std::vector<GaussPointData> gp_data; // Dados dos pontos de Gauss
    Pressure stress_xx, stress_yy, stress_xy; // Tensões médias no elemento
    double avg_plastic_strain; // Deformação plástica equivalente média no elemento

    Element(int id, int n1, int n2, int n3, int n4) :
        id(id), nodes({n1, n2, n3, n4}),
        stress_xx(Pa(0.0)), stress_yy(Pa(0.0)), stress_xy(Pa(0.0)),
        avg_plastic_strain(0.0) {
        // Inicializa os 4 pontos de Gauss
        gp_data.resize(4);
    }

    void reset_stresses() {
        stress_xx = Pa(0.0);
        stress_yy = Pa(0.0);
        stress_xy = Pa(0.0);
        avg_plastic_strain = 0.0;
    }

    void reset_gp_stresses() {
        for (auto& gp : gp_data) {
            gp.sigma_cauchy.setZero();
            gp.F_elastic = Matrix::Identity(2,2);
            gp.F_plastic = Matrix::Identity(2,2);
            gp.equivalent_plastic_strain = 0.0;
            // gp.current_yield_strength é atualizado pela função de plasticidade
        }
    }
};

struct StressStrainPoint {
    Strain strain;
    Pressure stress;
};

// --- Funções de Material ---
MaterialProperties get_material_properties_by_temperature(Temperature temp) {
    MaterialProperties props;
    if (temp.value >= HOT_PROCESS_TEMPERATURE_THRESHOLD.value) { // A quente
        props.youngs_modulus = Pa(150.0e9);  // Ex: Aço a quente - menor E
        props.poisson_ratio = 0.35;         // Ligeiramente maior nu
        props.yield_strength = Pa(80.0e6);   // Menor escoamento
        props.strength_coefficient = Pa(200.0e6); // Menor K
        props.hardening_exponent = 0.15;    // Menor n (menos endurecimento)
    } else { // A frio
        props.youngs_modulus = Pa(210.0e9); // Ex: Aço a frio - maior E
        props.poisson_ratio = 0.30;
        props.yield_strength = Pa(200.0e6);  // Maior escoamento
        props.strength_coefficient = Pa(400.0e6); // Maior K
        props.hardening_exponent = 0.25;    // Maior n (mais endurecimento)
    }
    return props;
}

// Funcao para calcular a 2ª Tensão de Piola-Kirchhoff de um tensor de Cauchy
// (para estado plano de tensão)
Matrix PK2_from_Cauchy(const Matrix& F_total, const Matrix& sigma_cauchy_matrix) {
    double J_total = F_total.determinant();
    Matrix F_total_inv = F_total.inverse();
    // S = J * F_inv * sigma_cauchy * F_inv_transpose
    return J_total * F_total_inv * sigma_cauchy_matrix * F_total_inv.transpose();
}


// Implementação do Algoritmo de Retorno ao Plano (Radial Return Mapping)
// para Plasticidade J2 (Von Mises) com Endurecimento Isotrópico (Lei de Hollomon)
// Assume Estado Plano de Tensão (Plane Stress)
void plasticity_stress_update(const Matrix& F_total, GaussPointData& gp_data,
                              const MaterialProperties& material_props, Matrix& C_ep_consistent) {

    double E_mod = material_props.youngs_modulus.value;
    double nu = material_props.poisson_ratio;
    double K = material_props.strength_coefficient.value;
    double n = material_props.hardening_exponent;
    double sigma_y0 = material_props.yield_strength.value;

    // Módulos elásticos (para estado plano de tensão)
    double G = E_mod / (2.0 * (1.0 + nu)); // Módulo de cisalhamento

    // Matriz constitutiva elástica para estado plano de tensão (Voigt notation [xx, yy, xy])
    // C_elastic_plane_stress mapeia d_strain_vec para d_stress_vec
    Matrix C_elastic_plane_stress(3,3);
    double factor = E_mod / (1.0 - nu * nu);
    C_elastic_plane_stress << factor, factor * nu, 0.0,
                              factor * nu, factor, 0.0,
                              0.0, 0.0, G; // C_ep_consistent(2,2) para cisalhamento xy

    // Gradiente de deformação elástico trial F_e_trial = F_total * F_p_inv
    Matrix F_plastic_inverse = gp_data.F_plastic.inverse();
    Matrix F_elastic_trial = F_total * F_plastic_inverse;

    // Calcular a deformação de Hencky (logarítmica) a partir de F_elastic_trial
    // Para elasticidade isotrópica, Hencky_strain = 0.5 * log(F_e_trial^T * F_e_trial)
    Matrix Be_trial = F_elastic_trial * F_elastic_trial.transpose();
    Eigen::SelfAdjointEigenSolver<Matrix> es(Be_trial);
    Vector lambda_e_squared = es.eigenvalues();
    Matrix R_spectral = es.eigenvectors(); // R são os auto-vetores

    Matrix log_sqrt_Be_trial(2,2);
    log_sqrt_Be_trial << 0.5 * std::log(lambda_e_squared(0)), 0.0,
                         0.0, 0.5 * std::log(lambda_e_squared(1));

    Matrix Hencky_strain_trial = R_spectral * log_sqrt_Be_trial * R_spectral.transpose();

    // Tensão de Cauchy "Trial" (em notação de Voigt: [sigma_xx, sigma_yy, sigma_xy])
    Vector epsilon_hencky_vec(3);
    epsilon_hencky_vec(0) = Hencky_strain_trial(0,0);
    epsilon_hencky_vec(1) = Hencky_strain_trial(1,1);
    epsilon_hencky_vec(2) = 2.0 * Hencky_strain_trial(0,1); // Shear strain gamma_xy = 2 * epsilon_xy

    Vector sigma_trial_vec = C_elastic_plane_stress * epsilon_hencky_vec;

    // Calcular a tensão de Von Mises "Trial"
    double sigma_x = sigma_trial_vec(0);
    double sigma_y = sigma_trial_vec(1);
    double tau_xy = sigma_trial_vec(2); // sigma_xy

    // Tensão desviatória "Trial"
    double sigma_m = (sigma_x + sigma_y) / 2.0; // Tensão hidrostática
    double s_x = sigma_x - sigma_m;
    double s_y = sigma_y - sigma_m;
    double s_xy = tau_xy;

    double sigma_vm_trial = std::sqrt(s_x*s_x + s_y*s_y + 2.0*s_xy*s_xy); // Von Mises para estado plano de tensão

    // Atualiza a tensão de escoamento com o endurecimento atual
    double current_yield = K * std::pow(initial_plastic_strain_offset + gp_data.equivalent_plastic_strain, n);
    current_yield = std::max(sigma_y0, current_yield);
    gp_data.current_yield_strength = Pa(current_yield);

    // --- Retorno ao Plano ---
    if (sigma_vm_trial <= current_yield + SMALL_VALUE_DOUBLE) { // Estado Elástico
        gp_data.sigma_cauchy = sigma_trial_vec;
        gp_data.F_elastic = F_elastic_trial;
        // F_plastic e equivalent_plastic_strain não mudam
        C_ep_consistent = C_elastic_plane_stress; // Matriz tangente elástica
    } else { // Estado Plástico - Retorno ao Plano
        double delta_gamma = 0.0; // Incremento do multiplicador plástico

        // Tensão desviatória trial (para calculo do r_trial)
        Vector s_trial_vec(3); // [s_x, s_y, s_xy]
        s_trial_vec << s_x, s_y, s_xy;

        double r_trial = std::sqrt(s_x*s_x + s_y*s_y + 2.0*s_xy*s_xy); // sqrt(2*J2)

        // Iteração de Newton-Raphson para delta_gamma
        double dg_prev = 0.0; // Chute inicial para delta_gamma
        double tolerance_dg = 1e-9;
        int max_dg_iter = 100;

        for (int iter_dg = 0; iter_dg < max_dg_iter; ++iter_dg) {
            double current_total_plastic_strain_trial = gp_data.equivalent_plastic_strain + dg_prev;
            if (current_total_plastic_strain_trial < 0) current_total_plastic_strain_trial = 0;

            double sigma_y_plastic = K * std::pow(initial_plastic_strain_offset + current_total_plastic_strain_trial, n);
            if (sigma_y_plastic < sigma_y0) sigma_y_plastic = sigma_y0;

            double f_dg = r_trial - 2.0 * G * dg_prev - sigma_y_plastic;

            double df_dg_ddg = -2.0 * G - n * K * std::pow(initial_plastic_strain_offset + current_total_plastic_strain_trial, n-1);
            if (std::abs(df_dg_ddg) < SMALL_VALUE_DOUBLE) {
                df_dg_ddg = -SMALL_VALUE_DOUBLE;
            }

            double delta_dg = -f_dg / df_dg_ddg;
            dg_prev += delta_dg;

            if (std::abs(delta_dg) < tolerance_dg) {
                break;
            }
        }
        delta_gamma = dg_prev;
        if (delta_gamma < 0.0) delta_gamma = 0.0;

        // Atualiza a deformação plástica equivalente acumulada
        gp_data.equivalent_plastic_strain += delta_gamma;

        // Atualiza a tensão de escoamento
        gp_data.current_yield_strength = Pa(K * std::pow(initial_plastic_strain_offset + gp_data.equivalent_plastic_strain, n));
        if (gp_data.current_yield_strength.value < sigma_y0) gp_data.current_yield_strength = Pa(sigma_y0);


        // Tensão de Cauchy plastificada (retornada ao plano de escoamento)
        double eta = (gp_data.current_yield_strength.value) / (r_trial + SMALL_VALUE_DOUBLE);
        Vector sigma_cauchy_plastic_vec(3);
        sigma_cauchy_plastic_vec(0) = eta * s_x + sigma_m;
        sigma_cauchy_plastic_vec(1) = eta * s_y + sigma_m;
        sigma_cauchy_plastic_vec(2) = eta * s_xy; // sigma_xy

        gp_data.sigma_cauchy = sigma_cauchy_plastic_vec;

        // Atualiza F_plastic e F_elastic para grandes deformações
        // Para uma atualização consistente, usa-se a exponencial de matriz.
        // d_epsilon_p = delta_gamma * (s_vec_norm / r_trial)
        // delta_epsilon_p_mat = delta_gamma * N_dev, onde N_dev = s_dev / norm(s_dev)

        Matrix s_dev_mat(2,2);
        s_dev_mat << s_x, s_xy,
                     s_xy, s_y;

        Matrix N_dev_mat = s_dev_mat / (r_trial + SMALL_VALUE_DOUBLE);

        // F_p_new = exp(delta_gamma * N_dev) * F_p_old
        // Para uma exponencial de matriz 2x2, Eigen tem um método para matrizes simétricas,
        // ou podemos usar a série de Taylor (para pequenos incrementos)
        // exp(X) approx I + X + X^2/2! + ...
        // Para esta implementação, vamos usar uma aproximação de ordem 1: exp(X) ~ I + X
        // Isso é uma simplificação comum em muitos códigos FEM quando o passo de tempo é pequeno.
        Matrix exp_plastic_increment = Matrix::Identity(2,2) + (delta_gamma * N_dev_mat);

        gp_data.F_plastic = exp_plastic_increment * gp_data.F_plastic;
        gp_data.F_elastic = F_total * gp_data.F_plastic.inverse();

        // --- Matriz Tangente Consistente (C_ep_consistent) ---
        double H_prime = n * K * std::pow(initial_plastic_strain_offset + gp_data.equivalent_plastic_strain, n-1);
        if (gp_data.equivalent_plastic_strain < SMALL_VALUE_DOUBLE) H_prime = 0.0;

        double M_factor_denominator = (2.0 * G + H_prime);
        if (std::abs(M_factor_denominator) < SMALL_VALUE_DOUBLE) M_factor_denominator = SMALL_VALUE_DOUBLE; // Evita div por zero
        double M = 2.0 * G * G / M_factor_denominator;

        // Tensão desviatória final (após retorno ao plano)
        double final_s_x = gp_data.sigma_cauchy(0) - (gp_data.sigma_cauchy(0) + gp_data.sigma_cauchy(1)) / 2.0;
        double final_s_y = gp_data.sigma_cauchy(1) - (gp_data.sigma_cauchy(0) + gp_data.sigma_cauchy(1)) / 2.0;
        double final_s_xy = gp_data.sigma_cauchy(2);

        double final_r = std::sqrt(final_s_x*final_s_x + final_s_y*final_s_y + 2.0*final_s_xy*final_s_xy);

        // n_dev (vetor unitário na direção desviatória da tensão)
        Vector n_dev_vec(3); // [nx, ny, nxy]
        n_dev_vec(0) = final_s_x / (final_r + SMALL_VALUE_DOUBLE);
        n_dev_vec(1) = final_s_y / (final_r + SMALL_VALUE_DOUBLE);
        n_dev_vec(2) = final_s_xy / (final_r + SMALL_VALUE_DOUBLE);

        Matrix n_dev_dyadic = n_dev_vec * n_dev_vec.transpose();

        C_ep_consistent = C_elastic_plane_stress - M * n_dev_dyadic;
    }
}

// --- Funções de MEF ---

// Obtém as coordenadas atuais dos nós de um elemento
std::vector<Vector> get_nodal_current_coordinates_for_element(const Element& element, const std::vector<Node>& nodes) {
    std::vector<Vector> current_coords(4);
    for (int i = 0; i < 4; ++i) {
        int global_node_id = element.nodes[i];
        current_coords[i] = Vector(2);
        current_coords[i] << nodes[global_node_id].current_x().value, nodes[global_node_id].current_y().value;
    }
    return current_coords;
}

// Calcula as funções de forma (shape functions) N e seus gradientes
// em relação a xi e eta para um elemento quadrilátero isoparamétrico.
void calculate_shape_functions_and_derivatives(double xi, double eta,
                                               double N[4], double dNdxi[4], double dNdeta[4]) {
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);

    dNdxi[0] = -0.25 * (1.0 - eta);
    dNdxi[1] = 0.25 * (1.0 - eta);
    dNdxi[2] = 0.25 * (1.0 + eta);
    dNdxi[3] = -0.25 * (1.0 + eta);

    dNdeta[0] = -0.25 * (1.0 - xi);
    dNdeta[1] = -0.25 * (1.0 + xi);
    dNdeta[2] = 0.25 * (1.0 + xi);
    dNdeta[3] = 0.25 * (1.0 + xi); // Corrigido aqui: dNdeta[3] = 0.25 * (1.0 - xi)
}

// Calcula as funções de forma N e seus gradientes em relação a X e Y (no espaço atual, deformado)
// para um dado ponto (x,y) e os nós do elemento (current_node_coords)
// Retorna (true, N_vec, B_matrix, detJ_current) se o ponto está dentro do elemento
// Retorna (false, ...) caso contrário (e detJ_current negativo)
bool calculate_N_and_B_at_point(const std::vector<Vector>& current_node_coords,
                                double x_target, double y_target,
                                Vector& N_vec_out, Matrix& B_matrix_out, double& detJ_current_out,
                                Vector& natural_coords_out) {

    // Iteração para encontrar as coordenadas naturais (xi, eta) do ponto (x_target, y_target)
    // Usando Newton-Raphson 2D
    natural_coords_out = Vector(2);
    natural_coords_out << 0.0, 0.0; // Chute inicial (centro do elemento natural)

    double tolerance = 1e-6;
    int max_iter = 50;

    double N_val[4], dNdxi_val[4], dNdeta_val[4];

    for (int iter = 0; iter < max_iter; ++iter) {
        double xi = natural_coords_out(0);
        double eta = natural_coords_out(1);

        calculate_shape_functions_and_derivatives(xi, eta, N_val, dNdxi_val, dNdeta_val);

        double x_current = 0.0, y_current = 0.0;
        for (int i = 0; i < 4; ++i) {
            x_current += N_val[i] * current_node_coords[i](0);
            y_current += N_val[i] * current_node_coords[i](1);
        }

        Vector R_vec(2); // Residual
        R_vec << x_current - x_target, y_current - y_target;

        if (R_vec.norm() < tolerance) {
            // Ponto encontrado, verificar se está dentro do elemento natural [-1, 1]
            if (std::abs(xi) <= 1.0 + SMALL_VALUE_DOUBLE && std::abs(eta) <= 1.0 + SMALL_VALUE_DOUBLE) {
                N_vec_out = Vector(4);
                N_vec_out << N_val[0], N_val[1], N_val[2], N_val[3];

                // Jacobian da transformação (de xi,eta para x,y no espaço deformado)
                Matrix J_current(2,2);
                J_current(0,0) = dNdxi_val[0]*current_node_coords[0](0) + dNdxi_val[1]*current_node_coords[1](0) + dNdxi_val[2]*current_node_coords[2](0) + dNdxi_val[3]*current_node_coords[3](0);
                J_current(0,1) = dNdxi_val[0]*current_node_coords[0](1) + dNdxi_val[1]*current_node_coords[1](1) + dNdxi_val[2]*current_node_coords[2](1) + dNdxi_val[3]*current_node_coords[3](1);
                J_current(1,0) = dNdeta_val[0]*current_node_coords[0](0) + dNdeta_val[1]*current_node_coords[1](0) + dNdeta_val[2]*current_node_coords[2](0) + dNdeta_val[3]*current_node_coords[3](0);
                J_current(1,1) = dNdeta_val[0]*current_node_coords[0](1) + dNdeta_val[1]*current_node_coords[1](1) + dNdeta_val[2]*current_node_coords[2](1) + dNdeta_val[3]*current_node_coords[3](1);

                detJ_current_out = J_current.determinant();
                if (detJ_current_out < SMALL_VALUE_DOUBLE) { // Elemento degenerado/invertido
                    return false;
                }
                Matrix J_current_inv = J_current.inverse();

                // Gradientes das funções de forma em relação a X e Y (no espaço deformado)
                B_matrix_out.setZero(3, 8); // Matriz B generalizada (para deformações)
                for (int k = 0; k < 4; ++k) {
                    double dNdX = J_current_inv(0,0) * dNdxi_val[k] + J_current_inv(0,1) * dNdeta_val[k];
                    double dNdY = J_current_inv(1,0) * dNdxi_val[k] + J_current_inv(1,1) * dNdeta_val[k];
                    B_matrix_out(0, 2*k) = dNdX;
                    B_matrix_out(1, 2*k+1) = dNdY;
                    B_matrix_out(2, 2*k) = dNdY;
                    B_matrix_out(2, 2*k+1) = dNdX;
                }
                return true;
            } else {
                return false; // Ponto fora do domínio natural [-1,1]
            }
        }

        // Jacobian da transformação (de xi,eta para x,y no espaço deformado)
        Matrix J_current_newton(2,2);
        J_current_newton(0,0) = dNdxi_val[0]*current_node_coords[0](0) + dNdxi_val[1]*current_node_coords[1](0) + dNdxi_val[2]*current_node_coords[2](0) + dNdxi_val[3]*current_node_coords[3](0);
        J_current_newton(0,1) = dNdxi_val[0]*current_node_coords[0](1) + dNdxi_val[1]*current_node_coords[1](1) + dNdxi_val[2]*current_node_coords[2](1) + dNdxi_val[3]*current_node_coords[3](1);
        J_current_newton(1,0) = dNdeta_val[0]*current_node_coords[0](0) + dNdeta_val[1]*current_node_coords[1](0) + dNdeta_val[2]*current_node_coords[2](0) + dNdeta_val[3]*current_node_coords[3](0);
        J_current_newton(1,1) = dNdeta_val[0]*current_node_coords[0](1) + dNdeta_val[1]*current_node_coords[1](1) + dNdeta_val[2]*current_node_coords[2](1) + dNdeta_val[3]*current_node_coords[3](1);

        detJ_current_out = J_current_newton.determinant();
        if (detJ_current_out < SMALL_VALUE_DOUBLE) { // Elemento degenerado
             return false;
        }

        Vector delta_natural_coords = J_current_newton.inverse() * R_vec;
        natural_coords_out -= delta_natural_coords;
    }
    return false; // Não convergiu
}


// Encontra o elemento que contém um determinado ponto (x_target, y_target) no espaço deformado
// Retorna um ponteiro para o elemento ou nullptr se não encontrado.
Element* find_element_containing_point(const std::vector<Element>& elements, const std::vector<Node>& nodes,
                                       double x_target, double y_target) {
    for (auto& element : elements) {
        std::vector<Vector> current_node_coords = get_nodal_current_coordinates_for_element(element, nodes);

        Vector N_vec_dummy;
        Matrix B_matrix_dummy;
        double detJ_dummy;
        Vector natural_coords_dummy;

        if (calculate_N_and_B_at_point(current_node_coords, x_target, y_target,
                                       N_vec_dummy, B_matrix_dummy, detJ_dummy, natural_coords_dummy)) {
            return &element;
        }
    }
    return nullptr;
}


// Interpola os dados dos pontos de Gauss de um elemento para um ponto arbitrário
// Isso geralmente é feito usando as funções de forma do elemento.
GaussPointData interpolate_gp_data_to_point(const Element& element,
                                            const std::vector<Vector>& current_node_coords,
                                            double x_target, double y_target) {

    GaussPointData interpolated_data;
    Vector N_vec;
    Matrix B_matrix;
    double detJ;
    Vector natural_coords;

    if (!calculate_N_and_B_at_point(current_node_coords, x_target, y_target, N_vec, B_matrix, detJ, natural_coords)) {
        // Isso não deveria acontecer se find_element_containing_point funcionou corretamente.
        // Retorna dados zerados ou padrão.
        std::cerr << "AVISO: Ponto fora do elemento durante interpolacao. Retornando dados zerados.\n";
        return interpolated_data;
    }

    // Interpolação direta dos valores de GP para o ponto
    // Para 4 GPs e 4 Ns, cada N está associado a um nó, não a um GP.
    // A interpolação de GP para ponto é geralmente feita mapeando os GPs para o espaço natural
    // e usando as funções de forma nos pontos de Gauss.
    // Para simplificar, faremos uma interpolação bilineiar baseada nas coordenadas naturais (xi, eta)
    // para cada GP.

    // Pontos de Gauss e pesos para interpolação
    const double gp_coord = 1.0 / std::sqrt(3.0);
    const double gauss_points_coords[4][2] = {
        {-gp_coord, -gp_coord},
        { gp_coord, -gp_coord},
        { gp_coord,  gp_coord},
        {-gp_coord,  gp_coord}
    };

    // N(xi, eta) para interpolar os valores dos GPs.
    // Usamos as coordenadas naturais do ponto de destino (natural_coords(0), natural_coords(1))
    // para interpolar *dos GPs*. Não é interpolação nodal aqui.
    // Isso é mais complexo do que simplesmente usar N_vec.
    // Uma forma simples é uma média ponderada pelos inverso da distância ou simplesmente a média.
    // Para interpolação mais rigorosa, precisamos das funções de forma em cada GP
    // e depois interpolar usando as funções de forma no ponto de interesse.

    // A interpolação de dados de GP para um ponto arbitrário é mais complexa do que N_vec.
    // Uma abordagem comum é primeiro extrapolar os dados dos GPs para os nós do elemento
    // e depois interpolar dos nós para o ponto. Isso "suaviza" os resultados.
    // Para simplificar aqui, faremos uma interpolação direta dos valores dos GPs para o ponto (xi, eta)
    // usando as funções de forma. Isso pode não ser numericamente consistente para tensores.

    // Neste caso, se temos a posição (xi, eta) do ponto no elemento natural,
    // podemos usar as funções de forma para interpolar DADOS DE NÓS para esse ponto.
    // Mas aqui temos DADOS DE GPs.
    // Uma aproximação para a interpolação de GP é usar os N_val no ponto (xi,eta)
    // para ponderar os GPs mais próximos ou fazer uma média simples.

    // A interpolação correta de valores de GP para um novo GP exige as funções de forma avaliadas nas coordenadas
    // naturais de cada GP no elemento de ORIGEM. E depois aplicar essas funções de forma aos valores dos GPs.

    // Para simplificar a demonstração e evitar uma biblioteca de interpolação complexa:
    // Farei uma média simples dos dados dos GPs do elemento de origem.
    // Em um sistema real, essa interpolação precisa ser mais rigorosa, especialmente para tensores como F_elastic e F_plastic.

    interpolated_data.sigma_cauchy.setZero();
    interpolated_data.F_elastic.setZero();
    interpolated_data.F_plastic.setZero();
    interpolated_data.equivalent_plastic_strain = 0.0;
    interpolated_data.current_yield_strength = Pa(0.0);

    for(const auto& gp : element.gp_data) {
        interpolated_data.sigma_cauchy += gp.sigma_cauchy;
        interpolated_data.F_elastic += gp.F_elastic;
        interpolated_data.F_plastic += gp.F_plastic;
        interpolated_data.equivalent_plastic_strain += gp.equivalent_plastic_strain;
        interpolated_data.current_yield_strength = Pa(interpolated_data.current_yield_strength.value + gp.current_yield_strength.value);
    }
    interpolated_data.sigma_cauchy /= element.gp_data.size();
    interpolated_data.F_elastic /= element.gp_data.size();
    interpolated_data.F_plastic /= element.gp_data.size();
    interpolated_data.equivalent_plastic_strain /= element.gp_data.size();
    interpolated_data.current_yield_strength = Pa(interpolated_data.current_yield_strength.value / element.gp_data.size());

    // Se o ponto não está dentro do elemento, os dados retornados não serão válidos.
    // A função `find_element_containing_point` garante que o ponto está no elemento.

    return interpolated_data;
}


// Discretiza um retângulo em elementos quadriláteros de 4 nós
void discretize_rectangle(Length width, Length height, int divisions_x, int divisions_y,
                          std::vector<Node>& nodes, std::vector<Element>& elements) {
    nodes.clear();
    elements.clear();
    next_node_id = 0; // Resetar IDs
    next_element_id = 0;

    // Cria os nós
    for (int j = 0; j <= divisions_y; ++j) {
        for (int i = 0; i <= divisions_x; ++i) {
            Length x_coord = m(static_cast<double>(i) / divisions_x * width.value);
            Length y_coord = m(static_cast<double>(j) / divisions_y * height.value);
            nodes.emplace_back(next_node_id++, x_coord, y_coord);
        }
    }

    // Cria os elementos
    for (int j = 0; j < divisions_y; ++j) {
        for (int i = 0; i < divisions_x; ++i) {
            int n1 = j * (divisions_x + 1) + i;
            int n2 = n1 + 1;
            int n3 = n1 + (divisions_x + 1) + 1;
            int n4 = n1 + (divisions_x + 1);
            elements.emplace_back(next_element_id++, n1, n2, n3, n4);
        }
    }
    std::cout << "Malha criada: " << nodes.size() << " nos, " << elements.size() << " elementos.\n";
}

// Calcula a matriz de rigidez tangente do elemento e o vetor de forças internas
// para grandes deformações (material + geométrico)
// Retorna a matriz de rigidez tangente do elemento
Matrix calculate_element_tangent_stiffness_and_internal_forces(Element& element,
                                                               const std::vector<Node>& nodes,
                                                               const MaterialProperties& material_props,
                                                               Vector& F_int_element, // Saída: Forças internas do elemento
                                                               Matrix& K_element_material_out, // Saída: Matriz de rigidez material
                                                               Matrix& K_element_geometric_out) { // Saída: Matriz de rigidez geométrica

    K_element_material_out.setZero(8, 8); // Matriz de rigidez material do elemento (Pa * m)
    K_element_geometric_out.setZero(8, 8); // Matriz de rigidez geométrica do elemento (Pa * m)
    F_int_element.setZero(8);             // Vetor de forças internas do elemento (N)

    // Coordenadas originais dos nós do elemento
    Vector x_coords(4), y_coords(4);
    // Coordenadas dos deslocamentos dos nós do elemento (para calcular o Grad_u)
    Vector u_coords(4), v_coords(4);
    // Coordenadas atuais (deformadas) dos nós do elemento
    Vector current_x_coords(4), current_y_coords(4);

    for (int i = 0; i < 4; ++i) {
        int global_node_id = element.nodes[i];
        x_coords(i) = nodes[global_node_id].x.value;
        y_coords(i) = nodes[global_node_id].y.value;
        u_coords(i) = nodes[global_node_id].ux.value; // Deslocamento acumulado em X
        v_coords(i) = nodes[global_node_id].uy.value; // Deslocamento acumulado em Y
        current_x_coords(i) = nodes[global_node_id].current_x().value;
        current_y_coords(i) = nodes[global_node_id].current_y().value;
    }

    // Pontos de Gauss e pesos para integração 2x2
    const double gp_coord = 1.0 / std::sqrt(3.0);
    const double gauss_points[4][2] = {
        {-gp_coord, -gp_coord},
        { gp_coord, -gp_coord},
        { gp_coord,  gp_coord},
        {-gp_coord,  gp_coord}
    };
    const double gp_weights[4] = {1.0, 1.0, 1.0, 1.0}; // Pesos para integração 2x2

    int gp_idx = 0;
    for (int i = 0; i < 2; ++i) { // Loop sobre xi
        for (int j = 0; j < 2; ++j) { // Loop sobre eta
            double xi = gauss_points[gp_idx][0];
            double eta = gauss_points[gp_idx][1];
            double w_xi = gp_weights[gp_idx]; // Peso em xi
            double w_eta = gp_weights[gp_idx]; // Peso em eta

            double N[4];
            double dNdxi[4];
            double dNdeta[4];
            calculate_shape_functions_and_derivatives(xi, eta, N, dNdxi, dNdeta);

            // Matriz Jacobiana J0 (para mapeamento do espaço natural para o espaço material - referencial)
            // J0 = [dx/dxi, dy/dxi; dx/deta, dy/deta]
            Matrix J0(2,2);
            J0(0,0) = dNdxi[0]*x_coords(0) + dNdxi[1]*x_coords(1) + dNdxi[2]*x_coords(2) + dNdxi[3]*x_coords(3);
            J0(0,1) = dNdxi[0]*y_coords(0) + dNdxi[1]*y_coords(1) + dNdxi[2]*y_coords(2) + dNdxi[3]*y_coords(3);
            J0(1,0) = dNdeta[0]*x_coords(0) + dNdeta[1]*x_coords(1) + dNdeta[2]*x_coords(2) + dNdeta[3]*x_coords(3);
            J0(1,1) = dNdeta[0]*y_coords(0) + dNdeta[1]*y_coords(1) + dNdeta[2]*y_coords(2) + dNdeta[3]*y_coords(3);

            double detJ0 = J0.determinant();
            if (detJ0 < SMALL_VALUE_DOUBLE) { // Para evitar Jacobiano negativo ou muito pequeno
                std::cerr << "ERRO: Determinante do Jacobiano de referencia muito pequeno ou negativo no elemento " << element.id << " gp " << gp_idx << ": " << detJ0 << "\n";
                detJ0 = SMALL_VALUE_DOUBLE; // Força um valor positivo pequeno para continuar
            }

            Matrix J0_inv = J0.inverse();

            // Derivadas das funções de forma em relação a X e Y (no espaço material/referencial)
            double dNdX[4];
            double dNdY[4];
            for (int k = 0; k < 4; ++k) {
                dNdX[k] = J0_inv(0,0) * dNdxi[k] + J0_inv(0,1) * dNdeta[k];
                dNdY[k] = J0_inv(1,0) * dNdxi[k] + J0_inv(1,1) * dNdeta[k];
            }

            // Gradiente de Deslocamento (Grad_u) [du/dX, du/dY; dv/dX, dv/dY]
            Matrix Grad_u(2,2);
            Grad_u(0,0) = (dNdX[0]*u_coords(0) + dNdX[1]*u_coords(1) + dNdX[2]*u_coords(2) + dNdX[3]*u_coords(3));
            Grad_u(0,1) = (dNdY[0]*u_coords(0) + dNdY[1]*u_coords(1) + dNdY[2]*u_coords(2) + dNdY[3]*u_coords(3));
            Grad_u(1,0) = (dNdX[0]*v_coords(0) + dNdX[1]*v_coords(1) + dNdX[2]*v_coords(2) + dNdX[3]*v_coords(3));
            Grad_u(1,1) = (dNdY[0]*v_coords(0) + dNdY[1]*v_coords(1) + dNdY[2]*v_coords(2) + dNdY[3]*v_coords(3));

            // Gradiente de Deformação F = I + Grad_u (Referencial para Atual)
            Matrix F_total = Matrix::Identity(2, 2) + Grad_u;

            // Pega os dados do ponto de Gauss atual
            GaussPointData& current_gp_data = element.gp_data[gp_idx];
            Matrix C_ep_consistent; // Matriz constitutiva elastoplástica tangente consistente

            // Chamada para a função de atualização de tensão e cálculo da tangente
            plasticity_stress_update(F_total, current_gp_data, material_props, C_ep_consistent);

            // Determinante de F_total (J)
            double J_total = F_total.determinant();
            if (std::abs(J_total) < SMALL_VALUE_DOUBLE) {
                J_total = SMALL_VALUE_DOUBLE; // Evita divisão por zero
            }

            // Tensão de Cauchy (sigma_cauchy_matrix) no formato matricial
            Matrix sigma_cauchy_matrix(2,2);
            sigma_cauchy_matrix << current_gp_data.sigma_cauchy(0), current_gp_data.sigma_cauchy(2),
                                   current_gp_data.sigma_cauchy(2), current_gp_data.sigma_cauchy(1);

            // 2ª Tensão de Piola-Kirchhoff (S)
            Matrix S_matrix = PK2_from_Cauchy(F_total, sigma_cauchy_matrix);

            // Vetor S_vec [S_xx, S_yy, S_xy] (para cálculo das forças internas)
            Vector S_vec(3);
            S_vec(0) = S_matrix(0,0);
            S_vec(1) = S_matrix(1,1);
            S_vec(2) = S_matrix(0,1);

            // Matriz B_L (para a parte material da rigidez e forças internas)
            // Mapeia deslocamentos para deformações de Green-Lagrange
            Matrix B_L(3, 8);
            B_L.setZero();
            for (int k = 0; k < 4; ++k) {
                B_L(0, 2*k) = dNdX[k];
                B_L(1, 2*k+1) = dNdY[k];
                B_L(2, 2*k) = dNdY[k];
                B_L(2, 2*k+1) = dNdX[k];
            }

            K_element_material_out += B_L.transpose() * C_ep_consistent * B_L * detJ0 * w_xi * w_eta;

            // Matriz S_hat (para a parte geométrica da rigidez)
            Matrix S_hat(4,4);
            S_hat.setZero();
            S_hat << S_vec(0), S_vec(2), 0.0, 0.0,
                     S_vec(2), S_vec(1), 0.0, 0.0,
                     0.0, 0.0, S_vec(0), S_vec(2),
                     0.0, 0.0, S_vec(2), S_vec(1);

            // Matriz G_matrix (para a parte geométrica da rigidez)
            Matrix G_matrix(4, 8);
            G_matrix.setZero();
            for (int k = 0; k < 4; ++k) {
                G_matrix(0, 2*k) = dNdX[k];
                G_matrix(1, 2*k) = dNdY[k];
                G_matrix(2, 2*k+1) = dNdX[k];
                G_matrix(3, 2*k+1) = dNdY[k];
            }

            K_element_geometric_out += G_matrix.transpose() * S_hat * G_matrix * detJ0 * w_xi * w_eta;

            F_int_element += B_L.transpose() * S_vec * detJ0 * w_xi * w_eta;

            gp_idx++;
        }
    }
    // Calcular a deformação plástica média do elemento para critério de refinamento
    double sum_peeq = 0.0;
    for (const auto& gp : element.gp_data) {
        sum_peeq += gp.equivalent_plastic_strain;
    }
    element.avg_plastic_strain = sum_peeq / element.gp_data.size();

    return K_element_material_out + K_element_geometric_out;
}

// Calcula as tensões de Cauchy nos elementos para grandes deformações
void calculate_element_stresses_large_deformation(std::vector<Element>& elements) {
    for (auto& element : elements) {
        Vector avg_sigma_cauchy = Vector::Zero(3);
        for (const auto& gp : element.gp_data) {
            avg_sigma_cauchy += gp.sigma_cauchy;
        }
        avg_sigma_cauchy /= element.gp_data.size();

        element.stress_xx = Pa(avg_sigma_cauchy(0)); // Atribui com tipo de unidade
        element.stress_yy = Pa(avg_sigma_cauchy(1));
        element.stress_xy = Pa(avg_sigma_cauchy(2));
    }
}

// --- Funções de Refinamento e Remalhagem ---

// H-refinement: divide elementos em 4 sub-elementos se o critério for satisfeito
void refine_mesh_h_refinement(std::vector<Node>& nodes, std::vector<Element>& elements,
                               double refinement_threshold_peeq,
                               const MaterialProperties& material_props) {
    std::vector<Element> new_elements;
    std::map<std::pair<int, int>, int> edge_mid_nodes; // Para gerenciar nós do meio da aresta (para não duplicar)

    std::cout << "\nIniciando h-refinement com threshold PEEQ = " << refinement_threshold_peeq << "...\n";

    for (const auto& element : elements) {
        if (element.avg_plastic_strain > refinement_threshold_peeq) {
            // Este elemento será refinado
            // 1. Criar nós no meio das arestas e no centro do elemento
            // 2. Criar 4 novos sub-elementos

            // IDs dos nós do elemento original
            int n0_id = element.nodes[0];
            int n1_id = element.nodes[1];
            int n2_id = element.nodes[2];
            int n3_id = element.nodes[3];

            // Coordenadas dos nós originais
            Length n0_x = nodes[n0_id].x; Length n0_y = nodes[n0_id].y;
            Length n1_x = nodes[n1_id].x; Length n1_y = nodes[n1_id].y;
            Length n2_x = nodes[n2_id].x; Length n2_y = nodes[n2_id].y;
            Length n3_x = nodes[n3_id].x; Length n3_y = nodes[n3_id].y;

            // Funcao para obter ou criar no no meio de uma aresta
            auto get_or_create_mid_node = [&](int node_a_id, int node_b_id) -> int {
                // Ordena os IDs para garantir que (a,b) e (b,a) resultem no mesmo ID de aresta
                std::pair<int, int> edge_key = {std::min(node_a_id, node_b_id), std::max(node_a_id, node_b_id)};
                if (edge_mid_nodes.count(edge_key)) {
                    return edge_mid_nodes[edge_key];
                } else {
                    Node& node_a = nodes[node_a_id];
                    Node& node_b = nodes[node_b_id];
                    Length mid_x = (node_a.x + node_b.x) / 2.0;
                    Length mid_y = (node_a.y + node_b.y) / 2.0;

                    int new_node_id = next_node_id++;
                    nodes.emplace_back(new_node_id, mid_x, mid_y);
                    // Interpolar deslocamentos do nó do meio
                    nodes[new_node_id].ux = (node_a.ux + node_b.ux) / 2.0;
                    nodes[new_node_id].uy = (node_a.uy + node_b.uy) / 2.0;
                    edge_mid_nodes[edge_key] = new_node_id;
                    return new_node_id;
                }
            };

            // Nós do meio das arestas
            int n01_id = get_or_create_mid_node(n0_id, n1_id);
            int n12_id = get_or_create_mid_node(n1_id, n2_id);
            int n23_id = get_or_create_mid_node(n2_id, n3_id);
            int n30_id = get_or_create_mid_node(n3_id, n0_id);

            // Nó central do elemento
            Length center_x = (n0_x + n1_x + n2_x + n3_x) / 4.0;
            Length center_y = (n0_y + n1_y + n2_y + n3_y) / 4.0;
            int center_node_id = next_node_id++;
            nodes.emplace_back(center_node_id, center_x, center_y);
            // Interpolar deslocamentos do nó central
            nodes[center_node_id].ux = (nodes[n0_id].ux + nodes[n1_id].ux + nodes[n2_id].ux + nodes[n3_id].ux) / 4.0;
            nodes[center_node_id].uy = (nodes[n0_id].uy + nodes[n1_id].uy + nodes[n2_id].uy + nodes[n3_id].uy) / 4.0;

            // Criar 4 novos elementos
            // Elemento 1 (inferior esquerdo)
            new_elements.emplace_back(next_element_id++, n0_id, n01_id, center_node_id, n30_id);
            // Elemento 2 (inferior direito)
            new_elements.emplace_back(next_element_id++, n01_id, n1_id, n12_id, center_node_id);
            // Elemento 3 (superior direito)
            new_elements.emplace_back(next_element_id++, center_node_id, n12_id, n2_id, n23_id);
            // Elemento 4 (superior esquerdo)
            new_elements.emplace_back(next_element_id++, n30_id, center_node_id, n23_id, n3_id);

            // Interpolar GP data para os novos elementos
            // Para cada novo elemento, seus GPs serão inicializados com a interpolação dos GPs do elemento pai.
            // Para simplicidade, vamos usar os dados médios do elemento pai.
            // Uma interpolação mais rigorosa precisaria mapear os GPs do novo elemento para as coordenadas naturais do elemento pai,
            // e depois usar as funções de forma do elemento pai para interpolar.
            GaussPointData avg_gp_data_from_parent;
            for(const auto& gp : element.gp_data) {
                avg_gp_data_from_parent.sigma_cauchy += gp.sigma_cauchy;
                avg_gp_data_from_parent.F_elastic += gp.F_elastic;
                avg_gp_data_from_parent.F_plastic += gp.F_plastic;
                avg_gp_data_from_parent.equivalent_plastic_strain += gp.equivalent_plastic_strain;
                avg_gp_data_from_parent.current_yield_strength = Pa(avg_gp_data_from_parent.current_yield_strength.value + gp.current_yield_strength.value);
            }
            avg_gp_data_from_parent.sigma_cauchy /= element.gp_data.size();
            avg_gp_data_from_parent.F_elastic /= element.gp_data.size();
            avg_gp_data_from_parent.F_plastic /= element.gp_data.size();
            avg_gp_data_from_parent.equivalent_plastic_strain /= element.gp_data.size();
            avg_gp_data_from_parent.current_yield_strength = Pa(avg_gp_data_from_parent.current_yield_strength.value / element.gp_data.size());

            for (int k = 0; k < 4; ++k) {
                new_elements[new_elements.size() - 4 + k].gp_data.assign(4, avg_gp_data_from_parent);
                // Garante que a tensão de escoamento inicial para os novos GPs seja consistente
                for (auto& gp : new_elements[new_elements.size() - 4 + k].gp_data) {
                    gp.current_yield_strength = material_props.yield_strength; // Reinicia ou usa a do parent
                }
            }

        } else {
            // Elemento não será refinado, adiciona-o à nova lista
            new_elements.push_back(element);
        }
    }
    elements = new_elements;
    // Renumerar IDs de nós para garantir sequencialidade após adição
    // Isso é importante para evitar lacunas e para o mapeamento direto em arrays/vectors.
    std::map<int, int> old_to_new_node_id_map;
    std::vector<Node> renumbered_nodes;
    next_node_id = 0; // Resetar para renumerar

    for (auto& node : nodes) {
        old_to_new_node_id_map[node.id] = next_node_id;
        node.id = next_node_id++;
        renumbered_nodes.push_back(node);
    }
    nodes = renumbered_nodes; // Atualiza a lista de nós com os renumerados

    // Atualiza a conectividade dos elementos com os novos IDs de nós
    for (auto& element : elements) {
        for (int i = 0; i < 4; ++i) {
            element.nodes[i] = old_to_new_node_id_map[element.nodes[i]];
        }
    }

    std::cout << "h-refinement concluido. Nova malha: " << nodes.size() << " nos, " << elements.size() << " elementos.\n";
}


// Remalhagem global: cria uma nova malha regular e interpola dados da malha antiga
void remesh_global_mesh(std::vector<Node>& nodes, std::vector<Element>& elements,
                        Length mesh_width, Length mesh_height,
                        int new_divisions_x, int new_divisions_y,
                        const MaterialProperties& material_props) {

    std::cout << "\nIniciando remalhagem global...\n";

    std::vector<Node> old_nodes = nodes;
    std::vector<Element> old_elements = elements;

    // 1. Criar a nova malha regular
    nodes.clear();
    elements.clear();
    next_node_id = 0; // Resetar IDs
    next_element_id = 0;

    discretize_rectangle(mesh_width, mesh_height, new_divisions_x, new_divisions_y, nodes, elements);

    // 2. Interpolar os dados da malha antiga para a nova malha
    // Para cada ponto de Gauss em cada novo elemento, encontrar o elemento antigo correspondente
    // e interpolar os dados de GP.

    const double gp_coord = 1.0 / std::sqrt(3.0);
    const double gauss_points[4][2] = {
        {-gp_coord, -gp_coord},
        { gp_coord, -gp_coord},
        { gp_coord,  gp_coord},
        {-gp_coord,  gp_coord}
    };

    int elements_with_data_interp = 0;
    for (auto& new_element : elements) {
        // Para cada GP do novo elemento, encontrar suas coordenadas no espaço referencial
        // e, em seguida, as coordenadas no espaço deformado.

        // Coordenadas originais (referencial) dos nós do NOVO elemento
        Vector new_el_ref_x_coords(4), new_el_ref_y_coords(4);
        for(int i=0; i<4; ++i) {
            new_el_ref_x_coords(i) = nodes[new_element.nodes[i]].x.value;
            new_el_ref_y_coords(i) = nodes[new_element.nodes[i]].y.value;
        }

        bool element_interp_success = false;
        new_element.gp_data.resize(4); // Garante 4 GPs

        for (int gp_idx = 0; gp_idx < 4; ++gp_idx) {
            double xi = gauss_points[gp_idx][0];
            double eta = gauss_points[gp_idx][1];

            double N_val_ref[4], dNdxi_val_ref[4], dNdeta_val_ref[4];
            calculate_shape_functions_and_derivatives(xi, eta, N_val_ref, dNdxi_val_ref, dNdeta_val_ref);

            // Calcula a coordenada (X_ref, Y_ref) no espaço referencial
            double X_ref_gp = 0.0, Y_ref_gp = 0.0;
            for (int i = 0; i < 4; ++i) {
                X_ref_gp += N_val_ref[i] * new_el_ref_x_coords(i);
                Y_ref_gp += N_val_ref[i] * new_el_ref_y_coords(i);
            }

            // Precisamos das coordenadas ATUAIS (deformadas) no ponto de Gauss.
            // Para isso, precisamos do deslocamento no ponto de Gauss.
            // Os deslocamentos dos nós da NOVA malha são zero a princípio.
            // Então, precisamos interpolar o deslocamento (ux, uy) para o ponto (X_ref_gp, Y_ref_gp)
            // da MALHA ANTIGA.

            // Uma abordagem mais robusta é interpolar diretamente o GRADIENTE DE DEFORMAÇÃO (F)
            // e a DEFORMAÇÃO PLÁSTICA EQUIVALENTE (PEEQ), e depois computar o resto.
            // Para esta simplificação, vamos interpolar os dados dos GPs da malha antiga para os nós da malha antiga,
            // e depois para as coordenadas de interesse na nova malha.
            // Isso requer uma busca de ponto em elemento.

            Element* source_element = find_element_containing_point(old_elements, old_nodes, X_ref_gp, Y_ref_gp);

            if (source_element) {
                // Obter as coordenadas ATUAIS dos nós do elemento fonte
                std::vector<Vector> source_current_node_coords = get_nodal_current_coordinates_for_element(*source_element, old_nodes);

                // Interpolar os dados do GP do elemento fonte para o ponto de gauss do novo elemento.
                // Como não temos um método direto para interpolar dados de GP entre elementos (apenas para o centro ou nós),
                // e para manter a complexidade gerenciável, vamos interpolar os dados médios do elemento fonte
                // para os GPs do novo elemento.
                // Uma solução mais precisa exigiria um mapeamento inverso das funções de forma
                // para encontrar o xi,eta do ponto de destino no elemento fonte.
                new_element.gp_data[gp_idx] = interpolate_gp_data_to_point(*source_element, source_current_node_coords, X_ref_gp, Y_ref_gp);
                element_interp_success = true;
            } else {
                // Se o ponto de Gauss da nova malha não cair em nenhum elemento da malha antiga,
                // ele provavelmente está fora do domínio deformado coberto pelos elementos antigos.
                // Isso pode acontecer se a remalhagem ocorrer após grande distorção, ou se a nova malha for maior.
                // Inicializa com material virgem.
                new_element.gp_data[gp_idx].sigma_cauchy.setZero();
                new_element.gp_data[gp_idx].F_elastic = Matrix::Identity(2,2);
                new_element.gp_data[gp_idx].F_plastic = Matrix::Identity(2,2);
                new_element.gp_data[gp_idx].equivalent_plastic_strain = 0.0;
                new_element.gp_data[gp_idx].current_yield_strength = material_props.yield_strength;
            }
        }
        if (!element_interp_success) {
            // Marca o elemento para não ser usado se a interpolação falhar para todos os GPs.
            // Para simplicidade, assumimos que se source_element foi encontrado para alguns GPs, ele funciona.
            // Em um sistema real, você faria um tratamento mais rigoroso.
        } else {
            elements_with_data_interp++;
        }
    }

    // Interpolar os deslocamentos dos nós da nova malha
    for (auto& new_node : nodes) {
        Element* source_element = find_element_containing_point(old_elements, old_nodes, new_node.x.value, new_node.y.value);
        if (source_element) {
            std::vector<Vector> source_current_node_coords = get_nodal_current_coordinates_for_element(*source_element, old_nodes);
            Vector N_vec; Matrix B_matrix; double detJ; Vector natural_coords;
            if (calculate_N_and_B_at_point(source_current_node_coords, new_node.x.value, new_node.y.value,
                                           N_vec, B_matrix, detJ, natural_coords)) {
                // Interpolar deslocamentos dos nós do elemento fonte para o novo nó
                new_node.ux = m(0.0);
                new_node.uy = m(0.0);
                for (int i = 0; i < 4; ++i) {
                    int old_node_id = source_element->nodes[i];
                    new_node.ux = new_node.ux + nodes[old_node_id].ux * N_vec(i); // Usar os Ns do ponto de interesse
                    new_node.uy = new_node.uy + nodes[old_node_id].uy * N_vec(i);
                }
            }
        }
        // Se o nó não for encontrado (fora da malha antiga), seus deslocamentos permanecem zero (condição inicial)
    }

    std::cout << "Remalhagem global concluida. Nova malha: " << nodes.size() << " nos, " << elements.size() << " elementos.\n";
}


// --- Solver Newton-Raphson com Contato e Melhores Práticas ---
// Retorna o número de iterações do NR para uso no passo adaptativo
int solve_fem_system_with_contact(std::vector<Node>& nodes, std::vector<Element>& elements,
                                  const MaterialProperties& material_props,
                                  Force total_applied_value_force_mode, // Recebe como Force
                                  Length mesh_width, Length initial_height, // Recebe como Length
                                  int divisions_x, int divisions_y, double friction_coefficient, int load_type_choice,
                                  int fixed_ux_node_id, Length current_value_at_step_target_disp_mode,
                                  double thickness) { // Nova variável para espessura

    int num_dofs = nodes.size() * 2;

    Vector U_current_step_total_NR_iter(num_dofs);
    for (size_t i = 0; i < nodes.size(); ++i) {
        U_current_step_total_NR_iter(i * 2) = nodes[i].ux.value; // Acessa o valor double
        U_current_step_total_NR_iter(i * 2 + 1) = nodes[i].uy.value;
    }

    std::vector<int> previous_contact_states(nodes.size());
    std::vector<int> current_contact_states(nodes.size());

    for (size_t i = 0; i < nodes.size(); ++i) {
        previous_contact_states[i] = nodes[i].contact_state;
        current_contact_states[i] = nodes[i].contact_state;
    }

    Matrix K_el_mat_temp, K_el_geom_temp;
    Vector F_int_el_temp;

    int actual_nr_iterations = 0;

    for (int iter = 0; iter < MAX_NR_ITERATIONS; ++iter) {
        actual_nr_iterations = iter + 1;

        // Salva o estado dos nós no início da iteração do NR para o line search
        std::vector<Length> nodes_ux_prev_iter(nodes.size());
        std::vector<Length> nodes_uy_prev_iter(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            nodes_ux_prev_iter[i] = nodes[i].ux;
            nodes_uy_prev_iter[i] = nodes[i].uy;
        }

        // Atualiza os nós com o deslocamento total acumulado para a iteração atual do NR
        for (size_t i = 0; i < nodes.size(); ++i) {
            nodes[i].ux = m(U_current_step_total_NR_iter(i * 2)); // Atribui com o tipo Length
            nodes[i].uy = m(U_current_step_total_NR_iter(i * 2 + 1));
        }

        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(elements.size() * 8 * 8);

        Vector F_internal_global = Vector::Zero(num_dofs); // Forças em Newtons

        for (auto& element : elements) {
            Matrix Kt_element = calculate_element_tangent_stiffness_and_internal_forces(element, nodes, material_props,
                                                                                       F_int_el_temp, K_el_mat_temp, K_el_geom_temp);

            for (int i_node = 0; i_node < 4; ++i_node) {
                for (int j_node = 0; j_node < 4; ++j_node) {
                    int global_row_x = element.nodes[i_node] * 2;
                    int global_row_y = element.nodes[i_node] * 2 + 1;
                    int global_col_x = element.nodes[j_node] * 2;
                    int global_col_y = element.nodes[j_node] * 2 + 1;

                    // Multiplica por 'thickness' para converter rigidez por unidade de comprimento para rigidez total
                    tripletList.emplace_back(global_row_x, global_col_x, Kt_element(2*i_node, 2*j_node) * thickness);
                    tripletList.emplace_back(global_row_x, global_col_y, Kt_element(2*i_node, 2*j_node+1) * thickness);
                    tripletList.emplace_back(global_row_y, global_col_x, Kt_element(2*i_node+1, 2*j_node) * thickness);
                    tripletList.emplace_back(global_row_y, global_col_y, Kt_element(2*i_node+1, 2*j_node+1) * thickness);
                }
            }
            for (int k_node = 0; k_node < 4; ++k_node) {
                F_internal_global(element.nodes[k_node] * 2) += F_int_el_temp(2*k_node) * thickness;
                F_internal_global(element.nodes[k_node] * 2 + 1) += F_int_el_temp(2*k_node + 1) * thickness;
            }
        }
        SparseMatrix K_tangent_global(num_dofs, num_dofs);
        K_tangent_global.setFromTriplets(tripletList.begin(), tripletList.end());


        Vector F_external_applied = Vector::Zero(num_dofs); // Forças em Newtons
        if (load_type_choice == 1) { // Carregamento por Força
            // Forca_por_no: N / (unidade) = N
            double force_per_node = total_applied_value_force_mode.value / (static_cast<double>(divisions_x) + 1.0); // Corrigido para double
            // A força é aplicada nos nós do topo
            for (const auto& node : nodes) {
                if (std::abs(node.y.value - initial_height.value) < SMALL_VALUE_DOUBLE) {
                    F_external_applied(node.id * 2 + 1) = -force_per_node; // Força para baixo (compressão)
                }
            }
        }

        bool any_contact_state_changed_in_this_iter = false;

        // Aplica as condições de contorno e contato
        for (auto& node : nodes) {
            // Restrição de UX para o nó fixo no centro da base
            if (node.id == fixed_ux_node_id) {
                // Impõe ux = 0.0
                K_tangent_global.coeffRef(node.id * 2, node.id * 2) = 1.0;
                for (int k = 0; k < num_dofs; ++k) {
                    if (k != node.id * 2) {
                        K_tangent_global.coeffRef(node.id * 2, k) = 0.0;
                        K_tangent_global.coeffRef(k, node.id * 2) = 0.0;
                    }
                }
                // Ajusta o vetor de forças para refletir o deslocamento prescrito (0.0)
                F_external_applied(node.id * 2) = -(node.ux.value);
            }

            // Condicao de contato na base (y=0)
            if ((node.y.value + node.uy.value) < CONTACT_GAP_TOLERANCE.value) { // Compara valores double
                double internal_fy = F_internal_global(node.id * 2 + 1);

                if (internal_fy < -TOLERANCE_PRESSURE.value) { // Tendência a descolar
                    current_contact_states[node.id] = 0; // FREE
                    if (previous_contact_states[node.id] != 0) {
                        any_contact_state_changed_in_this_iter = true;
                    }
                } else { // Em contato ou querendo entrar em contato
                    // Impõe restrição de deslocamento em Y (normal)
                    K_tangent_global.coeffRef(node.id * 2 + 1, node.id * 2 + 1) = 1.0;
                    for (int k = 0; k < num_dofs; ++k) {
                        if (k != node.id * 2 + 1) {
                            K_tangent_global.coeffRef(node.id * 2 + 1, k) = 0.0;
                            K_tangent_global.coeffRef(k, node.id * 2 + 1) = 0.0;
                        }
                    }
                    F_external_applied(node.id * 2 + 1) = -(node.y.value + node.uy.value);


                    double tangential_force_tendency = F_internal_global(node.id * 2);
                    double normal_force_magnitude = std::abs(nodes[node.id].fy_reaction.value);

                    if (normal_force_magnitude < SMALL_VALUE_DOUBLE && iter == 0) { // Se ainda não calculada, estima
                        normal_force_magnitude = std::max(0.0, F_internal_global(node.id * 2 + 1));
                    }
                    double max_friction_force = friction_coefficient * normal_force_magnitude;

                    if (std::abs(tangential_force_tendency) > max_friction_force + TOLERANCE_PRESSURE.value) {
                        // SLIP (Deslizando)
                        current_contact_states[node.id] = 2;
                        F_external_applied(node.id * 2) = -max_friction_force * (tangential_force_tendency / (std::abs(tangential_force_tendency) + SMALL_VALUE_DOUBLE));

                        K_tangent_global.coeffRef(node.id * 2, node.id * 2) = 0.0;
                        for (int k = 0; k < num_dofs; ++k) {
                            if (k != node.id * 2) {
                                K_tangent_global.coeffRef(node.id * 2, k) = 0.0;
                                K_tangent_global.coeffRef(k, node.id * 2) = 0.0;
                            }
                        }

                        if (previous_contact_states[node.id] != 2) {
                            any_contact_state_changed_in_this_iter = true;
                        }
                    } else {
                        // STICK (Grudado)
                        current_contact_states[node.id] = 1;
                        K_tangent_global.coeffRef(node.id * 2, node.id * 2) = 1.0;
                        for (int k = 0; k < num_dofs; ++k) {
                            if (k != node.id * 2) {
                                K_tangent_global.coeffRef(node.id * 2, k) = 0.0;
                                K_tangent_global.coeffRef(k, node.id * 2) = 0.0;
                            }
                        }
                        F_external_applied(node.id * 2) = -(node.ux.value);
                        if (previous_contact_states[node.id] != 1) {
                            any_contact_state_changed_in_this_iter = true;
                        }
                    }
                }
            } else { // Nó está acima da superfície de contato
                current_contact_states[node.id] = 0; // FREE
                if (previous_contact_states[node.id] != 0) {
                    any_contact_state_changed_in_this_iter = true;
                }
            }
            node.contact_state = current_contact_states[node.id];
        }

        // Condicao para carregamento por deslocamento prescrito no topo (load_type_choice == 2)
        if (load_type_choice == 2) {
            for (const auto& node : nodes) {
                if (std::abs(node.y.value - initial_height.value) < SMALL_VALUE_DOUBLE) { // Nó está no topo
                    K_tangent_global.coeffRef(node.id * 2 + 1, node.id * 2 + 1) = 1.0;
                    for (int k = 0; k < num_dofs; ++k) {
                        if (k != node.id * 2 + 1) {
                            K_tangent_global.coeffRef(node.id * 2 + 1, k) = 0.0;
                            K_tangent_global.coeffRef(k, node.id * 2 + 1) = 0.0;
                        }
                    }
                    F_external_applied(node.id * 2 + 1) = -(node.y.value + node.uy.value - current_value_at_step_target_disp_mode.value);
                }
            }
        }

        Vector F_residual = F_external_applied - F_internal_global;

        double residual_norm = F_residual.norm();
        double max_delta_U_norm = 0.0;

        Vector delta_U;

        Eigen::BiCGSTAB<SparseMatrix> solver;
        solver.setTolerance(TOLERANCE_PRESSURE.value);
        solver.setMaxIterations(MAX_NR_ITERATIONS);
        solver.compute(K_tangent_global);
        if(solver.info() != Eigen::Success) {
            std::cerr << "AVISO: Fatoracao do solver iterativo falhou. Retornando MAX_NR_ITERATIONS.\n";
            return MAX_NR_ITERATIONS;
        }
        delta_U = solver.solve(F_residual);
        if(solver.info() != Eigen::Success) {
            std::cerr << "AVISO: Solucao do sistema linear falhou. Retornando MAX_NR_ITERATIONS.\n";
            return MAX_NR_ITERATIONS;
        }

        // --- Line Search ---
        double alpha = 1.0;
        double min_alpha = 0.1;
        double acceptable_decrease_ratio = 0.8;

        auto calculate_residual_norm_for_alpha = [&](double current_alpha) {
            Vector temp_U = U_current_step_total_NR_iter + current_alpha * delta_U;

            std::vector<Node> temp_nodes = nodes;
            for(size_t i = 0; i < temp_nodes.size(); ++i) {
                temp_nodes[i].ux = m(temp_U(i * 2));
                temp_nodes[i].uy = m(temp_U(i * 2 + 1));
            }

            Vector temp_F_internal = Vector::Zero(num_dofs);
            Matrix K_el_mat_temp_ls, K_el_geom_temp_ls;

            for (auto& element_temp : elements) {
                Vector F_int_el_temp_ls;
                calculate_element_tangent_stiffness_and_internal_forces(element_temp, temp_nodes, material_props, F_int_el_temp_ls, K_el_mat_temp_ls, K_el_geom_temp_ls);
                for (int k_node = 0; k_node < 4; ++k_node) {
                    temp_F_internal(element_temp.nodes[k_node] * 2) += F_int_el_temp_ls(2*k_node) * thickness;
                    temp_F_internal(element_temp.nodes[k_node] * 2 + 1) += F_int_el_temp_ls(2*k_node + 1) * thickness;
                }
            }
            Vector temp_F_residual = F_external_applied - temp_F_internal;
            return temp_F_residual.norm();
        };

        double initial_residual_norm = F_residual.norm();
        double current_residual_norm_ls = initial_residual_norm;

        for (int ls_iter = 0; ls_iter < 10; ++ls_iter) {
            current_residual_norm_ls = calculate_residual_norm_for_alpha(alpha);
            if (current_residual_norm_ls <= initial_residual_norm * acceptable_decrease_ratio) {
                break;
            }
            alpha *= decrease_factor;
            if (alpha < min_alpha) {
                alpha = min_alpha;
                break;
            }
        }

        U_current_step_total_NR_iter += alpha * delta_U;

        max_delta_U_norm = (alpha * delta_U).norm();


        std::cout << "  NR Iter " << iter + 1 << ": Delta_U_norm = " << std::scientific << max_delta_U_norm
                  << ", Res_Norm = " << std::scientific << residual_norm << ", Alpha = " << std::scientific << alpha;
        if (any_contact_state_changed_in_this_iter) {
            std::cout << " (Estado de contato alterado)";
        }
        std::cout << "\n";

        if (max_delta_U_norm < TOLERANCE_DISPLACEMENT.value && residual_norm < TOLERANCE_PRESSURE.value && !any_contact_state_changed_in_this_iter) {
            std::cout << "Convergencia do Newton-Raphson alcancada em " << iter + 1 << " iteracoes.\n";
            break;
        }

        previous_contact_states = current_contact_states;

        if (iter == MAX_NR_ITERATIONS - 1) {
            std::cerr << "Aviso: Newton-Raphson NAO convergiu apos " << MAX_NR_ITERATIONS << " iteracoes. Max Delta U: " << std::scientific << max_delta_U_norm << ".\n";
            std::cerr << "  (Pode indicar problema de modelagem, tolerancia muito baixa ou numero insuficiente de iteracoes, ou problemas de convergencia de contato).\n";
        }
    }
    return actual_nr_iterations;
}

// Funcao para calcular as forcas de reacao nos nos que possuem restricoes
void calculate_nodal_reaction_forces(std::vector<Node>& nodes, std::vector<Element>& elements,
                                     const MaterialProperties& material_props,
                                     int divisions_x, int divisions_y, int fixed_ux_node_id,
                                     double friction_coefficient, int load_type_choice, Length initial_height,
                                     Length mesh_width, double thickness) {

    int num_dofs = nodes.size() * 2;
    Vector F_internal_global_for_reactions = Vector::Zero(num_dofs);

    for (auto& node : nodes) {
        node.fx_reaction = N(0.0);
        node.fy_reaction = N(0.0);
    }

    Matrix K_el_mat_dummy, K_el_geom_dummy;
    Vector F_int_el_temp_reactions;

    for (auto& element : elements) {
        calculate_element_tangent_stiffness_and_internal_forces(element, nodes, material_props,
                                                                 F_int_el_temp_reactions, K_el_mat_dummy, K_el_geom_dummy);

        for (int k = 0; k < 4; ++k) {
            F_internal_global_for_reactions(element.nodes[k] * 2) += F_int_el_temp_reactions(2*k) * thickness;
            F_internal_global_for_reactions(element.nodes[k] * 2 + 1) += F_int_el_temp_reactions(2*k + 1) * thickness;
        }
    }

    for (int i = 0; i < num_dofs; ++i) {
        int node_id = i / 2;
        bool is_x_dof = (i % 2 == 0);
        bool is_y_dof = (i % 2 == 1);

        if (is_y_dof && nodes[node_id].contact_state != 0 && std::abs(nodes[node_id].y.value) < SMALL_VALUE_DOUBLE) {
            nodes[node_id].fy_reaction = N(-F_internal_global_for_reactions(i));

            if (nodes[node_id].contact_state == 1) { // STICK
                nodes[node_id].fx_reaction = N(-F_internal_global_for_reactions(i-1));
            } else if (nodes[node_id].contact_state == 2) { // SLIP
                double normal_force_magnitude = std::abs(nodes[node_id].fy_reaction.value);
                double sign_tangential_force_tendency = (std::abs(F_internal_global_for_reactions(i-1)) > SMALL_VALUE_DOUBLE) ? (F_internal_global_for_reactions(i-1) / std::abs(F_internal_global_for_reactions(i-1))) : 0.0;
                nodes[node_id].fx_reaction = N(-friction_coefficient * normal_force_magnitude * sign_tangential_force_tendency);
            }
        }

        if (is_x_dof && node_id == fixed_ux_node_id) {
            nodes[node_id].fx_reaction = N(-F_internal_global_for_reactions(i));
        }

        if (load_type_choice == 2 && is_y_dof && std::abs(nodes[node_id].y.value - initial_height.value) < SMALL_VALUE_DOUBLE) {
            nodes[node_id].fy_reaction = N(-F_internal_global_for_reactions(i));
        }
    }
}

// --- Exportacao para VTU (ParaView) ---
void export_to_vtu(const std::string& filename, const std::vector<Node>& nodes, const std::vector<Element>& elements) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Erro: Nao foi possivel abrir o arquivo " << filename << " para escrita.\n";
        return;
    }

    ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    ofs << "  <UnstructuredGrid>\n";
    ofs << "    <Piece NumberOfPoints=\"" << nodes.size() << "\" NumberOfCells=\"" << elements.size() << "\">\n";

    ofs << "      <Points>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        ofs << std::fixed << std::setprecision(10) << node.current_x().value << " " // Exporta as coordenadas DEFORMADAS
            << std::fixed << std::setprecision(10) << node.current_y().value << " "
            << "0.0\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </Points>\n";

    ofs << "      <Cells>\n";
    ofs << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : elements) {
        ofs << element.nodes[0] << " " << element.nodes[1] << " " << element.nodes[2] << " " << element.nodes[3] << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 4;
    for (size_t i = 0; i < elements.size(); ++i) {
        ofs << offset << "\n";
        offset += 4;
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        ofs << "9\n"; // VTK_QUAD = 9
    }
    ofs << "        </DataArray>\n";
    ofs << "      </Cells>\n";

    ofs << "      <PointData>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Deslocamento\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        ofs << std::fixed << std::setprecision(10) << node.ux.value << " "
            << std::fixed << std::setprecision(10) << node.uy.value << " "
            << "0.0\n";
    }
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"Float64\" Name=\"Forca_X\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        ofs << std::fixed << std::setprecision(10) << node.fx_reaction.value << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Forca_Y\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        ofs << std::fixed << std::setprecision(10) << node.fy_reaction.value << "\n";
    }
    ofs << "        </DataArray>\n";

    ofs << "        <DataArray type=\"Int32\" Name=\"Contact_State\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& node : nodes) {
        ofs << node.contact_state << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </PointData>\n";

    ofs << "      <CellData>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Tensao_X\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& element : elements) {
        ofs << std::fixed << std::setprecision(10) << element.stress_xx.value << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Tensao_Y\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& element : elements) {
        ofs << std::fixed << std::setprecision(10) << element.stress_yy.value << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Tensao_Cisalhamento_XY\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& element : elements) {
        ofs << std::fixed << std::setprecision(10) << element.stress_xy.value << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"Deformacao_Plastica_Equivalente\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const auto& element : elements) {
        double total_peeq = 0.0;
        for (const auto& gp : element.gp_data) {
            total_peeq += gp.equivalent_plastic_strain;
        }
        ofs << std::fixed << std::setprecision(10) << total_peeq / element.gp_data.size() << "\n"; // Exporta a média dos GPs
    }
    ofs << "        </DataArray>\n";
    ofs << "      </CellData>\n";

    ofs << "    </Piece>\n";
    ofs << "  </UnstructuredGrid>\n";
    ofs << "</VTKFile>\n";
    ofs.close();
}

// Funcao para Simular e Exportar um Passo
void simulate_and_export_step(std::vector<Node>& nodes, std::vector<Element>& elements,
                              const MaterialProperties& material_props, // Ordem dos parâmetros ajustada
                              Force total_applied_value_force_mode, // Recebe Force
                              Length width, Length height_initial, // Recebe Length
                              int divisions_x, int divisions_y, double friction_coefficient, int choice,
                              int fixed_ux_node_id, Length current_value_at_step_target_disp_mode, double thickness,
                              int& nr_iterations_for_next_step_adjustment) { // Retorno por referência

    int nr_iterations_taken = solve_fem_system_with_contact(nodes, elements, material_props,
                                                             total_applied_value_force_mode, width, height_initial,
                                                             divisions_x, divisions_y, friction_coefficient, choice,
                                                             fixed_ux_node_id,
                                                             current_value_at_step_target_disp_mode, thickness);

    nr_iterations_for_next_step_adjustment = nr_iterations_taken;

    std::cout << "Deslocamentos apos o passo " << step_number << ":\n";
    for (const auto& node : nodes) {
        if (std::abs(node.y.value) < SMALL_VALUE_DOUBLE || std::abs(node.y.value - height_initial.value) < SMALL_VALUE_DOUBLE) {
            std::cout << "  Node " << node.id << " (x=" << node.x.value << ", y=" << node.y.value << "): ux = "
                      << std::scientific << node.ux.value << ", uy = " << std::scientific << node.uy.value
                      << ", Contact=" << node.contact_state << "\n";
        }
    }
    std::cout << "-------------------------------------------\n";

    calculate_element_stresses_large_deformation(elements);
    calculate_nodal_reaction_forces(nodes, elements, material_props, divisions_x, divisions_y, fixed_ux_node_id,
                                     friction_coefficient, choice, height_initial, width, thickness);

    std::stringstream ss;
    ss << "compressao_step_" << std::setw(4) << std::setfill('0') << step_number << ".vtu";
    export_to_vtu(ss.str(), nodes, elements);
}

// Funcao para exportar dados de tensao-deformacao para um arquivo TXT
void export_stress_strain_data(const std::string& filename, const std::vector<StressStrainPoint>& data) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Erro: Nao foi possivel abrir o arquivo " << filename << " para escrita dos dados de tensao-deformacao.\n";
        return;
    }

    ofs << "Deformacao_Engenharia,Tensao_Cauchy\n";
    for (const auto& point : data) {
        ofs << std::fixed << std::setprecision(10) << point.strain.value << ","
            << std::fixed << std::setprecision(10) << point.stress.value << "\n";
    }
    ofs.close();
    std::cout << "\nDados de tensao-deformacao exportados para " << filename << "\n";
}

// ==============================================================================
// === FUNÇÃO PRINCIPAL (main) ===
// ==============================================================================
int main() {
    std::vector<Node> nodes;
    std::vector<Element> elements;

    Length width = m(0.0);
    Length height_initial = m(0.0);
    int divisions_x = 0;
    int divisions_y = 0;
    Temperature process_temperature = C(0.0);
    Force total_applied_value_force_mode = N(0.0);
    Length total_applied_value_disp_mode = m(0.0);
    int initial_num_steps = 0;
    double friction_coefficient = 0.0;
    int choice = 0;
    double thickness = 1.0;

    // Parametros para refinamento e remalhagem
    double refinement_peeq_threshold = 0.1; // Limite de PEEQ para h-refinement (adimensional)
    int remesh_interval_steps = 10; // Remalhar a cada N passos
    int remesh_new_divisions_x = 0; // Novas divisões para remalhagem
    int remesh_new_divisions_y = 0;

    std::vector<StressStrainPoint> stress_strain_curve;

    std::cout << "--- Parametros da Malha ---\n";
    std::cout << "Digite a largura da malha (m): ";
    double input_width_val;
    std::cin >> input_width_val;
    width = m(input_width_val);

    std::cout << "Digite a altura inicial da malha (m): ";
    double input_height_val;
    std::cin >> input_height_val;
    height_initial = m(input_height_val);

    std::cout << "Digite a espessura da peca (m, para analise 2D assume-se 1.0 se nao especificado): ";
    std::cin >> thickness;
    if (thickness <= SMALL_VALUE_DOUBLE) {
        std::cerr << "Aviso: Espessura nao pode ser zero ou negativa. Usando 1.0m.\n";
        thickness = 1.0;
    }

    std::cout << "Digite o numero de divisoes na direcao X: ";
    std::cin >> divisions_x;

    std::cout << "Digite o numero de divisoes na direcao Y: ";
    std::cin >> divisions_y;

    std::cout << "\n--- Condicoes de Processo ---\n";
    std::cout << "Digite a temperatura de processo em Celsius (ex: 20 para frio, 800 para quente): ";
    double input_temp_val;
    std::cin >> input_temp_val;
    process_temperature = C(input_temp_val);

    if (process_temperature.value < 0.0) {
        std::cerr << "Aviso: Temperatura de processo nao pode ser negativa. Usando 0 Celsius.\n";
        process_temperature = C(0.0);
    }

    std::string process_regime_str;
    if (process_temperature.value >= HOT_PROCESS_TEMPERATURE_THRESHOLD.value) {
        process_regime_str = "A Quente";
    } else {
        process_regime_str = "A Frio";
    }

    MaterialProperties current_material_properties = get_material_properties_by_temperature(process_temperature);

    std::cout << "\n--- Modo de Carregamento ---\n";
    std::cout << "Deseja simular por:\n";
    std::cout << "1. Forca de compressao aplicada no topo (sentido Y negativo)\n";
    std::cout << "2. Altura final desejada\n";
    std::cout << "Escolha (1 ou 2): ";
    std::cin >> choice;

    while (choice != 1 && choice != 2) {
        std::cout << "Escolha invalida. Por favor, digite 1 ou 2: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin >> choice;
    }

    if (choice == 1) {
        std::cout << "Digite a forca de compressao TOTAL a ser aplicada (N). Um valor POSITIVO indica compressao: ";
        double input_force_val;
        std::cin >> input_force_val;
        total_applied_value_force_mode = N(input_force_val);
        if (total_applied_value_force_mode.value < 0.0) {
            std::cerr << "Aviso: Forca de compressao negativa inserida. Convertendo para valor positivo.\n";
            total_applied_value_force_mode = N(std::abs(total_applied_value_force_mode.value));
        }
    } else {
        std::cout << "Digite a altura final DESEJADA da malha (m): ";
        double input_final_height_val;
        std::cin >> input_final_height_val;
        total_applied_value_disp_mode = m(input_final_height_val);
        if (total_applied_value_disp_mode.value > height_initial.value) {
            std::cerr << "Aviso: A altura final desejada (" << total_applied_value_disp_mode.value << "m) e maior que a altura inicial (" << height_initial.value << "m). Isso resultara em tracao. Ajustando para altura inicial.\n";
            total_applied_value_disp_mode = height_initial;
        }
        if (total_applied_value_disp_mode.value < SMALL_VALUE_DOUBLE) {
            std::cerr << "Aviso: A altura final nao pode ser zero ou negativa. Ajustando para um valor pequeno acima de zero.\n";
            total_applied_value_disp_mode = m(SMALL_VALUE_DOUBLE);
        }
    }

    std::cout << "\nDigite o COEFICIENTE DE ATRITO (mu) nas superficies superior e inferior (0.0 para sem atrito, ex: 0.2, 0.5): ";
    std::cin >> friction_coefficient;

    if (friction_coefficient < 0.0) {
        std::cerr << "O coeficiente de atrito nao pode ser negativo. Usando 0.0.\n";
        friction_coefficient = 0.0;
    }

    std::cout << "Digite o numero INICIAL de passos de simulacao (sera ajustado adaptativamente, ex: 50): ";
    std::cin >> initial_num_steps;
    if (initial_num_steps <= 0) {
        std::cerr << "O numero de passos deve ser positivo. Usando 50 passos.\n";
        initial_num_steps = 50;
    }

    std::cout << "\n--- Configuracoes de Refinamento e Remalhagem ---\n";
    std::cout << "Deseja ativar h-refinement (dividir elementos com alta deformacao)? (1=Sim, 0=Nao): ";
    int enable_h_refinement_input;
    std::cin >> enable_h_refinement_input;
    bool enable_h_refinement = (enable_h_refinement_input == 1);

    if (enable_h_refinement) {
        std::cout << "  Digite o limite de Deformacao Plastica Equivalente para h-refinement (ex: 0.1 para 10%): ";
        std::cin >> refinement_peeq_threshold;
        if (refinement_peeq_threshold <= 0.0) {
            std::cerr << "  Aviso: Limite de PEEQ invalido. Usando 0.1.\n";
            refinement_peeq_threshold = 0.1;
        }
    }

    std::cout << "Deseja ativar remalhagem global periodica? (1=Sim, 0=Nao): ";
    int enable_remeshing_input;
    std::cin >> enable_remeshing_input;
    bool enable_remeshing = (enable_remeshing_input == 1);

    if (enable_remeshing) {
        std::cout << "  Remalhar a cada quantos passos? (ex: 10): ";
        std::cin >> remesh_interval_steps;
        if (remesh_interval_steps <= 0) {
            std::cerr << "  Aviso: Intervalo de remalhagem invalido. Usando 10.\n";
            remesh_interval_steps = 10;
        }
        std::cout << "  Nova malha: Numero de divisoes em X para remalhagem (ex: 40): ";
        std::cin >> remesh_new_divisions_x;
        std::cout << "  Nova malha: Numero de divisoes em Y para remalhagem (ex: 40): ";
        std::cin >> remesh_new_divisions_y;
        if (remesh_new_divisions_x <= 0 || remesh_new_divisions_y <= 0) {
            std::cerr << "  Aviso: Divisoes invalidas para remalhagem. Usando 40x40.\n";
            remesh_new_divisions_x = 40;
            remesh_new_divisions_y = 40;
        }
    }

    // --- Resumo dos Parametros ---
    std::cout << "\n--- Resumo dos Parametros ---\n";
    std::cout << "Largura da malha: " << width.value << "m\n";
    std::cout << "Altura inicial da malha: " << height_initial.value << "m\n";
    std::cout << "Espessura da peca: " << thickness << "m\n";
    std::cout << "Divisoes iniciais: " << divisions_x << "x" << divisions_y << "\n";
    std::cout << "Temperatura de Processo: " << process_temperature.value << " °C (" << process_regime_str << ")\n";
    std::cout << "Material: E=" << std::scientific << current_material_properties.youngs_modulus.value << " Pa, nu=" << current_material_properties.poisson_ratio
              << ", Yield=" << current_material_properties.yield_strength.value << " Pa, K=" << current_material_properties.strength_coefficient.value
              << " Pa, n=" << current_material_properties.hardening_exponent << "\n";
    std::cout << "Coeficiente de Atrito (mu): " << friction_coefficient << "\n";
    std::cout << "Numero INICIAL de passos para a animacao: " << initial_num_steps << "\n";

    if (choice == 1) {
        std::cout << "Forca de Compressao TOTAL Aplicada: " << total_applied_value_force_mode.value << " N\n";
    } else {
        std::cout << "Altura Final DESEJADA: " << total_applied_value_disp_mode.value << " m\n";
    }
    std::cout << "H-Refinement: " << (enable_h_refinement ? "Ativado" : "Desativado") << (enable_h_refinement ? " (Limite PEEQ: " + std::to_string(refinement_peeq_threshold) + ")" : "") << "\n";
    std::cout << "Remalhagem Global: " << (enable_remeshing ? "Ativado" : "Desativado") << (enable_remeshing ? " (A cada " + std::to_string(remesh_interval_steps) + " passos, para " + std::to_string(remesh_new_divisions_x) + "x" + std::to_string(remesh_new_divisions_y) + " divisoes)" : "") << "\n\n";


    // --- Discretização da Malha Inicial ---
    discretize_rectangle(width, height_initial, divisions_x, divisions_y, nodes, elements);

    // --- Inicialização dos Dados dos Pontos de Gauss para PLASTICIDADE ---
    for (auto& element : elements) {
        for (auto& gp : element.gp_data) {
            gp.current_yield_strength = current_material_properties.yield_strength;
            gp.F_elastic = Matrix::Identity(2,2);
            gp.F_plastic = Matrix::Identity(2,2);
        }
    }

    stress_strain_curve.push_back({dimensionless(0.0), Pa(0.0)});

    // --- Exportação do Estado Inicial (Passo 0) ---
    std::cout << "Exportando estado inicial (Passo 0)....\n";
    for (auto& node : nodes) {
        node.reset_displacements_and_forces();
        node.contact_state = 0;
    }

    int fixed_ux_node_id_main = -1;
    Length center_x_main = width / 2.0;
    double min_dist_x_main = std::numeric_limits<double>::max();

    for (const auto& node : nodes) {
        if (std::abs(node.y.value) < SMALL_VALUE_DOUBLE) {
            double current_dist_x = std::abs(node.x.value - center_x_main.value);
            if (current_dist_x < min_dist_x_main) {
                min_dist_x_main = current_dist_x;
                fixed_ux_node_id_main = node.id;
            }
        }
    }

    // Atualiza a tensão de escoamento nos GPs para as propriedades do material
    for (auto& element : elements) {
        for (auto& gp : element.gp_data) {
            gp.current_yield_strength = current_material_properties.yield_strength;
        }
    }

    calculate_element_stresses_large_deformation(elements);
    calculate_nodal_reaction_forces(nodes, elements, current_material_properties, divisions_x, divisions_y, fixed_ux_node_id_main,
                                     friction_coefficient, choice, height_initial, width, thickness);
    step_number = 0;
    export_to_vtu("compressao_step_0000.vtu", nodes, elements);

    // --- Loop Principal da Simulação Incremental com Passo Adaptativo ---
    std::cout << "Iniciando simulacao incremental com " << initial_num_steps << " passos iniciais.\n";

    double current_load_factor = 0.0;
    double total_load_target = 1.0;
    int last_nr_iterations = (int)TARGET_NR_ITERATIONS;

    while (current_load_factor < total_load_target - SMALL_VALUE_DOUBLE) {
        step_number++;

        double delta_load_factor_raw = (total_load_target - current_load_factor) / (initial_num_steps - step_number + 1.0);
        double adjustment_factor = TARGET_NR_ITERATIONS / static_cast<double>(last_nr_iterations);
        adjustment_factor = std::min(MAX_STEP_SIZE_FACTOR, std::max(MIN_STEP_SIZE_FACTOR, adjustment_factor));
        double delta_load_factor = delta_load_factor_raw * adjustment_factor;

        double next_load_factor = current_load_factor + delta_load_factor;
        if (next_load_factor > total_load_target) {
            next_load_factor = total_load_target;
        }

        Force current_force_at_step_target = N(0.0);
        Length current_height_at_step_target = m(0.0);

        if (choice == 1) {
            current_force_at_step_target = total_applied_value_force_mode * next_load_factor;
            std::cout << "Passo " << step_number << ": Fator de Carga = " << std::fixed << std::setprecision(4) << next_load_factor * 100.0 << "%, Forca = "
                      << std::fixed << std::setprecision(2) << current_force_at_step_target.value << " N\n";
        } else {
            current_height_at_step_target = height_initial + (total_applied_value_disp_mode - height_initial) * next_load_factor;
            std::cout << "Passo " << step_number << ": Fator de Carga = " << std::fixed << std::setprecision(4) << next_load_factor * 100.0 << "%, Altura alvo = "
                      << std::fixed << std::setprecision(6) << current_height_at_step_target.value << " m\n";
        }

        simulate_and_export_step(nodes, elements, current_material_properties,
                                 current_force_at_step_target, width, height_initial,
                                 divisions_x, divisions_y, friction_coefficient, choice,
                                 fixed_ux_node_id_main, current_height_at_step_target, thickness,
                                 last_nr_iterations);

        // --- Verificação e Aplicação de Refinamento/Remalhagem ---
        if (enable_h_refinement) {
            // H-refinement: verifica a cada passo (ou a cada N passos se preferir)
            // Note: O h-refinement atual pode causar elementos com geometrias ruins se feito repetidamente
            // sem remalhagem intermediária. É mais adequado para capturar gradientes em certas regiões.
            refine_mesh_h_refinement(nodes, elements, refinement_peeq_threshold, current_material_properties);
            // Após h-refinement, o número de nós/dofs muda, então o fixed_ux_node_id_main pode ser inválido.
            // É preciso re-localizá-lo ou garantir que o ID é persistente ou se refere a uma propriedade.
            // Para esta implementação, vamos re-localizar o nó central da base.
            fixed_ux_node_id_main = -1;
            min_dist_x_main = std::numeric_limits<double>::max();
            for (const auto& node : nodes) {
                if (std::abs(node.y.value) < SMALL_VALUE_DOUBLE) {
                    double current_dist_x = std::abs(node.x.value - center_x_main.value);
                    if (current_dist_x < min_dist_x_main) {
                        min_dist_x_main = current_dist_x;
                        fixed_ux_node_id_main = node.id;
                    }
                }
            }
        }

        if (enable_remeshing && (step_number % remesh_interval_steps == 0)) {
            // Remalhagem Global
            // O mesh_height precisa ser a altura ATUAL da peça para a remalhagem.
            // Para simplificar, vamos usar a altura inicial para remalhar numa grade de referência,
            // mas o ideal é usar a altura MÁXIMA alcançada ou a altura atual se ela for significativamente menor.
            // Usaremos a altura inicial (no referencial) e o `width` inicial.
            remesh_global_mesh(nodes, elements, width, height_initial, remesh_new_divisions_x, remesh_new_divisions_y, current_material_properties);

            // Após remalhagem, os deslocamentos acumulados são transferidos, mas as forças de reação e estados de contato
            // precisam ser reavaliados no próximo passo do solver.
            // É crucial que os IDs dos nós sejam remapeados corretamente para as BCs.
            // O fixed_ux_node_id_main também precisa ser re-localizado
            fixed_ux_node_id_main = -1;
            min_dist_x_main = std::numeric_limits<double>::max();
            for (const auto& node : nodes) {
                if (std::abs(node.y.value) < SMALL_VALUE_DOUBLE) {
                    double current_dist_x = std::abs(node.x.value - center_x_main.value);
                    if (current_dist_x < min_dist_x_main) {
                        min_dist_x_main = current_dist_x;
                        fixed_ux_node_id_main = node.id;
                    }
                }
            }
            // Zera as reações e estados de contato após remalhar, eles serão recalculados no próximo NR
            for(auto& node : nodes) {
                node.fx_reaction = N(0.0);
                node.fy_reaction = N(0.0);
                node.contact_state = 0;
            }
        }


        // Cálculo para a curva Tensão-Deformação
        Strain current_strain = dimensionless(0.0);
        Pressure current_stress = Pa(0.0);

        Length sum_uy_top = m(0.0);
        int count_top_nodes = 0;
        for (const auto& node : nodes) {
            if (std::abs(node.y.value - height_initial.value) < SMALL_VALUE_DOUBLE) {
                sum_uy_top = sum_uy_top + node.uy;
                count_top_nodes++;
            }
        }
        Length avg_uy_top = (count_top_nodes > 0) ? sum_uy_top / static_cast<double>(count_top_nodes) : m(0.0);
        current_strain = dimensionless(std::abs(avg_uy_top.value) / height_initial.value);

        Force total_normal_force_top = N(0.0);
        if (choice == 1) {
            total_normal_force_top = current_force_at_step_target;
        } else {
            for (const auto& node : nodes) {
                if (std::abs(node.y.value - initial_height.value) < SMALL_VALUE_DOUBLE) {
                    total_normal_force_top = total_normal_force_top + N(std::abs(node.fy_reaction.value));
                }
            }
        }
        current_stress = total_normal_force_top / (width * m(thickness));

        stress_strain_curve.push_back({current_strain, current_stress});

        current_load_factor = next_load_factor;

        if (step_number >= 500) {
            std::cerr << "AVISO: Numero maximo de passos (500) atingido no passo adaptativo. Simulacao encerrada.\n";
            break;
        }
    }

    std::cout << "\nSimulacao concluida. " << step_number + 1 << " arquivos VTU gerados para animacao no ParaView.\n";
    export_stress_strain_data("tensao_deformacao.txt", stress_strain_curve);

    return 0;
}
