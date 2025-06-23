#include <iostream>  // Para entrada e saída no console (cout, cerr)
#include <cmath>     // Para funções matemáticas como log e pow (potência)
#include <fstream>   // Para manipulação de arquivos (ofstream)
#include <iomanip>   // Para formatar a saída do arquivo (setprecision, fixed)
#include <vector>    // Para usar o container std::vector

// --- Classe Material ---
// Representa as propriedades do material sendo forjado.
// Atualmente implementa a Lei de Potência de Hollomon para endurecimento por deformação.
class Material {
public:
    double K;  // Coeficiente de resistência (Strength Coefficient) em Pascals (Pa)
    double n;  // Expoente de encruamento (Strain Hardening Exponent) - adimensional

    // Tensão de escoamento inicial: um valor de referência para pequenas deformações.
    // Calculado como a tensão de fluxo para uma pequena deformação (ex: 0.2%).
    double tensaoEscoamentoInicial;

    // Construtor: inicializa o material com os parâmetros K e n.
    Material(double K_val, double n_val)
        : K(K_val), n(n_val)
    {
        // Define a tensão de escoamento inicial como a tensão de fluxo para 0.2% de deformação
        // Isso ajuda o gráfico Tensão x Deformação a começar de um ponto mais realista.
        tensaoEscoamentoInicial = K * std::pow(0.002, n);
    }

    // Retorna a tensão de fluxo (true stress) para uma dada deformação real (true strain).
    // Usa a Lei de Potência de Hollomon: sigma = K * epsilon^n.
    // Para deformação zero, retorna a tensão de escoamento inicial.
    double getTensaoDeFluxo(double deformacaoReal) const {
        if (deformacaoReal <= 0) {
            return tensaoEscoamentoInicial;
        }
        return K * std::pow(deformacaoReal, n);
    }
};

// --- Classe PecaCilindrica ---
// Representa o tarugo cilíndrico que será forjado.
// Gerencia suas dimensões (altura e raio) e calcula propriedades como área e deformação.
class PecaCilindrica {
private:
    double alturaInicial;  // Altura original da peça em metros (m)
    double raioInicial;    // Raio original da peça em metros (m)
    double volumeInicial;  // Volume inicial da peça (constante durante a deformação plástica)

public:
    double alturaAtual;    // Altura atual da peça em metros (m)
    double raioAtual;      // Raio atual da peça em metros (m)

    // Construtor: inicializa a peça com suas dimensões originais.
    PecaCilindrica(double h0, double r0) :
        alturaInicial(h0),
        raioInicial(r0),
        alturaAtual(h0),
        raioAtual(r0)
    {
        // Calcula o volume inicial da peça (assumido constante na plasticidade)
        volumeInicial = M_PI * raioInicial * raioInicial * alturaInicial;
    }

    // Métodos "getter" para acessar as propriedades privadas (encapsulamento).
    double getAlturaInicial() const {
        return alturaInicial;
    }

    double getRaioInicial() const {
        return raioInicial;
    }

    // Atualiza as dimensões da peça com base em uma nova altura.
    // O raio é recalculado para manter o volume constante.
    void atualizarDimensoes(double novaAltura) {
        if (novaAltura <= 0) {
            std::cerr << "Erro: A altura não pode ser zero ou negativa." << std::endl;
            return;
        }

        alturaAtual = novaAltura;
        // Recalcula o raio para manter o volume constante: V = pi * R^2 * H
        // Novo R = sqrt(V / (pi * H))
        raioAtual = std::sqrt(volumeInicial / (M_PI * alturaAtual));
    }

    // --- CORREÇÃO AQUI: Método getAreaAtual() re-adicionado ---
    // Calcula a área da seção transversal atual da peça.
    double getAreaAtual() const {
        return M_PI * raioAtual * raioAtual;
    }

    // Calcula a deformação real (true strain) na direção da compressão.
    // epsilon = ln(H_inicial / H_atual)
    double getDeformacaoReal() const {
        if (alturaAtual <= 0) { // Evita log de zero ou negativo
            return 0.0;
        }
        return std::log(alturaInicial / alturaAtual);
    }
};

// --- Classe SimuladorForjamento ---
// Orquestra a simulação do processo de forjamento.
// Contém o material e a peça, e os métodos para executar e registrar a simulação.
class SimuladorForjamento {
private:
    Material material;    // Objeto Material usado na simulação
    PecaCilindrica peca;  // Objeto PecaCilindrica sendo forjada

    // Vetores para armazenar os resultados em cada passo da simulação,
    // para posterior exportação e visualização.
    std::vector<double> alturas;
    std::vector<double> raios;
    std::vector<double> deformacoes;
    std::vector<double> forcas;
    std::vector<double> tensoesDeFluxo;

public:
    // Construtor: inicializa o simulador com um material e uma peça.
    SimuladorForjamento(const Material& mat, const PecaCilindrica& p) :
        material(mat), peca(p) {}

    // Simula um único passo de forjamento para uma altura alvo específica.
    // Retorna a força de forjamento necessária para atingir essa altura.
    double simularPasso(double alturaAlvo) {
        peca.atualizarDimensoes(alturaAlvo); // Atualiza as dimensões da peça

        double deformacaoReal = peca.getDeformacaoReal();
        double tensaoDeFluxo = material.getTensaoDeFluxo(deformacaoReal);
        double areaAtual = peca.getAreaAtual();

        double forca = tensaoDeFluxo * areaAtual; // Força = Tensão * Área (modelo sem atrito)
        return forca;
    }

    // Executa a simulação completa do processo de forjamento,
    // do estado inicial até uma altura final em um número de passos.
    void simularProcesso(double alturaFinal, int numPassos) {
        // Validações básicas de entrada
        if (alturaFinal >= peca.getAlturaInicial()) {
            std::cerr << "Erro: A altura final deve ser menor que a altura inicial para forjamento." << std::endl;
            return;
        }
        if (numPassos <= 0) {
            std::cerr << "Erro: O número de passos deve ser positivo." << std::endl;
            return;
        }

        // Limpa os vetores de resultados para uma nova simulação
        alturas.clear();
        raios.clear();
        deformacoes.clear();
        forcas.clear();
        tensoesDeFluxo.clear();

        // Reinicia a peça para suas dimensões originais antes da simulação
        peca.alturaAtual = peca.getAlturaInicial();
        peca.raioAtual = peca.getRaioInicial();

        // Adiciona o estado inicial (Passo 0) aos vetores de resultados
        alturas.push_back(peca.getAlturaInicial());
        raios.push_back(peca.getRaioInicial());
        deformacoes.push_back(0.0); // Deformação real inicial é zero
        forcas.push_back(0.0);      // Força no início da deformação é zero
        tensoesDeFluxo.push_back(material.getTensaoDeFluxo(0.0)); // Tensão de fluxo inicial (tensão de escoamento)

        // Calcula o incremento de altura para cada passo
        double incrementoAltura = (peca.getAlturaInicial() - alturaFinal) / numPassos;

        std::cout << "--- Simulacao de Forjamento ---" << std::endl;
        std::cout << "Altura Inicial: " << peca.getAlturaInicial() * 1000 << " mm, Raio Inicial: " << peca.getRaioInicial() * 1000 << " mm" << std::endl;
        // Exibe os parâmetros do material de Hollomon
        std::cout << "Material (Hollomon): K = " << material.K / 1e6 << " MPa, n = " << material.n << std::endl;

        // Imprime o estado inicial no console
        std::cout << "Passo 0 (Estado Inicial): "
                  << "Altura = " << peca.getAlturaInicial() * 1000 << " mm, "
                  << "Raio = " << peca.getRaioInicial() * 1000 << " mm, "
                  << "Deformacao Real = " << deformacoes.back() << ", "
                  << "Forca = " << forcas.back() / 1000 << " kN" << ", "
                  << "Tensao de Fluxo = " << tensoesDeFluxo.back() / 1e6 << " MPa" << std::endl;

        // Loop principal da simulação, iterando pelos passos de deformação
        for (int i = 1; i <= numPassos; ++i) {
            double alturaParaSimular = peca.getAlturaInicial() - (i * incrementoAltura);
            // Garante que a altura não passe da altura final, especialmente no último passo
            if (alturaParaSimular < alturaFinal) alturaParaSimular = alturaFinal;

            // Simula o passo, o que atualiza peca.alturaAtual e peca.raioAtual
            double forcaNecessaria = simularPasso(alturaParaSimular);

            // Obtém os valores de deformação e tensão de fluxo após o passo
            double deformacaoRealAtual = peca.getDeformacaoReal();
            double tensaoDeFluxoAtual = material.getTensaoDeFluxo(deformacaoRealAtual);

            // Armazena os resultados atuais nos vetores para exportação
            alturas.push_back(peca.alturaAtual);
            raios.push_back(peca.raioAtual);
            deformacoes.push_back(deformacaoRealAtual);
            forcas.push_back(forcaNecessaria);
            tensoesDeFluxo.push_back(tensaoDeFluxoAtual);

            // Imprime os resultados do passo no console
            std::cout << "Passo " << i << ": "
                      << "Altura = " << peca.alturaAtual * 1000 << " mm, "
                      << "Raio = " << peca.raioAtual * 1000 << " mm, "
                      << "Deformacao Real = " << deformacaoRealAtual << ", "
                      << "Forca = " << forcaNecessaria / 1000 << " kN" << ", "
                      << "Tensao de Fluxo = " << tensaoDeFluxoAtual / 1e6 << " MPa" << std::endl;

            // Para o loop se a altura final foi atingida ou ultrapassada
            if (alturaParaSimular <= alturaFinal) break;
        }
        std::cout << "--- Fim da Simulacao ---" << std::endl;
    }

    // Exporta todos os resultados armazenados para um arquivo CSV.
    // O CSV pode ser lido pelo ParaView para gerar gráficos.
    void exportarResultadosCSV(const std::string& nomeArquivo) const {
        std::ofstream arquivoSaida(nomeArquivo);

        if (!arquivoSaida.is_open()) {
            std::cerr << "Erro ao abrir o arquivo: " << nomeArquivo << std::endl;
            return;
        }

        // Configura a precisão de saída para números de ponto flutuante no arquivo
        arquivoSaida << std::fixed << std::setprecision(6);

        // --- CABEÇALHO DA TENSÃO DE FLUXO EM 'MPa' ---
        arquivoSaida << "Altura (m),Raio (m),Deformacao Real,Forca (N),Tensao de Fluxo (MPa)\n";

        // Itera sobre os vetores de resultados e escreve cada linha no formato CSV
        for (size_t i = 0; i < alturas.size(); ++i) {
            arquivoSaida << alturas[i] << ","
                         << raios[i] << ","
                         << deformacoes[i] << ","
                         << forcas[i] << ","
                         << tensoesDeFluxo[i] / 1e6 << "\n"; // --- DIVIDINDO POR 1e6 PARA CONVERTER PARA MPa ---
        }

        arquivoSaida.close();
        std::cout << "Resultados exportados para " << nomeArquivo << std::endl;
    }
};

// --- Função Principal (main) para Teste ---
// Ponto de entrada do programa.
int main() {
    // Define as unidades base para os cálculos: metros (m) para comprimento,
    // Newtons (N) para força, e Pascals (Pa) para tensão.

    // Define o material como Aço SAE 1045, usando a Lei de Hollomon.
    // K = 700 MPa (700 * 10^6 Pa) é o coeficiente de resistência.
    // n = 0.15 é o expoente de encruamento.
    Material aco1045(700e6, 0.15);

    // Define a peça cilíndrica inicial (tarugo).
    // Altura inicial de 100 mm (0.1 m).
    // Raio inicial de 50 mm (0.05 m).
    PecaCilindrica tarugo(0.1, 0.05);

    // Cria uma instância do simulador de forjamento.
    SimuladorForjamento simulador(aco1045, tarugo);

    // Executa a simulação do forjamento.
    // A peça será forjada de 100 mm até uma altura final de 50 mm (0.05 m),
    // em 20 passos incrementais para uma curva mais suave.
    simulador.simularProcesso(0.05, 20);

    // Exporta os resultados da simulação para um arquivo CSV.
    // Este arquivo pode ser aberto no ParaView para gerar gráficos (Tensão x Deformação, Força x Altura, etc.).
    simulador.exportarResultadosCSV("resultados_forjamento.csv");

    system("PAUSE");
    return 0; // Indica que o programa terminou com sucesso
}
