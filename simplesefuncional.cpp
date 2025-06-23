#include <iostream>  // Para entrada e sa�da no console (cout, cerr)
#include <cmath>     // Para fun��es matem�ticas como log e pow (pot�ncia)
#include <fstream>   // Para manipula��o de arquivos (ofstream)
#include <iomanip>   // Para formatar a sa�da do arquivo (setprecision, fixed)
#include <vector>    // Para usar o container std::vector

// --- Classe Material ---
// Representa as propriedades do material sendo forjado.
// Atualmente implementa a Lei de Pot�ncia de Hollomon para endurecimento por deforma��o.
class Material {
public:
    double K;  // Coeficiente de resist�ncia (Strength Coefficient) em Pascals (Pa)
    double n;  // Expoente de encruamento (Strain Hardening Exponent) - adimensional

    // Tens�o de escoamento inicial: um valor de refer�ncia para pequenas deforma��es.
    // Calculado como a tens�o de fluxo para uma pequena deforma��o (ex: 0.2%).
    double tensaoEscoamentoInicial;

    // Construtor: inicializa o material com os par�metros K e n.
    Material(double K_val, double n_val)
        : K(K_val), n(n_val)
    {
        // Define a tens�o de escoamento inicial como a tens�o de fluxo para 0.2% de deforma��o
        // Isso ajuda o gr�fico Tens�o x Deforma��o a come�ar de um ponto mais realista.
        tensaoEscoamentoInicial = K * std::pow(0.002, n);
    }

    // Retorna a tens�o de fluxo (true stress) para uma dada deforma��o real (true strain).
    // Usa a Lei de Pot�ncia de Hollomon: sigma = K * epsilon^n.
    // Para deforma��o zero, retorna a tens�o de escoamento inicial.
    double getTensaoDeFluxo(double deformacaoReal) const {
        if (deformacaoReal <= 0) {
            return tensaoEscoamentoInicial;
        }
        return K * std::pow(deformacaoReal, n);
    }
};

// --- Classe PecaCilindrica ---
// Representa o tarugo cil�ndrico que ser� forjado.
// Gerencia suas dimens�es (altura e raio) e calcula propriedades como �rea e deforma��o.
class PecaCilindrica {
private:
    double alturaInicial;  // Altura original da pe�a em metros (m)
    double raioInicial;    // Raio original da pe�a em metros (m)
    double volumeInicial;  // Volume inicial da pe�a (constante durante a deforma��o pl�stica)

public:
    double alturaAtual;    // Altura atual da pe�a em metros (m)
    double raioAtual;      // Raio atual da pe�a em metros (m)

    // Construtor: inicializa a pe�a com suas dimens�es originais.
    PecaCilindrica(double h0, double r0) :
        alturaInicial(h0),
        raioInicial(r0),
        alturaAtual(h0),
        raioAtual(r0)
    {
        // Calcula o volume inicial da pe�a (assumido constante na plasticidade)
        volumeInicial = M_PI * raioInicial * raioInicial * alturaInicial;
    }

    // M�todos "getter" para acessar as propriedades privadas (encapsulamento).
    double getAlturaInicial() const {
        return alturaInicial;
    }

    double getRaioInicial() const {
        return raioInicial;
    }

    // Atualiza as dimens�es da pe�a com base em uma nova altura.
    // O raio � recalculado para manter o volume constante.
    void atualizarDimensoes(double novaAltura) {
        if (novaAltura <= 0) {
            std::cerr << "Erro: A altura n�o pode ser zero ou negativa." << std::endl;
            return;
        }

        alturaAtual = novaAltura;
        // Recalcula o raio para manter o volume constante: V = pi * R^2 * H
        // Novo R = sqrt(V / (pi * H))
        raioAtual = std::sqrt(volumeInicial / (M_PI * alturaAtual));
    }

    // --- CORRE��O AQUI: M�todo getAreaAtual() re-adicionado ---
    // Calcula a �rea da se��o transversal atual da pe�a.
    double getAreaAtual() const {
        return M_PI * raioAtual * raioAtual;
    }

    // Calcula a deforma��o real (true strain) na dire��o da compress�o.
    // epsilon = ln(H_inicial / H_atual)
    double getDeformacaoReal() const {
        if (alturaAtual <= 0) { // Evita log de zero ou negativo
            return 0.0;
        }
        return std::log(alturaInicial / alturaAtual);
    }
};

// --- Classe SimuladorForjamento ---
// Orquestra a simula��o do processo de forjamento.
// Cont�m o material e a pe�a, e os m�todos para executar e registrar a simula��o.
class SimuladorForjamento {
private:
    Material material;    // Objeto Material usado na simula��o
    PecaCilindrica peca;  // Objeto PecaCilindrica sendo forjada

    // Vetores para armazenar os resultados em cada passo da simula��o,
    // para posterior exporta��o e visualiza��o.
    std::vector<double> alturas;
    std::vector<double> raios;
    std::vector<double> deformacoes;
    std::vector<double> forcas;
    std::vector<double> tensoesDeFluxo;

public:
    // Construtor: inicializa o simulador com um material e uma pe�a.
    SimuladorForjamento(const Material& mat, const PecaCilindrica& p) :
        material(mat), peca(p) {}

    // Simula um �nico passo de forjamento para uma altura alvo espec�fica.
    // Retorna a for�a de forjamento necess�ria para atingir essa altura.
    double simularPasso(double alturaAlvo) {
        peca.atualizarDimensoes(alturaAlvo); // Atualiza as dimens�es da pe�a

        double deformacaoReal = peca.getDeformacaoReal();
        double tensaoDeFluxo = material.getTensaoDeFluxo(deformacaoReal);
        double areaAtual = peca.getAreaAtual();

        double forca = tensaoDeFluxo * areaAtual; // For�a = Tens�o * �rea (modelo sem atrito)
        return forca;
    }

    // Executa a simula��o completa do processo de forjamento,
    // do estado inicial at� uma altura final em um n�mero de passos.
    void simularProcesso(double alturaFinal, int numPassos) {
        // Valida��es b�sicas de entrada
        if (alturaFinal >= peca.getAlturaInicial()) {
            std::cerr << "Erro: A altura final deve ser menor que a altura inicial para forjamento." << std::endl;
            return;
        }
        if (numPassos <= 0) {
            std::cerr << "Erro: O n�mero de passos deve ser positivo." << std::endl;
            return;
        }

        // Limpa os vetores de resultados para uma nova simula��o
        alturas.clear();
        raios.clear();
        deformacoes.clear();
        forcas.clear();
        tensoesDeFluxo.clear();

        // Reinicia a pe�a para suas dimens�es originais antes da simula��o
        peca.alturaAtual = peca.getAlturaInicial();
        peca.raioAtual = peca.getRaioInicial();

        // Adiciona o estado inicial (Passo 0) aos vetores de resultados
        alturas.push_back(peca.getAlturaInicial());
        raios.push_back(peca.getRaioInicial());
        deformacoes.push_back(0.0); // Deforma��o real inicial � zero
        forcas.push_back(0.0);      // For�a no in�cio da deforma��o � zero
        tensoesDeFluxo.push_back(material.getTensaoDeFluxo(0.0)); // Tens�o de fluxo inicial (tens�o de escoamento)

        // Calcula o incremento de altura para cada passo
        double incrementoAltura = (peca.getAlturaInicial() - alturaFinal) / numPassos;

        std::cout << "--- Simulacao de Forjamento ---" << std::endl;
        std::cout << "Altura Inicial: " << peca.getAlturaInicial() * 1000 << " mm, Raio Inicial: " << peca.getRaioInicial() * 1000 << " mm" << std::endl;
        // Exibe os par�metros do material de Hollomon
        std::cout << "Material (Hollomon): K = " << material.K / 1e6 << " MPa, n = " << material.n << std::endl;

        // Imprime o estado inicial no console
        std::cout << "Passo 0 (Estado Inicial): "
                  << "Altura = " << peca.getAlturaInicial() * 1000 << " mm, "
                  << "Raio = " << peca.getRaioInicial() * 1000 << " mm, "
                  << "Deformacao Real = " << deformacoes.back() << ", "
                  << "Forca = " << forcas.back() / 1000 << " kN" << ", "
                  << "Tensao de Fluxo = " << tensoesDeFluxo.back() / 1e6 << " MPa" << std::endl;

        // Loop principal da simula��o, iterando pelos passos de deforma��o
        for (int i = 1; i <= numPassos; ++i) {
            double alturaParaSimular = peca.getAlturaInicial() - (i * incrementoAltura);
            // Garante que a altura n�o passe da altura final, especialmente no �ltimo passo
            if (alturaParaSimular < alturaFinal) alturaParaSimular = alturaFinal;

            // Simula o passo, o que atualiza peca.alturaAtual e peca.raioAtual
            double forcaNecessaria = simularPasso(alturaParaSimular);

            // Obt�m os valores de deforma��o e tens�o de fluxo ap�s o passo
            double deformacaoRealAtual = peca.getDeformacaoReal();
            double tensaoDeFluxoAtual = material.getTensaoDeFluxo(deformacaoRealAtual);

            // Armazena os resultados atuais nos vetores para exporta��o
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
    // O CSV pode ser lido pelo ParaView para gerar gr�ficos.
    void exportarResultadosCSV(const std::string& nomeArquivo) const {
        std::ofstream arquivoSaida(nomeArquivo);

        if (!arquivoSaida.is_open()) {
            std::cerr << "Erro ao abrir o arquivo: " << nomeArquivo << std::endl;
            return;
        }

        // Configura a precis�o de sa�da para n�meros de ponto flutuante no arquivo
        arquivoSaida << std::fixed << std::setprecision(6);

        // --- CABE�ALHO DA TENS�O DE FLUXO EM 'MPa' ---
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

// --- Fun��o Principal (main) para Teste ---
// Ponto de entrada do programa.
int main() {
    // Define as unidades base para os c�lculos: metros (m) para comprimento,
    // Newtons (N) para for�a, e Pascals (Pa) para tens�o.

    // Define o material como A�o SAE 1045, usando a Lei de Hollomon.
    // K = 700 MPa (700 * 10^6 Pa) � o coeficiente de resist�ncia.
    // n = 0.15 � o expoente de encruamento.
    Material aco1045(700e6, 0.15);

    // Define a pe�a cil�ndrica inicial (tarugo).
    // Altura inicial de 100 mm (0.1 m).
    // Raio inicial de 50 mm (0.05 m).
    PecaCilindrica tarugo(0.1, 0.05);

    // Cria uma inst�ncia do simulador de forjamento.
    SimuladorForjamento simulador(aco1045, tarugo);

    // Executa a simula��o do forjamento.
    // A pe�a ser� forjada de 100 mm at� uma altura final de 50 mm (0.05 m),
    // em 20 passos incrementais para uma curva mais suave.
    simulador.simularProcesso(0.05, 20);

    // Exporta os resultados da simula��o para um arquivo CSV.
    // Este arquivo pode ser aberto no ParaView para gerar gr�ficos (Tens�o x Deforma��o, For�a x Altura, etc.).
    simulador.exportarResultadosCSV("resultados_forjamento.csv");

    system("PAUSE");
    return 0; // Indica que o programa terminou com sucesso
}
