#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>    // para sqrt
using namespace std;

class particula {
public:
    double Forca_Positiva_X = 0;
    double Forca_Negativa_X = 0;

    double Forca_Positiva_Y = 0;
    double Forca_Negativa_Y = 0;

    // Somatório forças X
    double Somatorio_Forca_X() const {
        return Forca_Positiva_X + Forca_Negativa_X;
    }

    // Somatório forças Y
    double Somatorio_Forca_Y() const {
        return Forca_Positiva_Y + Forca_Negativa_Y;
    }

    // Força resultante (módulo do vetor força)
    double Forca_Resultante() const {
        double Fx = Somatorio_Forca_X();
        double Fy = Somatorio_Forca_Y();
        return sqrt(Fx*Fx + Fy*Fy);
    }
};

int main() {
    int n;
    double P;

    cout << "Tamanho da malha: ";
    cin >> n;

    cout << "Carga total (P): ";
    cin >> P;

    // Criando a malha de partículas n x n
    vector<vector<particula>> a(n, vector<particula>(n));

    // Dividir carga igualmente na direção Y e X (vamos supor só Y para começar)
    double b = P / n;

    // Aplicar forças na borda superior (i==0)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                a[i][j].Forca_Negativa_Y = -b;  // força para baixo
                // Exemplo: nenhuma força em X por enquanto
                a[i][j].Forca_Positiva_X = 0;
                a[i][j].Forca_Negativa_X = 0;
                a[i][j].Forca_Positiva_Y = 0;
            }
        }
    }

    // Exibir a matriz de forças
    cout << "\nMatriz das forças em Y (somatório):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j].Somatorio_Forca_Y() << "\t";
        }
        cout << "\n";
    }

    cout << "\nMatriz das forças em X (somatório):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j].Somatorio_Forca_X() << "\t";
        }
        cout << "\n";
    }

    cout << "\nMatriz das forças resultantes:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j].Forca_Resultante() << "\t";
        }
        cout << "\n";
    }

    // Exportar para CSV
    ofstream arquivo("forcas_particulas.csv");
    if (!arquivo) {
        cerr << "Erro ao abrir arquivo para escrita.\n";
        return 1;
    }

    arquivo << "i,j,Forca_Pos_X,Forca_Neg_X,Forca_Pos_Y,Forca_Neg_Y,Somatorio_X,Somatorio_Y,Resultante\n";

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            particula &p = a[i][j];
            arquivo << i << "," << j << ","
                    << p.Forca_Positiva_X << ","
                    << p.Forca_Negativa_X << ","
                    << p.Forca_Positiva_Y << ","
                    << p.Forca_Negativa_Y << ","
                    << p.Somatorio_Forca_X() << ","
                    << p.Somatorio_Forca_Y() << ","
                    << p.Forca_Resultante() << "\n";
        }
    }

    arquivo.close();
    cout << "\nDados exportados para 'forcas_particulas.csv'.\n";

    // --- Preparação para modelagem física futura ---
    // Exemplo de vetor de ligações entre partículas
    // (vamos declarar uma estrutura para as ligações)

    struct Ligacao {
        int i1, j1;   // coordenadas da primeira partícula
        int i2, j2;   // coordenadas da segunda partícula
        double rigidez; // constante da ligação (exemplo)
    };

    // Exemplo: criar uma ligação entre a partícula (0,0) e (0,1)
    vector<Ligacao> ligacoes;
    if (n > 1) {
        ligacoes.push_back({0,0,0,1, 1000.0});
        cout << "Criada ligação entre (0,0) e (0,1) com rigidez 1000.0\n";
    }

    // Aqui futuramente pode implementar cálculo de deformações, reações, etc.

    return 0;
}
