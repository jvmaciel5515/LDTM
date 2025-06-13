#include <stdio.h>
#include <iostream>
using namespace std;

// caracterizacao do objeto2D
class objeto2D
{
     
     public:
            
            double Posicao_X;
            
            double Posicao_Y;
            
};


int main()
{
     
     double h;
     
     double w;
     
     int i;
     
     int j;
     
     
     // ENTRADAS DO USUARIO
     
     // objeto 2D
     cout << " altura do objeto em milimetros: ";
     cin >> h;
     
     printf("\n");
     
     cout << " largura do objeto em milimetros: ";
     cin >> w;
     
     objeto2D objeto1[2][2];
     
     objeto1[0][0].Posicao_X = 0;
     objeto1[0][0].Posicao_Y = h;
     
     objeto1[0][1].Posicao_X = w;
     objeto1[0][1].Posicao_Y = h;
     
     objeto1[1][0].Posicao_X = 0;
     objeto1[1][0].Posicao_Y = 0;
     
     objeto1[1][1].Posicao_X = w;
     objeto1[1][1].Posicao_Y = 0;
     
     
     
     // imprime na tela
     printf("\n");
     
     for(i=0;i<2;i++)
     {
          printf("\n\n");
          for(j=0;j<2;j++)
          {
               
               cout << "(" << objeto1[i][j].Posicao_X;
               cout << "," << objeto1[i][j].Posicao_Y; 
               cout << ")";
               printf("\t");
          
          }
     
     }
     
     printf("\n\n\n");
     
     cout << " objeto 2D gerado com sucesso ";
     
     printf("\n\n\n");
     
     system("PAUSE");          
     
     return 0;
}
