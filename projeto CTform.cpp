// BIBLIOTECAS
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;



// MATRIZ DE PONTOS
class matriz_de_pontos
{     
     
     public:
            
           double Coordenada_X;     
     
           double Coordenada_Y;
           
};


// MATRIZ DE ELEMENTOS DE AREA
class matriz_de_areas
{

     public:
            
            double Dimensao_X;
            
            double Dimensao_Y;

            double Area();
            
};

double matriz_de_areas::Area()
{

     return Dimensao_X * Dimensao_Y;

}


// MATRIZ DE TENSOES
class matriz_de_tensoes
{


};


// MATRIZ DE DEFORMAÇÕES
class matriz_de_deformacoes
{
  
  
      
};      


// FUNÇÃO PRINCIPAL
int main()
{
     
     // VARIAVEIS
     double h;
     
     double w;
     
     double m;
     
     int i;
     
     int j;
     
     int n;
     
     int b;
     
     
     // MATRIZ DE PONTOS //////////////////////////////////////////////////////////////////
     
     // dimensoes da matriz de pontos
     cout << " Altura da matriz de pontos: ";
     cin >> h;
     
     printf("\n");
     
     cout << " Largura da matriz de pontos: ";
     cin >> w;
     
     printf("\n");
     
     cout << " Quantidade de pontos da matriz de pontos: ";
     cin >> n;
     
     m = pow( n , 0.5 );
     
     n = (int)m;
     
     b = n - 1;
     
     matriz_de_pontos matriz_1[n][n];
    
     for(i=0;i<n;i++)
     {
          
          for(j=0;j<n;j++)
               
               {
          
                    matriz_1[i][j].Coordenada_X = (j * w / b);
                    
                    matriz_1[i][j].Coordenada_Y = ((-1) * i * h / b) + h;
               
               }
     }     
     
     printf("\n\n\nmatriz de pontos:\n\n");
     
       for(i=0;i<n;i++)
       {
            
          printf("\n\n");
            
          for(j=0;j<n;j++)
          {
          
               cout << "(" << matriz_1[i][j].Coordenada_X;
               
               cout << "," << matriz_1[i][j].Coordenada_Y; 
               
               cout << ")";
               
               cout << "\t";
          
          }
       
     }
     
    // MATRIZ DE ELEMENTOS DE AREA //////////////////////////////////////////////////////////////////
  
    printf("\n\n\nmatriz de elementos de area:\n\n");
  
      
  
     printf("\n\n\n");
     
     system("PAUSE");          
     
     return 0;
}
