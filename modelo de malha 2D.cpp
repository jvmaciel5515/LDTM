#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;

// caracterizacao da malha
class malha
{
     
     public:
          
          // forças em Y
          double Forca_Positiva_Y;
          
          double Forca_Negativa_Y;
          
          double Somatorio_de_Forca_Y();
          
          
          // forças em X
          double Forca_Positiva_X;
          
          double Forca_Negativa_X;
          
          double Somatorio_de_Forca_X();
          
          
          // posicao inicial
          double Posicao_Inicial_Y;
          
          double Posicao_Inicial_X;
          
          
          // posicao final
          double Posicao_Final_Y;
          
          double Posicao_Final_X;
          
          
};

double malha::Somatorio_de_Forca_Y()
{

     return Forca_Positiva_Y + Forca_Negativa_Y;

}

double malha::Somatorio_de_Forca_X()
{

     return Forca_Positiva_X + Forca_Negativa_X;

}


int main()
{
    
     int i;
     
     int j;
     
     int n;
     
     double P;
     
     double N;
     
     double m;
     
     double h;
     
     double w;
     
     double b;
     
     
     // ENTRADAS DO USUARIO
     
     // dimensoes da malha
     cout << " Altura da malha: ";
     cin >> h;
     
     printf("\n");
     
     cout << " Largura da malha: ";
     cin >> w;
     
     printf("\n");
     
     cout << " Quantidade de pontos da malha: ";
     cin >> n;
     
     m = pow( n , 0.5 );
     
     n = (int)m;
     
     b = n - 1;
     
     malha malha1[n][n];
    
     for(i=0;i<n;i++)
     {
          
          for(j=0;j<n;j++)
               
               {
          
                    malha1[i][j].Posicao_Inicial_X = (j * w / b);
                    malha1[i][j].Posicao_Inicial_Y = ((-1) * i * h / b) + h;
               
               }
     }     
     
       for(i=0;i<n;i++)
       {
            
          printf("\n\n");
            
          for(j=0;j<n;j++)
          {
          
               cout << "(" << malha1[i][j].Posicao_Inicial_X;
               cout << "," << malha1[i][j].Posicao_Inicial_Y; 
               cout << ")";
               printf("\t");
          
          }
       
     }
     
     
     
     
     
     
     
     // forca
     printf("\n\n");
     cout << " Carga: ";
     cin >> P;
    
     N = P;
    
    for(i=0;i<n;i++)
    {
          
          for(j=0;j<n;j++)
          {
          
               malha1[i][j].Forca_Negativa_Y = ( (-P) / n );
               malha1[i][j].Forca_Positiva_Y = ( N / n );
               
          }
    }     
    
    //////////////////////////////////////////////////////////
      for(i=0;i<n;i++)
    {
          
          for(j=0;j<n;j++)
          {
          
               if(i == 0)
               {
               
                       malha1[i][j].Forca_Negativa_Y = ( (-P) / n );
                       malha1[i][j].Forca_Positiva_Y = 0;
               
               }
          }
    }     
    ///////////////////////////////////////////////////////////////
      for(i=0;i<n;i++)
    {
          
          for(j=0;j<n;j++)
          {
          
               if(i == (n-1))
               {
               
                       malha1[i][j].Forca_Negativa_Y = 0;
                       malha1[i][j].Forca_Positiva_Y = ( N / n );
               
               }
          }
    }     
    ///////////////////////////////////////////////////////////////// 
     
     
     for(i=0;i<n;i++)
     {
            
          printf("\n\n");
            
          for(j=0;j<n;j++)
          {
          
               cout << "\t" << malha1[i][j].Somatorio_de_Forca_Y();
          
          }
       
     }
     
     
     printf("\n\n\n");
     
     system("PAUSE");          
     
     return 0;
}
