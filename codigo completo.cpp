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
     
     
     // ENTRADAS DO USUARIO
     
     // dimensoes da malha
     cout << " altura da malha: ";
     cin >> h;
     
     cout << " largura da malha: ";
     cin >> w;
     
     cout << " Quantidade de pontos da malha: ";
     cin >> n;
     
     m = pow( n , 0.5 );
     
     n = (int)m;
     
     malha a[n][n];
    
     // forca
     printf("\n\n");
     cout << " Carga: ";
     cin >> P;
    
     N = P;
    
    for(i=0;i<n;i++)
    {
          
          for(j=0;j<n;j++)
          {
          
               a[i][j].Forca_Negativa_Y = ( (-P) / n );
               a[i][j].Forca_Positiva_Y = ( N / n );
               
          }
    }     
    
    //////////////////////////////////////////////////////////
      for(i=0;i<n;i++)
    {
          
          for(j=0;j<n;j++)
          {
          
               if(i == 0)
               {
               
                       a[i][j].Forca_Negativa_Y = ( (-P) / n );
                       a[i][j].Forca_Positiva_Y = 0;
               
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
               
                       a[i][j].Forca_Negativa_Y = 0;
                       a[i][j].Forca_Positiva_Y = ( N / n );
               
               }
          }
    }     
    ///////////////////////////////////////////////////////////////// 
     
     
     for(i=0;i<n;i++)
     {
            
          printf("\n\n");
            
          for(j=0;j<n;j++)
          {
          
               cout << "\t" << a[i][j].Somatorio_de_Forca_Y();
          
          }
       
     }
     
     
     printf("\n\n\n");
     
     system("PAUSE");          
     
     return 0;
}
