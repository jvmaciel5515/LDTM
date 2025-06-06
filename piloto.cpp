#include <stdio.h>
#include <iostream>
using namespace std;


class particula
{
     
     public:
          
          double Forca_Positiva_Y;
          
          double Forca_Negativa_Y;
     
          double Somatorio_Forca_Y();


};

double particula::Somatorio_Forca_Y()
{

     return Forca_Negativa_Y + Forca_Positiva_Y;

}

int main()
{
    
     int i,j;
     
     int n;
     
     double P;
     
     double b;
     
     
     cout << " Tamanho da malha: ";
     cin >> n;
     
     printf("\n\n");
     
     cout << " Carga: ";
     cin >> P;
     
     particula a[n][n];
     
     b = P / n;
    
    for(i=0;i<n;i++)
    {
          
          for(j=0;j<n;j++)
          {
          
               if(i == 0)
               {
                    
                    a[i][j].Forca_Negativa_Y = -b;
                    a[i][j].Forca_Positiva_Y = 0;
               
               }
               
               else
               {
                    
                    a[i][j].Forca_Negativa_Y = 0;
                    a[i][j].Forca_Positiva_Y = 0;
               
               }
          }
    }     
     
     
     for(i=0;i<n;i++)
     {
            
          printf("\n\n");
            
          for(j=0;j<n;j++)
          {
          
               cout << "\t" << a[i][j].Somatorio_Forca_Y();
          
          }
       
     }
     
     
     
     printf("\n\n\n");
     
     system("PAUSE");          
     
     return 0;
}
