#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
    // variaveis
    // variais do primeiro bloco de codigo
    int matriz_particula_estatica[3][3];
    int i,j;
    // variaveis do segundo bloco de codigo
    
    // inicio do primeiro bloco de código
    for(i=0;i<3;i++)
    {
         for(j=0;j<3;j++)
         {
              matriz_particula_estatica[i][j]=0;
         }
    }
     for(i=0;i<3;i++)
     {
          printf("\n\n");
          for(j=0;j<3;j++)
          {
               printf(" %d \t" , matriz_particula_estatica[i][j]);
          }
     }
    // fim do primeiro bloco de codigo
    
    // inicio do segundo bloco de codigo
    
    printf("\n\n\n");
    system("PAUSE");          
    return 0;
}