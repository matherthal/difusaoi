#include <stdio.h>

/*Se retorna 0 houve erro, senao nao*/
int decomp_lu(int n, float a[n-1][n-1], float b[n-1], float x[n-1])
{
	int i,j,k;
	float M, E;

	/*printf("Seja o sistema de equa��es lineares
	 * \na11 x1 + a12 x2 + ... + a1n xn = b1
	 * \na21 x1 + a22 x2 + ... + a2n xn = b2
	 * \n.\t.\t.\t.\t.
	 * \nan1 x1 + an2 x2 + ... + ann xn = b2
	 * \nde n equa��es a n inc�gnitas a ser resolvido. Digite n (n <= 10).");*/

	/*scanf("%d",&n);*/
	float m[n][n], l[n][n];
	float y[n], s[n];

	for(i=0;i<n;i++)		/*Inicializacao das matrizes*/
	{
		for(j=0;j<n;j++)
		{
			/*a[i][j] = 0;*/
			m[i][j] = 0.0;
			l[i][j] = 0.0;
		}
		/*b[i] = 0;*/
		y[i] = 0.0;
		s[i] = 0.0;
		/*x[i] = 0.0;*/
	}

	/*printf("\n\n DEBUG decomp_lh.c: ");	/*Impressao em arquivo e na tela das solucoes encontradas
	for (i=0; i<n; i++)
	{
		printf("\n");
		for (j=0; j<n; j++)
		{
			printf(" %f ", a[i][j]);
		}
	}*/

	/*printf("\nATEN��O: os elementos da diagonal principal devem ser n�o nulos!\n");*/
	for (i = 0; i < n; i++)	/*Verificacao de se ha elementos nulos na diagonal principal de a*/
		if (a[i][i] == 0)
		{
			printf("\n ERRO: Ha elementos nulos na diagonal principal de a:\n a[%d][%d] = %f",i,i,a[i][i]);
			return 0;
		}

	/*for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			LABEL:
			printf("Digite a[%d][%d]: ",i+1,j+1);
			scanf("%f",&a[i][j]);
			if(i==j && a[i][j]==0)
			{
				goto LABEL;
			}
		}
	}*/

	/*for(i=0;i<n;i++)
	{
		printf("Digite b[%d]: ",i+1);
		scanf("%f",&b[i]);
	}*/

	for(j=0;j<n-1;j++)	/*Escalonamento*/
	{
		for(i=j+1;i<n;i++)
		{
			M = a[i][j] / a[j][j];
			a[i][j] = 0;
			m[i][j] = M;
			for(k=j+1;k<n;k++)
			{
				a[i][k] = a[i][k] - (M * a[j][k]);
			}
		}
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
			{
				l[i][j] = 1;
			}
			else
			{
				if(i>j)
				{
					l[i][j] = m[i][j];
				}
			}
		}
	}

	y[0] = b[0];
	for(k=1;k<n;k++)
	{
		y[k] = b[k];

		for(j=k;j>=0;j--)
		{
			s[k] = s[k] + l[k][j-1] * y[j-1];
			y[k] = b[k] - s[k];
		}
	}

	for(i=n-1;i>=0;i--)
	{
		E=0;

		if(i!=(n-1))
		{
			for(k=n-1;k>=i+1;k--)
			{
				E = E + (a[i][k] * x[k]);
			}
		}
		x[i] = (y[i] - E) / a[i][i];
	}

	/*printf("\n");*
	for(i=0;i<n;i++)
	{
		printf("X[%d] = %5.5f\n",i+1,x[i]);
	}*/

	return 1;
}
