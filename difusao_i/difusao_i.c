/*
 * difusao_i.c
 *
 *  Created on: Jun 28, 2009
 *      Author: matheus, giulio
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "decomp_lu.h"

/* Programa de demonstracao que implementa o metodo de resolucao de equacao */
/* parabolica  por um metodo explicito */

float f(float x)
/* Funcao que da' a condicao inicial */
{
	if ((x > 0.0) && (x <= 0.3))
		return(100.0 * x);
	else
		return(0);
}

void Sair()
{
	printf("\n Sair da aplicacao");
	getchar();
	getchar();
}

int main()
{
	float t, t0, tfim, x, x0, xfim, h, k, alfa, lambda;

	int i, j, nx, nt;

	h = 0.001; 	/* Discretizacao do tempo */

	k = 1.0/6.0;  	/* Discretizacao do espaco */

	t0   = 0.0;
	tfim = 0.004;

	x0   = 0.0;
	xfim = 1.0;

	nt = (int) ((tfim - t0)/h) + 1;
	nx = (int) ((xfim - x0)/k) + 2;

	float* A[nx]; // A guarda a matriz tridiagonal dos coeficientes     A[linha][coluna]
	float	b[nx]; //a matriz igualdade guarda os phi(i,j)    (matriz b)
	float	X[nx]; //b é a matriz que guarda os phi(i+1,j)     (matriz X da equação: A*X = b)
	float	solucao[nx][nt];	//guarda os vetores b	solucao[linha][coluna]

	for (i = 0; i < nx; i++)
		A[i] = (float *)malloc(nx*sizeof(float));
	fprintf(stdout, "\n Numero de intervalos temporais %d e espaciais %d", nt, nx);

	/*Zerar a matriz solucao*/
	for(i = 0; i < nx; i++)
		for(j = 0; j < nt; j++)
			solucao[i][j] = 0.0;

	/* Parametros fisicos e variavel auxiliar */

	alfa = 1.0;

	lambda = alfa * h/(k * k);
	printf("\n k^2 = %f",lambda);

	x = x0;
	/* Condicao inicial*/
	printf("\n");
	for (j = 0; j < nx; j++)
	{
		b[j] = f(x);
		x += k;
	}
	/*b[0] += lambda*k;*/
	b[nx-1] += lambda*k;
	for (j = 0; j < nx; j++)
	{
		printf("\n b[%d] = %f",j,b[j]);
	}

	/* Carregar matriz diagonal da memória */
	for (i = 0; i < nx; i++)
	{
		/*fprintf(stdout, "\n i = %d",i);*/
		/*t += h;*/

		for (j = 0; j < nx; j++)
		{
			if (j == i)
				A[i][j] = (1.0 + 2.0 * lambda);
			else if ((j==(i+1)) || (j==(i-1)))
				A[i][j] = -lambda;
			else
				A[i][j] = 0.0;
		}

		/*A[nx - 1][nx - 1] -= lambda;*/
	}
	A[0][0] -= lambda;
	A[nx-1][nx-1] -= lambda;

	getchar();

	printf("\n\n A:");
	for (i = 0; i < nx; i++)		/*IMPRESSAO DA MATRIZ A*/
	{
		printf("\n");

		for (j = 0; j < nx; j++)
		{
			printf(" %f ",A[i][j]);
		}
	}

	/*Calcular os vetores de phi(i+1,j) e guardá-los na matriz solucao*/
	for (i = 0; i < nx; i++)
	{
		float *L[nx];

		printf("\n i = %d",i);
		decomp_lu(L,A,nx);
		for (j = 0; j < nx; j++)
		{
			float sol[nx];
			solve_lu(sol,L,A,b,nx);
			int k;
			for (k = 0; k < nx; k++)
			{
				solucao[k][j] = sol[k];
				b[k] = sol[k];
			}
			b[nx-1] += lambda*k;
			printf("\nsolucao[%d][%d] = %f",j,i,solucao[j][i]);
		}
		for (i = 0; i < nx; i++)
			free(L[i]);
	}
	for (i = 0; i < nx; i++)
		free(A[i]);

	printf("\n\n Solucoes encontradas: ");	/*Impressao em arquivo e na tela das solucoes encontradas*/
	for (i = 0; i < nx; i++)
	{
		printf("\n");
		for (j = 0; j < nt; j++)
			printf(" i=%d,j=%d %f ",i,j, solucao[i][j]);
	}

	FILE *outf = fopen("difusao_i.dat", "w");

	x = x0;
	for (i = 0; i < nx; i++)
	{
		fprintf(outf,"%f",x);
		for (j = 0; j < nt; j++)
		{
			fprintf(outf," %f",solucao[i][j]);
		}
		fprintf(outf,"\n");
		x += k;
	}
	fclose(outf);

	Sair();
	return 0;
}
