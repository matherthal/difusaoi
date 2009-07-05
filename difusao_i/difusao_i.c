/*
 * difusao_i.c
 *
 *  Created on: Jun 28, 2009
 *      Author: matheus
 */

#include <stdio.h>
#include <math.h>
#include "decomp_lu.h";

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

int main(void)
{
	FILE *outf;

	void Sair()
	{
		fclose(outf);
		printf("\n Sair da aplicacao");
		return 0;
	}

	float t, t0, tfim, x, x0, xfim, h, k, alfa, lambda;

	int i, j, nx, nt;

	if ((outf = fopen("difusao_i.dat", "w")) == NULL)
	{
		printf("\nProblemas na abertura do arquivo\n");
	}

	h = 0.001; 	/* Discretizacao do tempo */

	k = 1.0/6.0;  	/* Discretizacao do espaco */

	t0   = 0.0;
	tfim = 0.004;

	x0   = 0.0;
	xfim = 1.0;

	nt = (int) ((tfim - t0)/h);
	nx = (int) ((xfim - x0)/k);

	float A[nx-1][nx-1], // A guarda a matriz tridiagonal dos coeficientes     A[linha][coluna]
		b[nx-1], //a matriz igualdade guarda os phi(i,j)    (matriz b)
		X[nx-1], //b é a matriz que guarda os phi(i+1,j)     (matriz X da equação: A*X = b)
		solucao[nx-1][nt-1];	//guarda os vetores b	solucao[linha][coluna]

	fprintf(stdout, "\n Numero de intervalos temporais %d e espaciais %d", nt, nx);

	/*Zerar a matriz solucao*/
	for(i=0; i<nx; i++)
		for(j=0; j<nt; j++)
			solucao[i][j] = 0.0;

	/* Parametros fisicos e variavel auxiliar */

	alfa = 1.0;

	lambda = alfa * h/(k * k);
	fprintf(stdout, "\n k² = %f",lambda);

	x = x0;
	/* Condicao inicial*/
	printf("\n");
	for (j=0;j<nx;j++)
	{
		b[j] = f(x);
		printf("\n b[%d] = %f",j,b[j]);

		x += k;
	}
	b[nx-1] += lambda*k;

	/* Carregar matriz diagonal da memória */
	for (i= 0; i<nx; i++)
	{
		/*fprintf(stdout, "\n i = %d",i);*/
		/*t += h;*/

		for (j=0; j<nx; j++)
		{
			if (j==i)
				A[i][j] = (1.0 + 2.0 * lambda);
			else if ((j==(i+1)) || (j==(i-1)))
				A[i][j] = -lambda;
			else
				A[i][j] = 0.0;
		}

		/*A[nx - 1][nx - 1] -= lambda;*/
	}
	A[0][0] -= lambda;

	printf("\n\n A:");
	for (i= 0; i<nx; i++)		/*IMPRESSAO DA MATRIZ A*/
	{
		printf("\n");

		for (j=0; j<nx; j++)
		{
			printf(" %f ",A[i][j]);
		}
	}

	/*Calcular os vetores de phi(i+1,j) e guardá-los na matriz solucao*/
	for (i=0; i<nt; i++)
	{
		printf("\n i = %d",i);
		if (0 == decomp_lu(nx,A,b,X))	/*Calculo X*/
			Sair();
		for (j=0; j<nx; j++)
		{
			b[j] = X[j];	/*Atualizo o valor de b para ser encontrado o novo X*/
			solucao[j][i] = X[j];	/*Guardo em solucao os valores encontrados para X*/
			printf("\nsolucao[%d][%d] = %f",j,i,solucao[j][i]);
			fprintf(outf,"\nsolucao[%d][%d] = %f",j,i,solucao[j][i]);
		}
	}
	/**********************/
	printf("\n\nsolucao[1][0] = %f\n",solucao[1][0]);

	printf("\n\n Solucoes encontradas: ");	/*Impressao em arquivo e na tela das solucoes encontradas*/
	for (i=0; i<nx; i++)
	{
		printf("\n");
		for (j=0; j<nt; j++)
		{
			fprintf(outf, "%f ", solucao[i][j]);
			printf(" i=%d,j=%d,%f ",i,j, solucao[i][j]);
		}
		fprintf(outf,"\n");
	}

	/*x = x0;
	for (j = 0; j < nx; j++)
	{
		fprintf(matrizphi, "%f ", x);
		fprintf(stdout, "%f ", x);
		for (i = 0; i < nt; i++)
		{
			fprintf(matrizphi, "%f ", A[i][j]);
			fprintf(stdout, "%f ", A[i][j]);
		}
		fprintf(matrizphi, "\n");
		fprintf(stdout, "\n");

		x += k;
	}*/
	Sair();
}
