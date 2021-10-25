#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "routines.h"
#include <math.h>


matrix* mmake(int rows, int cols)
{
    matrix* ptm;
    double** row_pointers;
    int i;

    /* Check that no dimension is zero. */
    if ((rows == 0) || (cols == 0)) puts("Invalid dimension error in mmake!");


    /* Memory allocation for matrix structure. */
    ptm = (matrix*)malloc(sizeof(matrix));

    /* Memory for array of pointers to rows. */;
    row_pointers = (double**)malloc(rows * sizeof(double*));

    /* Memory for all rows, initialize row pointers. */
    row_pointers[0] = (double*)malloc(rows * cols * sizeof(double));
    for (i = 1; i < rows; i++) {
        row_pointers[i] = row_pointers[i - 1] + cols;
    }

    /* Check if last allocation was ok ! */
    if (!row_pointers[0]) {
        puts("Memory allocation error in mmmake!");
    }

    ptm->row = rows;             /* Initialize matrix structure */
    ptm->col = cols;
    ptm->mat = row_pointers;     /* Pointer to row pointers     */

    return ptm;           /* Return pointer to matrix structure */
}

void mfree(matrix* ptm)
{
    /* Deallocate rows */
    free(ptm->mat[0]);

    /* Deallocate row pointer array */
    free(ptm->mat);

    /* Deallocate matrix structure */
    free(ptm);
}

void multmat_mat(matrice a, matrice b, matrice s)
{
    int i, j, k;
    for (i = 1; i < rang; i++)
    {
        for (j = 1; j < rang; j++)
        {
            s[i][j] = 0.0;
            for (k = 1; k < rang; k++)
            {
                s[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void mulmat_const(matrice a, double factor, matrice s)
{
    int i, j;
    for (i = 1; i < rang; i++)
    {
        for (j = 1; j < rang; j++)
        {
            s[i][j] = factor * a[i][j];
        }
    }
}

void multmat_vect(matrice a, vecteur b, vecteur s)
{
    int l, m;
    double temp;

    for (l = 1; l < rang; l++)
    {
        temp = 0.0;
        for (m = 1; m < rang; m++)
            temp = temp + a[l][m] * b[m];
        s[l] = temp;
    }
}

void sousvect_vect(vecteur c, vecteur d, vecteur s)
{
    int i;

    for (i = 1; i < rang; i++)
        s[i] = c[i] - d[i];
}

void addvect_vect(vecteur c, vecteur d, vecteur s)
{
    int i;

    for (i = 1; i < rang; i++)
        s[i] = c[i] + d[i];
}

void multvect_const(vecteur a, double b, vecteur s)
{
    int i;

    for (i = 1; i < rang; i++)
        s[i] = a[i] * b;
}

void inv(matrice L, matrice inv_L)
{
    int i, j;
    double r, pivot;
    double det = 0;
    //static double det = 0;
    int ligcour, ligpivot;

    for (i = 1; i < rang; i++)
        for (j = 1; j < rang; j++)
            if (i == j) inv_L[i][j] = 1.0; else inv_L[i][j] = 0.0;

    for (ligcour = 1; ligcour < rang; ligcour++)
    {
        ligpivot = ligcour;
        while ((fabs(L[ligpivot][ligcour]) < 1.0e-12) && (ligpivot < 20))
        {
            ligpivot++;
            det = -det;
        }
        if ((fabs(L[ligpivot][ligcour]) < 1.0e-12) && (ligpivot == 20))
            det = 0;
        else
        {
            for (j = 1; j < rang; j++)
            {
                r = L[ligcour][j];
                L[ligcour][j] = L[ligpivot][j];
                L[ligpivot][j] = r;
                r = inv_L[ligcour][j];
                inv_L[ligcour][j] = inv_L[ligpivot][j];
                inv_L[ligpivot][j] = r;
            }
            pivot = L[ligcour][ligcour];
            for (j = ligcour + 1; j < rang; j++)
                L[ligcour][j] /= pivot;
            for (j = 1; j < rang; j++)
                inv_L[ligcour][j] /= pivot;
            det *= pivot;
            L[ligcour][ligcour] = 1.0;
            for (i = 1; i < rang; i++)
                if (i != ligcour)
                {
                    r = L[i][ligcour];
                    for (j = ligcour + 1; j < rang; j++)
                        L[i][j] -= r * L[ligcour][j];
                    for (j = 1; j < rang; j++)
                        inv_L[i][j] -= r * inv_L[ligcour][j];
                    L[i][ligcour] = 0.0;
                }
        }
    }
}