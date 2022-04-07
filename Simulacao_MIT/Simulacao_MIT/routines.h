#ifndef routinesH
#define routinesH

typedef struct {           /* Matrix structure for C library  */
    int row;               /* These are meant to be "private" */
    int col;               /* and should only be accessed via */
    double** mat;          /* the "member functions" below.   */
} matrix;

const int rang = 35;

typedef double vecteur[rang];
typedef double matrice[rang][rang];

matrix* mmake(int rows, int cols);
void mfree(matrix* ptm);
void multmat_mat(matrice a, matrice b, matrice s);
void mulmat_const(matrice a, double factor, matrice s);
void multmat_vect(matrice a, vecteur b, vecteur s);
void sousvect_vect(vecteur c, vecteur d, vecteur s);
void addvect_vect(vecteur c, vecteur d, vecteur s);
void inv(matrice L, matrice inv_L);
void multvect_const(vecteur a, double b, vecteur s);


#endif