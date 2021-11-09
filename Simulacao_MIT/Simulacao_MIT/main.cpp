/*  Máquina assíncrona trifásica */
/*  Modelo para componente fundamental (sem harmônicas espaciais)  */
/*  Simulação  de defeito no estator (curto-circuito entre espiras) e assimetrias no rotor*/
#define _CRT_SECURE_NO_WARNINGS
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "routines.h"

#define pi     3.14159265358979323846
#define NiterG 300000    //130000;4000;300000(com conversor) /* Número total de pontos da simulação */
#define Niter  40       //20 ;200;40(com conversor)         /* Número loops antes de pegar um  ponto */
#define dt     0.000001 //0.0000025//0.000005 //0.000001 (com conversor) /* DeltaT de simulação dividido por Niter  0.00001 com Niter = 10 para ter Fa = 10kHz -freq. de amostragem*/
#define Umax   311  //1600 //439.82//1600 //311=220*sqrt(2) 129*sqrt(2)=179,61 /*Tensão máxima nos bornes da carga*/

int i, j;

// Dados de Placa do Motor de Indução Trifásico WEG 60
//  (0.75 kW / 1.0 cv, lab.)
//  Classe N
// wn = 1720 rpm = 180.118 rad/s - escorregamento = 4.4444%
// Ip/In= 7.2
// 220/380 V - 3.02/1.75 A
// cos = 0.82    rend. = 79.5 %
// S = 1150.8 VA  P = 943.6351 W
// Cn = 4.1639 Nm

// Dados obtidos a partir dos ensaios classicos realizados pelo Clayton
// Maq. weg 60 (1 cv, lab.)
const double Rs = 12.5;         /*  resistência do estator    */
const double Lls = 0.0383;      /*  Indutância de dispersão do estator  lsl  */
const double Lms = 0.4506;      /*  Indutância principal (própia) do estator  Lsp  */
const double Lsr = 0.4506;      /*  Indutância mútua estator-rotor  Msr */
const double Rr = 8.9;           /*  resistência do rotor   */
const double Llr = 0.0383;      /*  Indutância de dispersão do rotor  lrl */
const double Lmr = 0.4506;      /*  Indutância principal do rotor  Lrp */
const double Lsc = 1.5*Lms+Lls; /*  Indutância cíclica do estator   */
const double Lrc = 1.5*Lmr+Llr; /*  Indutância cíclica do rotor   */
//const double Lsc = 1.5*Lsp+lsl; /*  Indutância cíclica do estator   */
//const double Lrc = 1.5*Lrp+lrl; /*  Indutância cíclica do rotor   */
//const double p = 1.0;           /*  Número de pares de pólos   */
const double p = 1.0;           /*  Número de pares de pólos   */
const double Jt = 0.023976 * 1.;  /*  Momento de Inércia   */
const double fv = 1. * 0.0014439; /*  coeficiente de atrito dinâmico   */

// **********************************************************************

// Maq. weg 60 (1 cv, lab.)
// Dados de catálogo e ensaios a vazio e curto-circuito realizados pelo Claudio
//const double Rs = 7.65;         /*  resistência do estator    */
//const double lsl = 0.0241;      /*  Indutância de dispersão do estator   */
//const double Lsp = 0.4709;      /*  Indutância principal (própia) do estator    */
//const double Msr = 0.4709;      /*  Indutância mútua estator-rotor   */
//const double Rr = 8.166;        /*  resistência do rotor   */
//const double lrl = 0.0257;      /*  Indutância de dispersão do rotor   */
//const double Lrp = 0.4709;      /*  Indutância principal do rotor   */
//const double Lsc = 1.5*Lsp+lsl; /*  Indutância cíclica do estator   */
//const double Lrc = 1.5*Lrp+lrl; /*  Indutância cíclica do rotor   */
//const double p = 2.0;           /*  Número de pares de pólos   */
////const double Jt = 0.00294;      /*  Momento de Inércia do motor  */
//const double Jt = 0.023976;  /*  Momento de Inércia   motor+carga   */
////const double fv = 0.0008169;    /*  coeficiente de atrito dinâmico do motor  */
//const double fv = 0.0014439; /*  coeficiente de atrito dinâmico  motor+carga  */

// **********************************************************************

double oms, Cr, teta, va1, vb1, vc1, isa1, isb1, isc1, ira, irb, irc, t, om, te, Rbf;

/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii*/

double rqi3 = 0.5774;
//
double vsdef, isdef, vl13, vl15, vl42, vl62, vx, ix, zx;
double vsa_charge1, vsa_charge2, vsa_charge3, vn1, vn1n2;
double vsa_charge4, vsa_charge5, vsa_charge6, vn2;
double ps_alim, ps_mot;

matrice L, dL, A, b, inv_L, B, R;
vecteur x, U, BU, Ad;

/* Método específico de inversão de Matriz SVD */
vecteur w;
matrice u, v, inva;

FILE* fichier;

matrix* results;

void calcul_de_R(void)
{
	double c;
	double sinteta;
	double sinteta_pos;
	double sinteta_neg;

	c = -Lsr * x[7] * p; //TODO: Verificar
	sinteta = sin(p * teta);
	sinteta_pos = sin(p * teta + 2.0 * pi / 3.0);
	sinteta_neg = sin(p * teta - 2.0 * pi / 3.0);

	R[1][1] = R[2][2] = R[3][3] = Rs;
	R[4][4] = R[5][5] = R[6][6] = Rr;
	R[1][2] = R[1][3] = R[1][7] = R[1][8] = 0.0;
	R[2][1] = R[2][3] = R[2][7] = R[2][8] = 0.0;
	R[3][1] = R[3][2] = R[3][7] = R[3][8] = 0.0;
	R[4][5] = R[4][6] = R[4][7] = R[4][8] = 0.0;
	R[5][4] = R[5][6] = R[5][7] = R[5][8] = 0.0;
	R[6][4] = R[6][5] = R[6][7] = R[6][8] = 0.0;
	R[7][1] = R[7][2] = R[7][3] = R[7][8] = 0.0;
	R[8][1] = R[8][2] = R[8][3] = R[8][4] = R[8][5] = R[8][6] = R[8][8] = 0.0;
	//R[1][4] = R[2][5] = R[3][6] = R[4][1] = R[5][2] = R[6][3] = c * sinteta;
	//R[1][5] = R[2][6] = R[3][4] = R[4][3] = R[5][1] = R[6][2] = c * sinteta_neg;
	//R[1][6] = R[2][4] = R[3][5] = R[4][2] = R[5][3] = R[6][1] = c * sinteta_pos;
	//R[7][4] = Lsr * p * (x[1] * sinteta + x[2] * sinteta_pos + x[3] * sinteta_neg);
	//R[7][5] = Lsr * p * (x[1] * sinteta_neg + x[2] * sinteta + x[3] * sinteta_pos);
	//R[7][6] = Lsr * p * (x[1] * sinteta_pos + x[2] * sinteta_neg + x[3] * sinteta);
	R[1][4] = R[2][5] = R[3][6] = R[4][1] = R[5][2] = R[6][3] = c * sinteta;
	R[1][5] = R[2][6] = R[3][4] = R[4][3] = R[5][1] = R[6][2] = c * sinteta_pos;
	R[1][6] = R[2][4] = R[3][5] = R[4][2] = R[5][3] = R[6][1] = c * sinteta_neg;
	R[7][4] = Lsr * p * (x[1] * sinteta + x[2] * sinteta_neg + x[3] * sinteta_pos);
	R[7][5] = Lsr * p * (x[1] * sinteta_pos + x[2] * sinteta + x[3] * sinteta_neg);
	R[7][6] = Lsr * p * (x[1] * sinteta_neg + x[2] * sinteta_pos + x[3] * sinteta);
	R[7][7] = fv;
	R[8][7] = -1.0;
}

void calcul_de_L(void)
{
	double costeta;
	double costeta_neg;
	double costeta_pos;

	costeta = cos(p * teta);
	costeta_neg = cos(p * teta - 2.0 * pi / 3.0);
	costeta_pos = cos(p * teta + 2.0 * pi / 3.0);

	/*L[1][1] = L[2][2] = L[3][3] = Lls + Lms;*/
	L[1][1] = L[2][2] = L[3][3] = Lsc;
	//L[1][2] = L[1][3] = L[2][3] = L[2][1] = L[3][1] = L[3][2] = -1.0 / 2.0 * Lms;
	L[1][2] = L[1][3] = L[2][3] = L[2][1] = L[3][1] = L[3][2] = 0;
	L[1][4] = L[2][5] = L[3][6] = L[4][1] = L[5][2] = L[6][3] = Lsr * costeta;
	L[1][5] = L[2][6] = L[3][4] = L[4][3] = L[5][1] = L[6][2] = Lsr * costeta_pos;
	L[1][6] = L[2][4] = L[3][5] = L[4][2] = L[5][3] = L[6][1] = Lsr * costeta_neg;
	//L[4][4] = L[5][5] = L[6][6] = Llr + Lmr;
	L[4][4] = L[5][5] = L[6][6] = Lrc;
	/*L[4][5] = L[4][6] = L[5][4] = L[5][6] = L[6][4] = L[6][5] = -1.0 / 2.0 * Lmr;*/
	L[4][5] = L[4][6] = L[5][4] = L[5][6] = L[6][4] = L[6][5] = 0;
	L[7][1] = L[7][2] = L[7][3] = L[7][4] = L[7][5] = L[7][6] = L[7][8] = 0.0;
	L[8][1] = L[8][2] = L[8][3] = L[8][4] = L[8][5] = L[8][6] = L[8][7] = 0.0;
	L[1][7] = L[2][7] = L[3][7] = L[4][7] = L[5][7] = L[6][7] = 0.0;
	L[1][8] = L[2][8] = L[3][8] = L[4][8] = L[5][8] = L[6][8] = 0.0;
	L[7][7] = Jt;
	L[8][8] = 1.0;

}

// Geração da Tensão de Alimentação
void calcul_de_U(void)
{
	double v1, v3, v5, v7, vdc, vdc2, fs;
	fs = 60.0;
	oms = 2 * pi * fs;

	va1 = Umax * cos(oms * t);
	vb1 = Umax * cos(oms * t - 2.0 * pi / 3.0);
	vc1 = Umax * cos(oms * t + 2.0 * pi / 3.0);

	vsdef = 0.0;

	// Entrada do vetor de Tensão
	U[1] = va1;
	U[2] = vb1;
	U[3] = vc1;
	U[4] = U[5] = U[6] = 0.0;  // Tensões Rotóricas
	U[7] = -Cr;            // Conjugado de Carga (resistente)
	U[8] = 0.0;
}


/*****************************************************************************/


void calcul_vecteur(double d[rang], double c[rang])
{
	// As funções de inversão, multiplicação e subtração de matrizes estão no programa routines
	calcul_de_L();          /* Calcula-se a matriz L */

	calcul_de_R();          /* Matriz R */

	inv(L, inv_L);  /* Inversão de L */

	// Método do Pivô de Gauss (foi usado aqui pq é mais rápido)

	multmat_mat(inv_L, R, A);
	multmat_vect(A, d, Ad);  // A * x
	multmat_vect(inv_L, U, BU);  // B * U   B = inv_L;
	sousvect_vect(BU, Ad, c);
}

// Inicialização das Variáveis
void init()
{

	int l, c;

	for (l = 0; l < rang; l++) x[l] = 0;
	for (l = 0; l < 8; l++)               // Inicializa-se R com 0
		for (c = 0; c < 8; c++)
			R[l][c] = 0.0;
	for (l = 0; l < rang; l++)             // Inicializa-se b e A com 0
		for (c = 0; c < rang; c++)
		{
			b[l][c] = 0.0;
			A[l][c] = 0.0;
		}

	t = 0.0;
	calcul_de_U();
	om = 0.0;
	Cr = 0.0;
	teta = 0.0;
	Rbf = 0.0;
	isa1 = x[1] = 0.0;
	isb1 = x[2] = 0.0;
	isc1 = x[3] = 0.0;
	ira = x[4] = 0.0;
	irb = x[5] = 0.0;
	irc = x[6] = 0.0;
	om = x[7] = 0.0;
	teta = x[8] = 0.0;
}

// Método Runge-kutta
void mas()
{
	vecteur Y0, dx, Y1, Y2, Y3, Y4, tempo, K1, K2, K3, K4, tempo_2, tempo_3;
	int l;

	calcul_de_U();

	calcul_vecteur(x, Y0);
	multvect_const(Y0, dt, K1);
	multvect_const(K1, 0.5, dx);
	addvect_vect(x, dx, Y1);

	calcul_vecteur(Y1, tempo);
	multvect_const(tempo, dt, K2);
	multvect_const(K2, 0.5, dx);
	addvect_vect(x, dx, Y3);

	calcul_vecteur(Y3, tempo);
	multvect_const(tempo, dt, K3);
	addvect_vect(x, K3, Y4);

	calcul_vecteur(Y4, tempo);
	multvect_const(tempo, dt, K4);

	multvect_const(K2, 2, tempo);
	for (l = 1; l < rang; l++) K2[l] = tempo[l];

	multvect_const(K3, 2, tempo);
	for (l = 1; l < rang; l++) K3[l] = tempo[l];

	addvect_vect(K1, K2, tempo);
	addvect_vect(K3, K4, tempo_2);
	addvect_vect(tempo, tempo_2, tempo_3);

	multvect_const(tempo_3, 1 / 6.0, dx);
	addvect_vect(x, dx, tempo);

	for (l = 1; l < rang; l++) x[l] = tempo[l];

	isa1 = x[1];
	isb1 = x[2];
	isc1 = x[3];
	ira = x[4];
	irb = x[5];
	irc = x[6];
	om = x[7];
	teta = x[8];

	/* Cálculo do conjugado duas formas equivalentes*/

	te = -(R[7][4] * ira + R[7][5] * irb + R[7][6] * irc);

	// Para determinação da tensão entre os neutros
	// Cálculo da tensao nos bornes da carga e tensao de neutro
	/*
	double somme_dLI1, somme_L_dI1;
	double somme_dLI2, somme_L_dI2;
	double somme_dLI3, somme_L_dI3;

	somme_dLI1 = 0.0 ;
	somme_L_dI1 = 0.0 ;
	somme_dLI2 = 0.0 ;
	somme_L_dI2 = 0.0 ;
	somme_dLI3 = 0.0 ;
	somme_L_dI3 = 0.0 ;

		//FASE-1
	   for(l=5;l<=7;l++)
	   somme_L_dI1 +=  A[1][l] * dx[l];
	   for(l=5;l<=7;l++)
	   somme_dLI1 += b[1][l] * x[l] ;
	   somme_L_dI1 = somme_L_dI1 + A[1][1]*dx[1]+ A[1][4]*dx[4];
	   vsa_charge1 = b[1][1] * x[1] + b[1][4]*x[4] + somme_dLI1 + somme_L_dI1/hm2;
	   //vsa_charge =  Rs * isa +  somme_dLI * om  + somme_L_dI / dt  ;

		//FASE-3
	   for(l=5;l<=7;l++)
	   somme_L_dI2 +=  A[2][l] * dx[l];
	   for(l=5;l<=7;l++)
	   somme_dLI2 += b[2][l] * x[l] ;
	   somme_L_dI2 = somme_L_dI2 + A[2][2]*dx[2] + A[2][4]*dx[4];
	   vsa_charge2 = b[2][2] * x[2] + somme_dLI2 + somme_L_dI2/hm2;

	   //FASE-5
	   for(l=5;l<=7;l++)
	   somme_L_dI3 +=  A[3][l] * dx[l];
	   for(l=5;l<=7;l++)
	   somme_dLI3 += b[3][l] * x[l] ;
	   somme_L_dI3 = somme_L_dI3 + A[3][3]*dx[3]+ A[3][4]*dx[4];
	   vsa_charge3 = b[3][3] * x[3] + somme_dLI3 + somme_L_dI3/hm2;

	   //vn1 = vsa_charge1 - va1;
	   vn1 = (vsa_charge1 + vsa_charge2 + vsa_charge3)/3.;
	*/
}

/* Programa principal */
int wmain()
{

	init();    /* x est rempli à 0 */
	int l, c;
	FILE* fic;




	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	/* denominação do arquivo de saída */
	//fic=fopen("s000c30c100bi11_30A.dat","w");
	fic = fopen("joao_pedro.dat", "w");
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	results = mmake(NiterG, 7);

	results->mat[i][0] = t;
	results->mat[i][1] = om;
	results->mat[i][2] = isa1;
	results->mat[i][3] = isb1;
	results->mat[i][4] = isc1;
	results->mat[i][5] = isdef;
	results->mat[i][6] = te;


	// Partida com carga para acelerar a simulação

	Cr = 4.1;// Conjugado de carga nominal       = 4.1 Nm <<=======


  //  Cr=0.0;// Conjugado a vazio (desconectado) = 0.0 Nm <<=======
  //  Cr=0.5;// Conjugado a vazio (conectado) = 0.5 Nm <<=======
	for (i = 1; i < NiterG; i++)
	{
		for (j = 1; j <= Niter; j++)
		{
			mas();
			t = t + dt;
		}

		results->mat[i][0] = t;
		results->mat[i][1] = om;
		results->mat[i][2] = isa1;
		results->mat[i][3] = isb1;
		results->mat[i][4] = isc1;
		results->mat[i][5] = isdef;
		results->mat[i][6] = te;   // conjugado eletromagnetico


		if (_kbhit() && (_getch() == 27)) i = NiterG;

	}

	puts("acabou !");

	for (l = 1; l < 8; l++) printf("%lf \n", x[l]);
	printf("t= %lf\n", t);


	/* Aramazenamento dos resultados no arquivo de saída  */

	for (i = 0; i < NiterG; i++)
	{
		fprintf(fic, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",
			results->mat[i][0],
			results->mat[i][1],
			results->mat[i][2],
			results->mat[i][3],
			results->mat[i][4],
			results->mat[i][5],
			results->mat[i][6]);
	}

	fclose(fic);
	mfree(results);

}
