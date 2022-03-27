/*  M�quina ass�ncrona trif�sica */
/*  Modelo para componente fundamental (sem harm�nicas espaciais)  */
/*  Simula��o  de defeito no estator (curto-cir3uito entre espir1s) e assimetrias no rotor*/
#define _CRT_SECURE_NO_WARNINGS
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "routines.h"

#define pi     3.14159265358979323846
#define NiterG 100000    //130000;4000;300000(com conversor) /* N�mero total de pontos da simula��o */
#define Niter  10       //20 ;200;40(com conversor)         /* N�mero loops antes de pegar um  ponto */
#define dt     0.00001 //0.0000025//0.000005 //0.000001 (com conversor) /* DeltaT de simula��o dividido por Niter  0.00001 com Niter = 10 para ter Fa = 10kHz -freq. de amostragem*/
#define Umax   311  //1600 //439.82//1600 //311=220*sqrt(2) 129*sqrt(2)=179,61 /*Tens�o m�xima nos bornes da carga*/

int i, j;

// Dados de Placa do Motor de Indu��o Trif�sico WEG 60
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
const double Rs = 12.5;         /*  resist�ncia do estator    */
//const double Lls = 0.0383;      /*  Indut�ncia de dispers�o do estator  lsl  */
const double Lls = 0.007;
double Lms = 0.4506;      /*  Indut�ncia principal (pr�pia) do estator  Lsp  */
double Lsr = 0.4506;      /*  Indut�ncia m�tua estator-rotor  Msr */
const double Rr = 8.9;           /*  resist�ncia do rotor   */
const double Llr = 0.0383;      /*  Indut�ncia de dispers�o do rotor  lrl */
const double Lmr = 0.4506;      /*  Indut�ncia principal do rotor  Lrp */
const double Lsc = 1.5* Lms + Lls; /*  Indut�ncia c�clica do estator   */
const double Lrc = 1.5* Lmr + Llr; /*  Indut�ncia c�clica do rotor   */
//const double p = 1.0;           /*  N�mero de pares de p�los   */
const double p = 2.0;           /*  N�mero de pares de p�los   */
const double Jt = 0.023976 * 1.;  /*  Momento de In�rcia   */
const double fv = 1. * 0.0014439; /*  coeficiente de atrito din�mico   */

//Parâmetros do motor para novas equações (Obtidos a partir de artigo pesquisado)

const double Ns = 156; /* Número de voltas da bobina do estator */
const double rs = 1.5; /* Resitencia do estator */
const double Ls1 = 0.007; /* Indutância de dispersão do estator */

const double nb = 6; /* Número de barras do rotor */
const double Rb = 0.000096940036 ; /*  resistência das barras   */
const double Re = 0.000005 ;	/* resistencia do endring */
const double Lb = 0.00000028 ; /* Auto indutância da barra do rotor */
const double Le = 0.000000036; /* Auto indutância do endring */

const double r = 0.070 ; /* Raio medio do entreferro */
const double g = 0.00028 ; /* Entreferro */
const double l = 0.120 ; /* Comprimento efetivo do rotor */

const double mu  = 0.00000125663;



// **********************************************************************

// Maq. weg 60 (1 cv, lab.)
// Dados de cat�logo e ensaios a vazio e curto-cir3uito realizados pelo Claudio
//const double Rs = 7.65;         /*  resist�ncia do estator    */
//const double lsl = 0.0241;      /*  Indut�ncia de dispers�o do estator   */
//const double Lsp = 0.4709;      /*  Indut�ncia principal (pr�pia) do estator    */
//const double Msr = 0.4709;      /*  Indut�ncia m�tua estator-rotor   */
//const double Rr = 8.166;        /*  resist�ncia do rotor   */
//const double lrl = 0.0257;      /*  Indut�ncia de dispers�o do rotor   */
//const double Lrp = 0.4709;      /*  Indut�ncia principal do rotor   */
//const double Lsc = 1.5*Lsp+lsl; /*  Indut�ncia c�clica do estator   */
//const double Lrc = 1.5*Lrp+lrl; /*  Indut�ncia c�clica do rotor   */
//const double p = 2.0;           /*  N�mero de pares de p�los   */
////const double Jt = 0.00294;      /*  Momento de In�rcia do motor  */
//const double Jt = 0.023976;  /*  Momento de In�rcia   motor+carga   */
////const double fv = 0.0008169;    /*  coeficiente de atrito din�mico do motor  */
//const double fv = 0.0014439; /*  coeficiente de atrito din�mico  motor+carga  */

// **********************************************************************

double oms, Cr, teta, va1, vb1, vc1, isa1, isb1, isc1, ir1, ir2, ir3, ir4, ir5, ir6, ire, t, om, te, Rbf;

/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii*/

double rqi3 = 0.5774;
//
double vsdef, isdef, vl13, vl15, vl42, vl62, vx, ix, zx;
double vsa_charge1, vsa_charge2, vsa_charge3, vn1, vn1n2;
double vsa_charge4, vsa_charge5, vsa_charge6, vn2;
double ps_alim, ps_mot;

double Xa1, Xa2, Xa3, Xa4, Xa5, Xa6;
double Xb1, Xb2, Xb3, Xb4, Xb5, Xb6;
double Xc1, Xc2, Xc3, Xc4, Xc5, Xc6;
double K1, K2, K3, K4, K5, K6;

matrice dL, A, b, inv_L, B;
vecteur x, U, BU, Ad;

double R_r[7][7];
double L_r[7][7];
double L_sr[3][7];
double dL_sr[3][7];
double L_rs[7][3];
double dL_rs[7][3];
double R_s[3][3];
double L_s[3][3];

double L[13][13];
double R[13][13];

int tam = 7;

/* M�todo espec�fico de invers�o de Matriz SVD */
vecteur w;
matrice u, v, inva;

FILE* fichier;

matrix* results;

//void init_angles(void) 
//{
//
//	Xa1 = p * (teta + (pi / 6));
//	Xa2 = p * (teta + (pi / 3) + (pi / 6));
//	Xa3 = p * (teta + 2 * (pi / 3) + (pi / 6));
//	Xa4 = p * (teta + 3 * (pi / 3) + (pi / 6));
//	Xa5 = p * (teta + 4 * (pi / 3) + (pi / 6));
//	Xa6 = p * (teta + 5 * (pi / 3) + (pi / 6));
//
//	Xb1 = p * (teta + (pi / 6)) - 2*pi/3;
//	Xb2 = p * (teta + (pi / 3) + (pi / 6)) - 2 * pi / 3;
//	Xb3 = p * (teta + 2 * (pi / 3) + (pi / 6)) - 2 * pi / 3;
//	Xb4 = p * (teta + 3 * (pi / 3) + (pi / 6)) - 2 * pi / 3;
//	Xb5 = p * (teta + 4 * (pi / 3) + (pi / 6)) - 2 * pi / 3;
//	Xb6 = p * (teta + 5 * (pi / 3) + (pi / 6)) - 2 * pi / 3;
//
//	Xc1 = p * (teta + (pi / 6)) + 2 * pi / 3;
//	Xc2 = p * (teta + (pi / 3) + (pi / 6)) + 2 * pi / 3;
//	Xc3 = p * (teta + 2 * (pi / 3) + (pi / 6)) + 2 * pi / 3;
//	Xc4 = p * (teta + 3 * (pi / 3) + (pi / 6)) + 2 * pi / 3;
//	Xc5 = p * (teta + 4 * (pi / 3) + (pi / 6)) + 2 * pi / 3;
//	Xc6 = p * (teta + 5 * (pi / 3) + (pi / 6)) + 2 * pi / 3;
//}

void assemble_Rr(void)
{

	double a = 2 * (Rb + Re);
	double b = -Rb;
	double c = -Re;
	double d = tam * Re;

	for (int m = 0; m < tam; m++) {
		for (int n = 0; n < tam; n++) {
			if (m == 0 && n == tam - 2) {
				R_r[m][n] = b;
			}
			else
				if (m == tam - 2 && n == 0) {
					R_r[m][n] = b;
				}
				else if (m == n && m != tam - 1) {
					R_r[m][n] = a;
				}
				else if (m + 1 == n && n != tam - 1) {
					R_r[m][n] = b;
				}
				else if (m - 1 == n && m != tam - 1) {
					R_r[m][n] = b;
				}
				else if (n == tam - 1) {
					R_r[m][n] = c;
				}
				else if (m == tam - 1) {
					R_r[m][n] = c;
				}
				else {
					R_r[m][n] = 0;
				}
		}
	}
	R_r[tam - 1][tam - 1] = d;
}

void assemble_Lr(void)
{

	double alpha = 2 * pi / tam;
	double Lkk = (mu * l * r / g) * (1 - alpha / (2 * pi)) * alpha;
	double Lki = -((mu * l * r) / g) * ((alpha * alpha) / (2 * pi));
	double Ld = Lkk + 2 * (Le + Lb);
	double Lt = Lki - Lb;

	double a = Ld;
	double b = Lt;
	double c = -Le;
	double d = (tam - 1) * Le;

	for (int m = 0; m < tam; m++) {
		for (int n = 0; n < tam; n++) {
			if (m == 0 && n == tam - 2) {
				L_r[m][n] = b;
			}
			else
				if (m == tam - 2 && n == 0) {
					L_r[m][n] = b;
				}
				else if (m == n && m != tam - 1) {
					L_r[m][n] = a;
				}
				else if (m + 1 == n && n != tam - 1) {
					L_r[m][n] = b;
				}
				else if (m - 1 == n && m != tam - 1) {
					L_r[m][n] = b;
				}
				else if (n == tam - 1) {
					L_r[m][n] = c;
				}
				else if (m == tam - 1) {
					L_r[m][n] = c;
				}
				else {
					L_r[m][n] = Lki;
				}
		}
	}
	L_r[tam - 1][tam - 1] = d;
}

void assemble_Lsr() { //TODO: corrigir 

	double alpha = 2 * pi / (tam - 1);

	for (int m = 0; m < 3; m++) {
		for (int n = 0; n < tam - 1; n++) {
			if (m == 0) {
				L_sr[m][n] = Lsr * cos(p * (teta + n * alpha + alpha / 2));
			}
			if (m == 1) {
				L_sr[m][n] = Lsr * cos(p * (teta + n * alpha + alpha / 2) - 2 * pi / 3);
			}
			if (m == 2) {
				L_sr[m][n] = Lsr * cos(p * (teta + n * alpha + alpha / 2) + 2 * pi / 3);
			}
		}
	}
}

void assemble_dLsr() {

	double alpha = 2 * pi / (tam - 1);
	double K = -p * Lsr * om;

	for (int m = 0; m < 3; m++) {
		for (int n = 0; n < tam - 1; n++) {
			if (m == 0) {
				dL_sr[m][n] = K * sin(p * (teta + n * alpha + alpha / 2));
			}
			if (m == 1) {
				dL_sr[m][n] = K * sin(p * (teta + n * alpha + alpha / 2) - 2 * pi / 3);
			}
			if (m == 2) {
				dL_sr[m][n] = K * sin(p * (teta + n * alpha + alpha / 2) + 2 * pi / 3);
			}
		}
	}
}

void assemble_Lrs(void) {
	for (int m = 0; m < tam; m++) {
		for (int n = 0; n < 3; n++) {
			L_rs[m][n] = L_sr[n][m];
		}
	}
}

void assemble_dLrs(void) {
	for (int m = 0; m < tam; m++) {
		for (int n = 0; n < 3; n++) {
			dL_rs[m][n] = dL_sr[n][m];
		}
	}
}

void assemble_R_s(void) {
	R_s[0][0] = R_s[1][1] = R_s[2][2] = Rs;
}

void assemble_L_s(void) {
	L_s[0][0] = L_s[1][1] = L_s[2][2] = Lls + Lms;
	L_s[0][1] = L_s[0][2] = L_s[1][2] = L_s[1][0] = L_s[2][0] = L_s[2][1] = -0.5 * Lms;
}

//void calcul_de_R(void)
//{
//	double c;
//	double sinteta;
//	double sinteta_pos;
//	double sinteta_neg;
//
//	c = -Lsr * x[8] * p; //TODO: Verificar
//	sinteta = sin(p * teta);
//	sinteta_pos = sin(p * teta + 2.0 * pi / 3.0);
//	sinteta_neg = sin(p * teta - 2.0 * pi / 3.0);
//
//	K1 = Lsr * p * (x[1] * sin(Xa1) + x[2] * sin(Xb1) + x[3] * sin(Xc1));
//	K2 = Lsr * p * (x[1] * sin(Xa2) + x[2] * sin(Xb2) + x[3] * sin(Xc2));
//	K3 = Lsr * p * (x[1] * sin(Xa3) + x[2] * sin(Xb3) + x[3] * sin(Xc3));
//	K4 = Lsr * p * (x[1] * sin(Xa4) + x[2] * sin(Xb4) + x[3] * sin(Xc4));
//	K5 = Lsr * p * (x[1] * sin(Xa5) + x[2] * sin(Xb5) + x[3] * sin(Xc5));
//	K6 = Lsr * p * (x[1] * sin(Xa6) + x[2] * sin(Xb6) + x[3] * sin(Xc6));
//
//	
//	R[1][1] = R[2][2] = R[3][3] = rs;
//
//	R[4][4] = R[5][5] = R[6][6] = R[7][7] = R[8][8] = R[9][9] = 2*(Rb + Re);
//
//
//	R[1][4] = R[4][1] = c * sin(p*(teta));
//	R[1][5] = R[5][1] = c * sin(p*teta+2*pi/3);
//	R[1][6] = R[6][1] = c * sin(p*teta-2*pi/3.0);
//
//	R[1][4] = R[4][1] = c * sin(Xa1);
//	R[1][5] = R[5][1] = c * sin(Xa2);
//	R[1][6] = R[6][1] = c * sin(Xa3);
//	R[1][7] = R[7][1] = c * sin(Xa4);
//	R[1][8] = R[8][1] = c * sin(Xa5);
//	R[1][9] = R[9][1] = c * sin(Xa6);
//
//	R[2][4] = R[4][2] = c * sin(Xb1);
//	R[2][5] = R[5][2] = c * sin(Xb2);
//	R[2][6] = R[6][2] = c * sin(Xb3);
//	R[2][7] = R[7][2] = c * sin(Xb4);
//	R[2][8] = R[8][2] = c * sin(Xb5);
//	R[2][9] = R[9][2] = c * sin(Xb6);
//
//	R[3][4] = R[4][3] = c * sin(Xc1);
//	R[3][5] = R[5][3] = c * sin(Xc2);
//	R[3][6] = R[6][3] = c * sin(Xc3);
//	R[3][7] = R[7][3] = c * sin(Xc4);
//	R[3][8] = R[8][3] = c * sin(Xc5);
//	R[3][9] = R[9][3] = c * sin(Xc6);
//
//	R[4][5] = R[4][9] = R[5][4] = R[5][6] = R[6][5] = R[6][7] = R[7][6] = R[7][8] = R[8][7] = R[8][9] = R[9][4] = R[9][8] = -Rb;
//
//	R[10][4] = R[10][5] = R[10][6] = R[10][7] = R[10][8] = R[10][9] = -Re;
//	R[4][10] = R[5][10] = R[6][10] = R[7][10] = R[8][10] = R[9][10] = -Re;
//	R[11][4] = K1;
//	R[11][5] = K2;
//	R[11][6] = K3;
//	R[11][7] = K4;
//	R[11][8] = K5;
//	R[11][9] = K6;
//
//
//	R[10][10] = nb*Re;
//	R[11][11] = fv;
//	R[12][11] = -1;
//}

void assemble_L(void) {

	assemble_L_s();
	assemble_Lsr();
	assemble_Lrs();
	assemble_Lr();
	int m, n;

	for ( m = 1; m <= 3; m++) {
		for ( n = 1; n <= 3; n++) {
			L[m][n] = L_s[m - 1][n - 1];
		}
	}
	for (m = 1; m <= 3; m++) {
		for (n = 4; n <= tam + 3; n++) {
			L[m][n] = L_sr[m - 1][n - 4];
		}
	}
	for (m = 4; m <= tam + 3; m++) {
		for (n = 1; n <= 3; n++) {
			L[m][n] = L_rs[m - 4][n - 1];
		}
	}
	for (m = 4; m <= tam + 3; m++) {
		for (n = 4; n <= tam + 3; n++) {
			L[m][n] = L_r[m - 4][n - 4];
		}
	}
	L[tam + 3 + 1][tam + 3 + 1] = Jt;
	L[tam + 3 + 2][tam + 3 + 2] = 1.0;
}

void assemble_R(void) {

	assemble_R_s();
	assemble_dLsr();
	assemble_dLrs();
	assemble_Rr();
	double alpha = 2 * pi / tam;
	int m, n;

	for ( m = 1; m <= 3; m++) {
		for ( n = 1; n <= 3; n++) {
			R[m][n] = R_s[m - 1][n - 1];
		}
	}
	for (m = 1; m <= 3; m++) {
		for (n = 4; n <= tam + 3; n++) {
			R[m][n] = dL_sr[m - 1][n - 4];
		}
	}
	for (m = 4; m <= tam + 3; m++) {
		for (n = 1; n <= 3; n++) {
			R[m][n] = dL_rs[m - 4][n - 1];
		}
	}
	for (m = 4; m <= tam + 3; m++) {
		for (n = 4; n <= tam + 3; n++) {
			R[m][n] = R_r[m - 4][n - 4];
		}
	}

	for (n = 4; n <= tam + 2; n++) {
		double X1 = p * (teta + n - 4 * alpha + alpha / 2);
		double X2 = p * (teta + n - 4 * alpha + alpha / 2) - 2 * pi / 3;
		double X3 = p * (teta + n - 4 * alpha + alpha / 2) + 2 * pi / 3;
		R[tam + 4][n] = Lsr * p * (isa1 * sin(X1) + isb1 * sin(X2) + isc1 * sin(X3));
	}
	R[tam + 3 + 1][tam + 3 + 1] = fv;
	R[tam + 3 + 2][tam + 3 + 1] = -1.0;
}

//void calcul_de_L(void)
//{
//	double costeta;
//	double costeta_neg;
//	double costeta_pos;
//	double Lkk;
//	double Ld, Lt, Lki;
//
//	Lkk = (mu*l*r/g)*(1-(pi/3)/(2*pi))*(pi/3);
//	Ld = Lkk + 2 * (Le + Lb);
//	Lki = -((mu * l * r) / g) * (((pi / 3) * (pi / 3)) / (2 * pi));
//	Lt = Lki - Lb;
//
//	L[1][1] = L[2][2] = L[3][3] = Lls + Lms;
//	L[1][2] = L[1][3] = L[2][3] = L[2][1] = L[3][1] = L[3][2] = -0.5 * Lms;
//
//	L[1][4] = L[4][1] = Lsr * cos(Xa1);
//	L[1][5] = L[5][1] = Lsr * cos(Xa2);
//	L[1][6] = L[6][1] = Lsr * cos(Xa3);
//	L[1][7] = L[7][1] = Lsr * cos(Xa4);
//	L[1][8] = L[8][1] = Lsr * cos(Xa5);
//	L[1][9] = L[9][1] = Lsr * cos(Xa6);
//
//	L[2][4] = L[4][2] = Lsr * cos(Xb1);
//	L[2][5] = L[5][2] = Lsr * cos(Xb2);
//	L[2][6] = L[6][2] = Lsr * cos(Xb3);
//	L[2][7] = L[7][2] = Lsr * cos(Xb4);
//	L[2][8] = L[8][2] = Lsr * cos(Xb5);
//	L[2][9] = L[9][2] = Lsr * cos(Xb6);
//
//	L[3][4] = L[4][3] = Lsr * cos(Xc1);
//	L[3][5] = L[5][3] = Lsr * cos(Xc2);
//	L[3][6] = L[6][3] = Lsr * cos(Xc3);
//	L[3][7] = L[7][3] = Lsr * cos(Xc4);
//	L[3][8] = L[8][3] = Lsr * cos(Xc5);
//	L[3][9] = L[9][3] = Lsr * cos(Xc6);
//
//	L[4][4] = L[5][5] = L[6][6] = L[7][7] = L[8][8] = L[9][9] = Ld;
//
//	L[4][5] = L[4][9] = L[5][4] = L[5][6] = L[6][5] = L[6][7] = L[7][6] = L[7][8] = L[8][7] = L[8][9] = L[9][4] = L[9][8] = Lt;
//
//	L[4][6] = L[4][7] = L[4][8] = L[5][7] = L[5][8] = L[5][9] = L[6][4] = L[6][8] = L[6][9] = L[7][4] = L[7][5] = L[7][9] = L[8][4] = L[8][5] = L[8][6] = L[9][5] = L[9][6] = L[9][7] = Lki;
//
//	L[4][10] = L[5][10] = L[6][10] = L[7][10] = L[8][10] = L[9][10] = -Le;
//	L[10][4] = L[10][5] = L[10][6] = L[10][7] = L[10][8] = L[10][9] = -Le;
//
//	L[10][1] = L[10][2] = L[10][3] = L[11][1] = L[11][2] = L[11][3] = L[12][1] = L[12][2] = L[12][3] = 0.0;
//
//	L[11][4] = L[11][5] = L[11][6] = L[11][7] = L[11][8] = L[11][9] = L[11][10] = 0.0;
//	L[12][4] = L[12][5] = L[12][6] = L[12][7] = L[12][8] = L[12][9] = L[12][10] = 0.0;
//	L[10][10] = nb*Le;
//
//	L[1][10] = L[2][10] = L[3][10] = 0;
//	L[1][11] = L[2][11] = L[3][11] = L[4][11] = L[5][11] = L[6][11] = L[7][11] = L[8][11] = L[9][11] = L[10][11] = L[12][11] = 0;
//	L[1][12] = L[2][12] = L[3][12] = L[4][12] = L[5][12] = L[6][12] = L[7][12] = L[8][12] = L[9][12] = L[10][12] = L[11][12] = 0;
//
//
//	L[11][11] = Jt;
//	L[12][12] = 1.0;
//
//}

// Gera��o da Tens�o de Alimenta��o
void calcul_de_U(void)
{
	double v1, v3, v5, v7, vdc, vdc2, fs;
	fs = 50.0;
	oms = 2 * pi * fs;

	va1 = Umax * cos(oms * t);
	vb1 = Umax * cos(oms * t + 2.0 * pi / 3.0);
	vc1 = Umax * cos(oms * t - 2.0 * pi / 3.0);

	vsdef = 0.0;

	// Entrada do vetor de Tens�o
	U[1] = va1;
	U[2] = vb1;
	U[3] = vc1;
	U[4] = U[5] = U[6] = U[7] = U[8] = U[9] = 0.0;  // Tens�es Rot�ricas
	U[10] = 0.0; // Tensão endring
	U[11] = -Cr;            // Conjugado de Carga (resistente)
	U[12] = 0.0;
}


/*****************************************************************************/


void calcul_vecteur(double d[rang], double c[rang])
{
	// As fun��es de invers�o, multiplica��o e subtra��o de matrizes est�o no programa routines
	//calcul_de_L();          /* Calcula-se a matriz L */

	//calcul_de_R();          /* Matriz R */

	assemble_L();
	assemble_R();

	inv(L, inv_L);  /* Invers�o de L */

	// M�todo do Piv� de Gauss (foi usado aqui pq � mais r�pido)

	//multmat_mat(inv_L, R, A);
	//multmat_vect(A, d, Ad);  // A * x
	//multmat_vect(inv_L, U, BU);  // B * U   B = inv_L;
	//sousvect_vect(BU, Ad, c);

	multmat_vect(R, d, BU);
	sousvect_vect(U, BU, BU);
	multmat_vect(inv_L, BU, c);


}

// Inicializa��o das Vari�veis
void init()
{

	int l, c;

	for (l = 0; l < rang; l++) x[l] = 0;
	for (l = 0; l < rang; l++)               // Inicializa-se R com 0
		for (c = 0; c < rang; c++)
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

	isa1 = x[1] = 0;
	isb1 = x[2] = 0;
	isc1 = x[3] = 0;
	ir1 = x[4] = 0;
	ir2 = x[5] = 0;
	ir3 = x[6] = 0;
	ir4 = x[7] = 0;
	ir5 = x[8] = 0;
	ir6 = x[9] = 0;
	ire = x[10] = 0;
	om = x[11] = 0;
	teta = x[12] = 0;

}

// M�todo Runge-kutta
void mas()
{
	//vecteur Y0, dx, Y1, Y2, Y3, Y4, tempo, K1, K2, K3, K4, tempo_2, tempo_3;
	//int l;

	//calcul_de_U();

	//calcul_vecteur(x, Y0);
	//multvect_const(Y0, dt, K1);
	//multvect_const(K1, 0.5, dx);
	//addvect_vect(x, dx, Y1);

	//calcul_vecteur(Y1, tempo);
	//multvect_const(tempo, dt, K2);
	//multvect_const(K2, 0.5, dx);
	//addvect_vect(x, dx, Y3);

	//calcul_vecteur(Y3, tempo);
	//multvect_const(tempo, dt, K3);
	//addvect_vect(x, K3, Y4);

	//calcul_vecteur(Y4, tempo);
	//multvect_const(tempo, dt, K4);

	//multvect_const(K2, 2, tempo);
	//for (l = 1; l < rang; l++) K2[l] = tempo[l];

	//multvect_const(K3, 2, tempo);
	//for (l = 1; l < rang; l++) K3[l] = tempo[l];

	//addvect_vect(K1, K2, tempo);
	//addvect_vect(K3, K4, tempo_2);
	//addvect_vect(tempo, tempo_2, tempo_3);

	//multvect_const(tempo_3, 1 / 6.0, dx);
	//addvect_vect(x, dx, tempo);

	//for (l = 1; l < rang; l++) x[l] = tempo[l];

	vecteur Y0, dx, Y1, Y2, Y3, tempo;
	int l;

	calcul_de_U();

	calcul_vecteur(x, Y0);
	multvect_const(Y0, dt / 2.0, dx);
	addvect_vect(x, dx, Y1);

	calcul_vecteur(Y1, tempo);
	for (l = 1; l < rang; l++) Y1[l] = tempo[l];
	multvect_const(Y1, dt / 2.0, dx);
	addvect_vect(x, dx, Y2);

	calcul_vecteur(Y2, tempo);
	for (l = 1; l < rang; l++) Y2[l] = tempo[l];
	multvect_const(Y2, dt, dx);
	addvect_vect(x, dx, Y3);

	calcul_vecteur(Y3, tempo);
	for (l = 1; l < rang; l++) Y3[l] = tempo[l];

	multvect_const(Y1, 2.0, tempo);
	for (l = 1; l < rang; l++) Y1[l] = tempo[l];

	multvect_const(Y2, 2.0, tempo);
	for (l = 1; l < rang; l++) Y2[l] = tempo[l];
	for (l = 1; l < rang; l++) dx[l] = Y0[l];
	addvect_vect(dx, Y1, tempo);
	for (l = 1; l < rang; l++) dx[l] = tempo[l];
	addvect_vect(dx, Y2, tempo);
	for (l = 1; l < rang; l++) dx[l] = tempo[l];
	addvect_vect(dx, Y3, tempo);
	for (l = 1; l < rang; l++) dx[l] = tempo[l];

	multvect_const(dx, dt / 6.0, tempo);
	for (l = 1; l < rang; l++) dx[l] = tempo[l];
	addvect_vect(x, dx, tempo);

	for (l = 1; l < rang; l++) x[l] = tempo[l];

	isa1 = x[1];
	isb1 = x[2];
	isc1 = x[3];
	ir1 = x[4];
	ir2 = x[5];
	ir3 = x[6];
	ir4 = x[7];
	ir5 = x[8];
	ir6 = x[9];
	ire = x[10];
	om = x[11];
	teta = x[12];

	/* C�lculo do conjugado duas formas equivalentes*/

	te = -(R[11][4] * ir1 + R[11][5] * ir2 + R[11][6] * ir3 + R[11][7] * ir4 + R[11][8] * ir5 + R[11][9] * ir6);

	// Para determina��o da tens�o entre os neutros
	// C�lculo da tensao nos bornes da carga e tensao de neutro
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

	Lms = (Ns / (2 * p)) * (Ns / (2 * p)) * ((pi * mu * l * r) / g);
	Lsr = 4 / pi * Lms / Ns * sin(p * pi / 6);


	init();    /* x est rempli � 0 */
	int l, c;
	FILE* fic;




	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	/* denomina��o do arquivo de sa�da */
	//fic=fopen("s000c30c100bi11_30A.dat","w");
	fic = fopen("joao_pedro_6_barras.dat", "w");
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	results = mmake(NiterG, 11);

	results->mat[i][0] = t;
	results->mat[i][1] = om;
	results->mat[i][2] = isa1;
	results->mat[i][3] = isb1;
	results->mat[i][4] = isc1;
	results->mat[i][5] = isdef;
	results->mat[i][6] = te;
	results->mat[i][7] = teta;
	results->mat[i][8] = ir1;
	results->mat[i][9] = ir2;
	results->mat[i][10] = ir3;


	// Partida com carga para acelerar a simula��o

	Cr = 2.38;// Conjugado de carga nominal       = 4.1 Nm <<=======
	//Cr = 1;

  //  Cr=0.0;// Conjugado a vazio (desconectado) = 0.0 Nm <<=======
    //Cr=1;// Conjugado a vazio (conectado) = 0.5 Nm <<=======
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
		results->mat[i][7] = teta;
		results->mat[i][8] = ir1;
		results->mat[i][9] = ir2;
		results->mat[i][10] = ir3;

		if (_kbhit() && (_getch() == 27)) i = NiterG;

	}

	puts("acabou !");

	for (l = 1; l < 11; l++) printf("%lf \n", x[l]);
	printf("t= %lf\n", t);


	/* Aramazenamento dos resultados no arquivo de sa�da  */

	for (i = 0; i < NiterG; i++)
	{
		fprintf(fic, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",
			results->mat[i][0],
			results->mat[i][1],
			results->mat[i][2],
			results->mat[i][3],
			results->mat[i][4],
			results->mat[i][5],
			results->mat[i][6],
			results->mat[i][7],
			results->mat[i][8],
			results->mat[i][9],
			results->mat[i][10]);
	}

	fclose(fic);
	mfree(results);

}
