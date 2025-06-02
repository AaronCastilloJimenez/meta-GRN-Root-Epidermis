/*
This program is an update to the simulation of the coupling of genetic regulatory networks in a two-dimensional domain.
This domain is composed of a web of cells representing epidermal cells of the root of Arabidopsis thaliana.
Each cell contains a network composed of the proteins wer, myb23, gl3/egl3, ttg1, cpc, try, etc1, ttg2, gl2, zfp5, wrky75, scm;
the hormones ethylene, aux, and cytokinins; and the activator AC and inhibitor IC protein complexes.
The cells communicate by simulating cpc and gl3/egl3 diffusion.
*/

#include <stdio.h>
#include <math.h>   // For pow() and sqrt()
#include <stdlib.h> // For malloc(), free(), rand(), srand()
#include <time.h>   // For time()

// --- Constants ---
#define GRID_SIZE 24     // Dimension of the cell grid
#define MAX_ITERATIONS 85 // Default number of iterations per simulation

// --- Global parameters (user-defined or fixed conditions) ---
float DifI;         // Diffusion parameter for CPC
float DifA;         // Diffusion parameter for GL3/EGL3
int degcpc;         // Degradation rate for CPC
int deggl3egl3;     // Degradation rate for GL3/EGL3
int citoquininas;   // Cytokinins level
int phos;           // Phosphate level
int salt;           // Salt stress level
int CO2;            // CO2 level

// --- Function Prototypes ---
void Iniciales(int werx[][GRID_SIZE], int myb23x[][GRID_SIZE], int gl3egl3x[][GRID_SIZE],
               int ttg1x[][GRID_SIZE], int cpcx[][GRID_SIZE], int tryx[][GRID_SIZE],
               int etc1x[][GRID_SIZE], int ttg2x[][GRID_SIZE], int gl2x[][GRID_SIZE],
               int zfp5x[][GRID_SIZE], int scmx[][GRID_SIZE], int jkdx[][GRID_SIZE],
               int *wrky75x_ptr, int *etx_ptr, int *auxx_ptr,
               float cpcn[][GRID_SIZE], float gl3egl3n[][GRID_SIZE],
               int gl3egl3Local[][GRID_SIZE], int gl3egl3Dif[][GRID_SIZE],
               int cpcLocal[][GRID_SIZE], int cpcDif[][GRID_SIZE]);

int H(float x);

void difundecpc(float cpcn[][GRID_SIZE], int cpcx[][GRID_SIZE], int cpcDif[][GRID_SIZE]);
void difundegl3egl3(float gl3egl3n[][GRID_SIZE], int gl3egl3x[][GRID_SIZE], int gl3egl3Dif[][GRID_SIZE]);

void Itera(int ite, int posicion[][GRID_SIZE],
           int werx[][GRID_SIZE], int myb23x[][GRID_SIZE], int gl3egl3x[][GRID_SIZE],
           int ttg1x[][GRID_SIZE], int cpcx[][GRID_SIZE], int tryx[][GRID_SIZE],
           int etc1x[][GRID_SIZE], int ttg2x[][GRID_SIZE], int gl2x[][GRID_SIZE],
           int zfp5x[][GRID_SIZE], int scmx[][GRID_SIZE], int jkdx[][GRID_SIZE],
           int *wrky75x_ptr, int *etx_ptr, int *auxx_ptr,
           int wery[][GRID_SIZE], int myb23y[][GRID_SIZE], int gl3egl3y[][GRID_SIZE],
           int ttg1y[][GRID_SIZE], int cpcy[][GRID_SIZE], int tryy[][GRID_SIZE],
           int etc1y[][GRID_SIZE], int ttg2y[][GRID_SIZE], int gl2y[][GRID_SIZE],
           int zfp5y[][GRID_SIZE], int scmy[][GRID_SIZE], int jkdy[][GRID_SIZE],
           int *wrky75y_ptr, int *ety_ptr, int *auxy_ptr,
           int AC[][GRID_SIZE], int IC[][GRID_SIZE],
           int cpcLocal[][GRID_SIZE], int gl3egl3Local[][GRID_SIZE],
           int cpcDif[][GRID_SIZE], int gl3egl3Dif[][GRID_SIZE],
           int cell[][GRID_SIZE], int *Hectopic_ptr, int *NHectopic_ptr);

void imprimeC(int cell[][GRID_SIZE]);
void imprimedata(int cell[][GRID_SIZE]);


// --- Function Implementations ---

void Iniciales(int werx[][GRID_SIZE], int myb23x[][GRID_SIZE], int gl3egl3x[][GRID_SIZE],
               int ttg1x[][GRID_SIZE], int cpcx[][GRID_SIZE], int tryx[][GRID_SIZE],
               int etc1x[][GRID_SIZE], int ttg2x[][GRID_SIZE], int gl2x[][GRID_SIZE],
               int zfp5x[][GRID_SIZE], int scmx[][GRID_SIZE], int jkdx[][GRID_SIZE],
               int *wrky75x_ptr, int *etx_ptr, int *auxx_ptr,
               float cpcn[][GRID_SIZE], float gl3egl3n[][GRID_SIZE],
               int gl3egl3Local[][GRID_SIZE], int gl3egl3Dif[][GRID_SIZE],
               int cpcLocal[][GRID_SIZE], int cpcDif[][GRID_SIZE]) {
    int i, j;
    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            werx[i][j] = rand() % 2;
            myb23x[i][j] = 0;
            gl3egl3x[i][j] = rand() % 2;
            gl3egl3Local[i][j] = rand() % 2;
            gl3egl3Dif[i][j] = rand() % 2;
            ttg1x[i][j] = 1; // Fixed
            cpcx[i][j] = rand() % 2;
            cpcLocal[i][j] = rand() % 2;
            cpcDif[i][j] = rand() % 2;
            tryx[i][j] = rand() % 2;
            etc1x[i][j] = rand() % 2;
            ttg2x[i][j] = rand() % 2;
            gl2x[i][j] = rand() % 2;
            zfp5x[i][j] = rand() % 2;
            scmx[i][j] = rand() % 2;
            jkdx[i][j] = rand() % 2;
            cpcn[i][j] = 0;
            gl3egl3n[i][j] = 0;
        }
    }
    *wrky75x_ptr = 1;
    *etx_ptr = 0;
    *auxx_ptr = 1;
}

// Threshold function
int H(float x) {
    if (x <= 0.25) return 0;
    if (x > 0.25 && x <= 1) return 1;
    return 2;
}

void difundecpc(float cpcn[][GRID_SIZE], int cpcx[][GRID_SIZE], int cpcDif[][GRID_SIZE]) {
    int i, j;
    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            int i_plus_1 = (i + 1) % GRID_SIZE;
            int i_minus_1 = (i - 1 + GRID_SIZE) % GRID_SIZE;
            int j_plus_1 = (j + 1) % GRID_SIZE;
            int j_minus_1 = (j - 1 + GRID_SIZE) % GRID_SIZE;

            cpcn[i][j] = (float)cpcx[i][j] + DifI * (
                (float)cpcx[i_plus_1][j] +
                (float)cpcx[i_minus_1][j] +
                (float)cpcx[i][j_plus_1] +
                (float)cpcx[i][j_minus_1] -
                4.0 * (float)cpcx[i][j]
            );
            cpcDif[i][j] = H(cpcn[i][j]);
        }
    }
}

void difundegl3egl3(float gl3egl3n[][GRID_SIZE], int gl3egl3x[][GRID_SIZE], int gl3egl3Dif[][GRID_SIZE]) {
    int i, j;
    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            int i_plus_1 = (i + 1) % GRID_SIZE;
            int i_minus_1 = (i - 1 + GRID_SIZE) % GRID_SIZE;
            int j_plus_1 = (j + 1) % GRID_SIZE;
            int j_minus_1 = (j - 1 + GRID_SIZE) % GRID_SIZE;
            
            gl3egl3n[i][j] = (float)gl3egl3x[i][j] + DifA * (
                (float)gl3egl3x[i_plus_1][j] +
                (float)gl3egl3x[i_minus_1][j] +
                (float)gl3egl3x[i][j_plus_1] +
                (float)gl3egl3x[i][j_minus_1] -
                4.0 * (float)gl3egl3x[i][j]
            );
            gl3egl3Dif[i][j] = H(gl3egl3n[i][j]);
        }
    }
}

void Itera(int ite, int posicion[][GRID_SIZE],
           int werx[][GRID_SIZE], int myb23x[][GRID_SIZE], int gl3egl3x[][GRID_SIZE],
           int ttg1x[][GRID_SIZE], int cpcx[][GRID_SIZE], int tryx[][GRID_SIZE],
           int etc1x[][GRID_SIZE], int ttg2x[][GRID_SIZE], int gl2x[][GRID_SIZE],
           int zfp5x[][GRID_SIZE], int scmx[][GRID_SIZE], int jkdx[][GRID_SIZE],
           int *wrky75x_ptr, int *etx_ptr, int *auxx_ptr,
           int wery[][GRID_SIZE], int myb23y[][GRID_SIZE], int gl3egl3y[][GRID_SIZE],
           int ttg1y[][GRID_SIZE], int cpcy[][GRID_SIZE], int tryy[][GRID_SIZE],
           int etc1y[][GRID_SIZE], int ttg2y[][GRID_SIZE], int gl2y[][GRID_SIZE],
           int zfp5y[][GRID_SIZE], int scmy[][GRID_SIZE], int jkdy[][GRID_SIZE],
           int *wrky75y_ptr, int *ety_ptr, int *auxy_ptr,
           int AC[][GRID_SIZE], int IC[][GRID_SIZE],
           int cpcLocal[][GRID_SIZE], int gl3egl3Local[][GRID_SIZE],
           int cpcDif[][GRID_SIZE], int gl3egl3Dif[][GRID_SIZE],
           int cell[][GRID_SIZE], int *Hectopic_ptr, int *NHectopic_ptr) {
    
    int i, j;
    int myb;
    int coin, rvar;

    int current_wrky75x = *wrky75x_ptr;
    int current_etx = *etx_ptr;
    int current_auxx = *auxx_ptr;

    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            // IC (Inhibitor Complex)
            if (cpcx[i][j] == 0 && tryx[i][j] == 0 && etc1x[i][j] == 0) IC[i][j] = 0;
            else if (cpcx[i][j] == 0 && tryx[i][j] == 0 && etc1x[i][j] == 1) IC[i][j] = 0;
            else if (cpcx[i][j] == 0 && tryx[i][j] == 1 && etc1x[i][j] == 0) IC[i][j] = 0;
            else if (cpcx[i][j] == 0 && tryx[i][j] == 1 && etc1x[i][j] == 1) IC[i][j] = 1;
            else if (cpcx[i][j] == 0 && tryx[i][j] == 2 && etc1x[i][j] == 0) IC[i][j] = 1;
            else if (cpcx[i][j] == 0 && tryx[i][j] == 2 && etc1x[i][j] == 1) IC[i][j] = 2;
            else if (cpcx[i][j] == 1 && tryx[i][j] == 0 && etc1x[i][j] == 0) IC[i][j] = 1;
            else if (cpcx[i][j] == 1 && tryx[i][j] == 0 && etc1x[i][j] == 1) IC[i][j] = 1;
            else if (cpcx[i][j] == 1 && tryx[i][j] == 1 && etc1x[i][j] == 0) IC[i][j] = 1;
            else if (cpcx[i][j] == 1 && tryx[i][j] == 1 && etc1x[i][j] == 1) IC[i][j] = 1;
            else if (cpcx[i][j] == 1 && tryx[i][j] == 2 && etc1x[i][j] == 0) IC[i][j] = 2;
            else if (cpcx[i][j] == 1 && tryx[i][j] == 2 && etc1x[i][j] == 1) IC[i][j] = 2;
            else if (cpcx[i][j] == 2 && tryx[i][j] == 0 && etc1x[i][j] == 0) IC[i][j] = 1;
            else if (cpcx[i][j] == 2 && tryx[i][j] == 0 && etc1x[i][j] == 1) IC[i][j] = 2;
            else if (cpcx[i][j] == 2 && tryx[i][j] == 1 && etc1x[i][j] == 0) IC[i][j] = 2;
            else if (cpcx[i][j] == 2 && tryx[i][j] == 1 && etc1x[i][j] == 1) IC[i][j] = 2;
            else if (cpcx[i][j] == 2 && tryx[i][j] == 2 && etc1x[i][j] == 0) IC[i][j] = 2;
            else if (cpcx[i][j] == 2 && tryx[i][j] == 2 && etc1x[i][j] == 1) IC[i][j] = 2;
            else IC[i][j] = 0;

            myb = werx[i][j] + myb23x[i][j];

            // AC (Activator Complex) - MBW (MYB-bHLH-WDR)
            AC[i][j] = 0; 
            if (gl3egl3x[i][j] == 0 || myb == 0) AC[i][j] = 0;
            else if (gl3egl3x[i][j] == 1 && ttg1x[i][j] == 0) AC[i][j] = 0;
            else if (gl3egl3x[i][j] == 1 && ttg1x[i][j] == 1 && myb == 1) AC[i][j] = 0;
            else if (gl3egl3x[i][j] == 1 && ttg1x[i][j] == 1 && myb == 2) {
                if (ttg2x[i][j] == 0 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 0) AC[i][j] = 1;
            } else if (gl3egl3x[i][j] == 1 && ttg1x[i][j] == 1 && myb == 3) {
                if (ttg2x[i][j] == 0 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 1) AC[i][j] = 1;
            } else if (gl3egl3x[i][j] == 2 && ttg1x[i][j] == 1 && myb == 1) {
                if (ttg2x[i][j] == 0 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 1) AC[i][j] = 1;
            } else if (gl3egl3x[i][j] == 2 && ttg1x[i][j] == 1 && myb == 2) {
                if (ttg2x[i][j] == 0 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 1) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 2) AC[i][j] = 1;
            } else if (gl3egl3x[i][j] == 2 && ttg1x[i][j] == 1 && myb == 3) {
                if (ttg2x[i][j] == 0 && IC[i][j] == 0) AC[i][j] = 1;
                else if (ttg2x[i][j] == 0 && IC[i][j] == 1) AC[i][j] = 1;
                else if (ttg2x[i][j] == 0 && IC[i][j] == 2) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 0) AC[i][j] = 2;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 1) AC[i][j] = 1;
                else if (ttg2x[i][j] == 1 && IC[i][j] == 2) AC[i][j] = 1;
            }
           
            if (scmx[i][j] == 0) wery[i][j] = 2;
            else {
                if (scmx[i][j] == 1 && salt == 1) wery[i][j] = 1;
                else wery[i][j] = 0;
            }

            if (AC[i][j] == 0) ttg2y[i][j] = 0;
            else ttg2y[i][j] = 1;

            myb23y[i][j] = 0; 
            if (ite >= 10) { 
                if (AC[i][j] == 1 && ite < 25) myb23y[i][j] = 1; 
                else if (AC[i][j] == 2 && ite < 25) myb23y[i][j] = 1;
                else if (AC[i][j] > 0 && ite >=25) myb23y[i][j] = 1; 
            }

            etc1y[i][j] = 0; 
            if (AC[i][j] == 0 && current_auxx == 2) etc1y[i][j] = 1;
            else if (AC[i][j] == 1 && current_auxx >= 1) etc1y[i][j] = 1;
            else if (AC[i][j] == 2) etc1y[i][j] = 1;

            if (AC[i][j] == 0) {
                if (cpcx[i][j] == 0) gl3egl3Local[i][j] = 1;
                else gl3egl3Local[i][j] = 2; 
            } else if (AC[i][j] == 1) {
                if (cpcx[i][j] == 0) gl3egl3Local[i][j] = 0;
                else gl3egl3Local[i][j] = 1; 
            } else { 
                gl3egl3Local[i][j] = 0;
            }

            if (AC[i][j] == 2 || (AC[i][j] != 0 && (current_auxx == 2 || (current_auxx != 0 && (current_wrky75x == 0 || zfp5x[i][j] == 2))))) {
                cpcLocal[i][j] = 2;
            } else if ((AC[i][j] == 0 && current_auxx == 2) || (AC[i][j] == 1 && (current_auxx == 1 && zfp5x[i][j] < 2 && current_wrky75x < 2))) {
                cpcLocal[i][j] = 1;
            } else {
                cpcLocal[i][j] = 0;
            }

            if (current_etx == 0) zfp5y[i][j] = citoquininas; 
            else if (current_etx == 1) {
                if (citoquininas == 0) zfp5y[i][j] = 1;
                else zfp5y[i][j] = 2; 
            }

            if ((AC[i][j] == 2 && current_auxx != 0) || (AC[i][j] != 0 && current_auxx != 0 && current_wrky75x == 0)) {
                tryy[i][j] = 2;
            } else if (AC[i][j] == 1 && current_auxx == 1 && current_wrky75x == 1) {
                tryy[i][j] = 1;
            } else {
                tryy[i][j] = 0;
            }

            if (phos == 0) *wrky75y_ptr = 2;
            else if (phos == 1) *wrky75y_ptr = 1;
            else *wrky75y_ptr = 0; 

            if (CO2 == 0) {
                if (phos == 0) {
                    if (current_wrky75x == 2) *auxy_ptr = 2; else *auxy_ptr = 1;
                } else if (phos == 1) *auxy_ptr = 1;
                else *auxy_ptr = 0; 
            } else { 
                if (phos == 2) *auxy_ptr = 1; else *auxy_ptr = 2;
            }

            if (phos == 0) *ety_ptr = 1;
            else *ety_ptr = CO2;

            jkdy[i][j] = posicion[i][j]; 

            if (jkdy[i][j] == 1 && IC[i][j] > 0 && AC[i][j] < 2) scmy[i][j] = 1;
            else scmy[i][j] = 0;

            if (AC[i][j] == 0) gl2y[i][j] = 0;
            else gl2y[i][j] = 1; 

            cell[i][j] = (gl2y[i][j] == 0) ? 1 : 0;

            coin = rand() % 10;
            rvar = (coin < deggl3egl3) ? 1 : 0;
            gl3egl3y[i][j] = (gl3egl3Local[i][j] + gl3egl3Dif[i][j] - rvar);
            if (gl3egl3y[i][j] > 2) gl3egl3y[i][j] = 2;
            if (gl3egl3y[i][j] < 0) gl3egl3y[i][j] = 0;

            coin = rand() % 10;
            rvar = (coin < degcpc) ? 1 : 0;
            cpcy[i][j] = (cpcLocal[i][j] + cpcDif[i][j] - rvar);
            if (cpcy[i][j] > 2) cpcy[i][j] = 2;
            if (cpcy[i][j] < 0) cpcy[i][j] = 0;
        } 
    } 

    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            werx[i][j] = wery[i][j];
            myb23x[i][j] = myb23y[i][j];
            gl3egl3x[i][j] = gl3egl3y[i][j];
            cpcx[i][j] = cpcy[i][j];
            tryx[i][j] = tryy[i][j];
            etc1x[i][j] = etc1y[i][j];
            ttg2x[i][j] = ttg2y[i][j];
            gl2x[i][j] = gl2y[i][j];
            scmx[i][j] = scmy[i][j];
            zfp5x[i][j] = zfp5y[i][j];
        }
    }
    *wrky75x_ptr = *wrky75y_ptr;
    *auxx_ptr = *auxy_ptr;
    *etx_ptr = *ety_ptr;

    if (ite == (MAX_ITERATIONS - 1)) {
        int current_Hectopic = 0;
        int current_NHectopic = 0;
        for (i = 0; i < GRID_SIZE; i++) {
            for (j = 0; j < GRID_SIZE; j++) {
                if ((posicion[i][j] == 0) && (cell[i][j] == 1)) {
                    current_Hectopic++;
                }
                if ((posicion[i][j] == 1) && (cell[i][j] == 0)) {
                    current_NHectopic++;
                }
            }
        }
        *Hectopic_ptr = current_Hectopic;
        *NHectopic_ptr = current_NHectopic;
    }
}

void imprimeC(int cell[][GRID_SIZE]) {
    int k, l;
    for (k = 0; k < GRID_SIZE; k++) {
        for (l = 0; l < GRID_SIZE; l++) {
            printf("%d ", cell[k][l]);
        }
        printf("\n"); 
    }
    printf("\n");
}

void imprimedata(int cell[][GRID_SIZE]) {
    FILE *dat = fopen("WT_patron.dat", "w");
    int k, l;
    if (dat) {
        for (k = 0; k < GRID_SIZE; k++) {
            for (l = 0; l < GRID_SIZE; l++) {
                fprintf(dat, "%d ", cell[k][l]);
            }
            fprintf(dat, "\n");
        }
        fprintf(dat, "\n");
        fclose(dat);
    } else {
        perror("Error opening WT_patron.dat for writing");
    }
}

int main() {
    int posicion[GRID_SIZE][GRID_SIZE];
    int werx[GRID_SIZE][GRID_SIZE], wery[GRID_SIZE][GRID_SIZE];
    int myb23x[GRID_SIZE][GRID_SIZE], myb23y[GRID_SIZE][GRID_SIZE];
    int gl3egl3x[GRID_SIZE][GRID_SIZE], gl3egl3y[GRID_SIZE][GRID_SIZE];
    int ttg1x[GRID_SIZE][GRID_SIZE], ttg1y[GRID_SIZE][GRID_SIZE]; 
    int cpcx[GRID_SIZE][GRID_SIZE], cpcy[GRID_SIZE][GRID_SIZE];
    int tryx[GRID_SIZE][GRID_SIZE], tryy[GRID_SIZE][GRID_SIZE];
    int etc1x[GRID_SIZE][GRID_SIZE], etc1y[GRID_SIZE][GRID_SIZE];
    int ttg2x[GRID_SIZE][GRID_SIZE], ttg2y[GRID_SIZE][GRID_SIZE];
    int gl2x[GRID_SIZE][GRID_SIZE], gl2y[GRID_SIZE][GRID_SIZE];
    int zfp5x[GRID_SIZE][GRID_SIZE], zfp5y[GRID_SIZE][GRID_SIZE];
    int scmx[GRID_SIZE][GRID_SIZE], scmy[GRID_SIZE][GRID_SIZE];
    int jkdx[GRID_SIZE][GRID_SIZE], jkdy[GRID_SIZE][GRID_SIZE]; 

    int wrky75x, wrky75y;
    int etx, ety;
    int auxx, auxy;

    int AC[GRID_SIZE][GRID_SIZE], IC[GRID_SIZE][GRID_SIZE];
    float cpcn[GRID_SIZE][GRID_SIZE], gl3egl3n[GRID_SIZE][GRID_SIZE];
    int cpcLocal[GRID_SIZE][GRID_SIZE], gl3egl3Local[GRID_SIZE][GRID_SIZE];
    int cpcDif[GRID_SIZE][GRID_SIZE], gl3egl3Dif[GRID_SIZE][GRID_SIZE];
    
    int cell[GRID_SIZE][GRID_SIZE];

    int ite, contador, Numsim;
    int Hectopic = 0, NHectopic = 0; 
    float sumaH = 0, sumaEH = 0;
    float PromPelosTotal, PromEH;
    
    float *Ptotal_por_simulacion = NULL;
    float *Hectopic_por_simulacion = NULL;

    int total_H_positions = 0; 
    int total_N_positions = 0; 

    int i, j; 

    srand(time(NULL));

    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            if (j % 3 == 0) {
                posicion[i][j] = 1; 
                total_H_positions++;
            } else {
                posicion[i][j] = 0; 
                total_N_positions++;
            }
        }
    }
    
    printf("\n--- Genetic Regulatory Network Simulation ---\n");
    printf("Grid size is %d x %d.\n", GRID_SIZE, GRID_SIZE);
    printf("Predefined H-positions: %d, N-positions: %d, Total cells: %d\n", total_H_positions, total_N_positions, GRID_SIZE*GRID_SIZE);

    printf("\nEnter the parameters for the simulation:\n");
    printf("Diffusion parameter for CPC (e.g., 0.05): ");
    scanf("%f", &DifI);
    printf("Diffusion parameter for EGL3/GL3 (e.g., 0.01): ");
    scanf("%f", &DifA);
    printf("Degradation rate for CPC (integer, e.g., 9): ");
    scanf("%d", &degcpc);
    printf("Degradation rate for EGL3/GL3 (integer, e.g., 1): ");
    scanf("%d", &deggl3egl3);
    printf("Phosphate value (0=Low, 1=Normal, 2=High): ");
    scanf("%d", &phos);
    citoquininas = 1; 
    salt = 0;
    CO2 = 0;

    printf("Enter the total number of simulations to run: ");
    scanf("%d", &Numsim);

    if (Numsim > 0) {
        Ptotal_por_simulacion = (float *)malloc(Numsim * sizeof(float));
        Hectopic_por_simulacion = (float *)malloc(Numsim * sizeof(float));
        if (Ptotal_por_simulacion == NULL || Hectopic_por_simulacion == NULL) {
            perror("Error allocating memory for simulation results");
            // Free any potentially allocated memory before exiting
            if (Ptotal_por_simulacion != NULL) free(Ptotal_por_simulacion);
            if (Hectopic_por_simulacion != NULL) free(Hectopic_por_simulacion);
            return 1; 
        }
    } else {
        printf("Number of simulations must be greater than 0.\n");
        return 1;
    }

    for (contador = 1; contador <= Numsim; contador++) {
        printf("\n--- Running Simulation %d of %d ---\n", contador, Numsim);
        Iniciales(werx, myb23x, gl3egl3x, ttg1x, cpcx, tryx, etc1x, ttg2x, gl2x,
                  zfp5x, scmx, jkdx, &wrky75x, &etx, &auxx,
                  cpcn, gl3egl3n, gl3egl3Local, gl3egl3Dif, cpcLocal, cpcDif);

        for (ite = 0; ite < MAX_ITERATIONS; ite++) {
            Itera(ite, posicion,
                  werx, myb23x, gl3egl3x, ttg1x, cpcx, tryx, etc1x, ttg2x, gl2x,
                  zfp5x, scmx, jkdx, &wrky75x, &etx, &auxx,
                  wery, myb23y, gl3egl3y, ttg1y, cpcy, tryy, etc1y, ttg2y, gl2y,
                  zfp5y, scmy, jkdy, &wrky75y, &ety, &auxy,
                  AC, IC, cpcLocal, gl3egl3Local, cpcDif, gl3egl3Dif,
                  cell, &Hectopic, &NHectopic);
            
            difundecpc(cpcn, cpcx, cpcDif);
            difundegl3egl3(gl3egl3n, gl3egl3x, gl3egl3Dif);
        }

        int Ptotal = (total_H_positions - NHectopic) + Hectopic;

        printf("Results for Simulation %d:\n", contador);
        printf("  Hairs in N-positions (ectopic): %d\n", Hectopic);
        printf("  Non-hairs in H-positions (missing): %d\n", NHectopic);
        printf("  Total hairs observed: %d\n", Ptotal);
        
        sumaH += Ptotal;
        sumaEH += Hectopic; 

        Ptotal_por_simulacion[contador - 1] = (float)Ptotal;
        Hectopic_por_simulacion[contador - 1] = (float)Hectopic;
    }

    PromPelosTotal = (Numsim > 0) ? sumaH / Numsim : 0;
    PromEH = (Numsim > 0) ? sumaEH / Numsim : 0;

    printf("\n--- Overall Averages after %d Simulation(s) ---\n", Numsim);
    printf("Mean total hairs: %.1f\n", PromPelosTotal);
    printf("Mean ectopic hairs (hairs in N-position): %.1f\n", PromEH);

    if (Numsim > 0) {
        float suma_cuadrados_dif_Ptotal = 0;
        float suma_cuadrados_dif_Hectopic = 0;
        int k; 

        for (k = 0; k < Numsim; k++) {
            suma_cuadrados_dif_Ptotal += pow(Ptotal_por_simulacion[k] - PromPelosTotal, 2);
            suma_cuadrados_dif_Hectopic += pow(Hectopic_por_simulacion[k] - PromEH, 2);
        }

        float varianza_Ptotal = suma_cuadrados_dif_Ptotal / Numsim; 
        float desv_estandar_Ptotal = sqrt(varianza_Ptotal);

        float varianza_Hectopic = suma_cuadrados_dif_Hectopic / Numsim; 
        float desv_estandar_Hectopic = sqrt(varianza_Hectopic);

        printf("Standard deviation of total hairs: %.2f\n", desv_estandar_Ptotal);
        printf("Standard deviation of ectopic hairs: %.2f\n", desv_estandar_Hectopic);
    }
    
    printf("\nFinal spatial pattern from the last simulation:\n(1 = hair, 0 = non-hair)\n");
    imprimeC(cell);
    imprimedata(cell);

    if (Ptotal_por_simulacion != NULL) {
        free(Ptotal_por_simulacion);
    }
    if (Hectopic_por_simulacion != NULL) {
        free(Hectopic_por_simulacion);
    }

    return 0;
}
