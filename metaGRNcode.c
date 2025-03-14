/* This program is an update to the simulation of the coupling of genetic regulatory networks in a two-dimensional domain.This domain is composed of a web of cells 
representing epidermal cells of the root of Arabidopsis thaliana.
Each cell contains a network composed of the proteins wer, myb23, gl3/egl3, ttg1, cpc, try, etc1, ttg2, gl2, zfp5, wrky75, scm;
the hormones ethylene, aux, and cytokinins; and the activator AC and inhibitor IC protein complexes (the latter does not occur 
in nature; it is included to facilitate the program). The cells communicate by simulating cpc and gl3/egl3 diffusion. */

	/*
The parameters are:


	DifA=0.01; // diff of gl3/egl3
	DifI=0.05; // diff of cpc  
  degcpc = 9 // degradation of cpc
  deggl3egl3= 1 // degradatio of gl3/egl3

*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>


int i,j,k,l,w,u,v,m,n,d,f,coin,rvar;
int NHectopic,Hectopic,pelos, Ptotal, NPtotal,Numsim;
float sumaH=0;
float sumaEH=0;
float PromPelosTotal; //
float PromEH;
int contador;
int ite;
int veces=85; 
int citoquininas;
int phos;
int salt;
int CO2;

float vector[1000],num,suma,promedio, sigma,desv_standar;



int posicion[24][24]; 
int werx[24][24];
int myb23x[24][24];
int gl3egl3x[24][24];
int ttg1x[24][24];
int cpcx[24][24];
int tryx[24][24];
int etc1x[24][24];
int ttg2x[24][24];
int gl2x[24][24];
int zfp5x[24][24];
int scmx[24][24];
int wrky75x;
int jkdx[24][24];


int etx;
int auxx;

int wery[24][24];
int myb23y[24][24];
int gl3egl3y[24][24]; 
int ttg1y[24][24];
int cpcy[24][24]; 
int tryy[24][24];
int etc1y[24][24];
int ttg2y[24][24];
int gl2y[24][24];
int zfp5y[24][24];
int scmy[24][24];
int wrky75y;
int jkdy[24][24];

int ety;
int auxy;

int myb; 

int AC[24][24];
int IC[24][24];

float cpcn[24][24]; 
float gl3egl3n[24][24];


int cpcLocal[24][24];
int gl3egl3Local[24][24];

int cpcDif[24][24];
int gl3egl3Dif[24][24];

int cel=24; 
float DifI;// 
float DifA;// 


int degcpc;
int deggl3egl3;

int cell[24][24];


int Iniciales(){
   
	for(i=0;i<cel;i++){
		for(j=0;j<cel;j++){
			werx[i][j]=rand()%2;
			myb23x[i][j]=0;
			gl3egl3x[i][j]=rand()%2;
			gl3egl3Local[i][j]=rand()%2;
			gl3egl3Dif[i][j]=rand()%2;
			ttg1x[i][j]=1; //Lo fijamos
			cpcx[i][j]=rand()%2;
			cpcLocal[i][j]=rand()%2;
			cpcDif[i][j]=rand()%2;
			tryx[i][j]=rand()%2;
			etc1x[i][j]=rand()%2;
			ttg2x[i][j]=rand()%2;
			gl2x[i][j]=rand()%2;
			zfp5x[i][j]=rand()%2;
			scmx[i][j]=rand()%2;
			jkdx[i][j]=rand()%2;
			wrky75x=1;
			etx=0;
			auxx=1; 
			cpcn[i][j]=0; 
			gl3egl3n[i][j]=0; 


		}} //Cierras for's
} //Cierra iniciales


int H(float x){ 
	int y;
	if(x<=0.25)y=0;
	if(0.25<x&&x<=1)y=1;
	if(x>1)y=2;
	return y;
}


difundecpc(){ 
	for(i=0;i<cel;i++){
		for(j=0;j<cel;j++){

		

			if(i==0){
				if(j==0){
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i+1][j]+cpcx[i][j+1]+cpcx[i][cel-1]-3*cpcx[i][j]);
				}
				else if(j==(cel-1)){
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i][j-1]+cpcx[i+1][j]+cpcx[i][0]-3*cpcx[i][j]);
				}
				else{
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i+1][j]+cpcx[i][j+1]+cpcx[i][j-1]-3*cpcx[i][j]);
				}
			}

			else if(i==(cel-1)){
				if(j==0){
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i-1][j]+cpcx[i][j+1]+cpcx[i][cel-1]-3*cpcx[i][j]);
				}
				else if(j==(cel-1)){
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i-1][j]+cpcx[i][j-1]+cpcx[i][0]-3*cpcx[i][j]);
				}
				else{
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i-1][j]+cpcx[i][j+1]+cpcx[i][j-1]-3*cpcx[i][j]);
				}
			}

			else if(j==0&&i!=(cel-1)&&i!=0){
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i-1][j]+cpcx[i][j+1]+cpcx[i+1][j]+cpcx[i][cel-1]-4*cpcx[i][j]);
			}

			else if(j==(cel-1)&&i!=(cel-1)&&i!=0){
				cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i-1][j]+cpcx[i+1][j]+cpcx[i][j-1]+cpcx[i][0]-4*cpcx[i][j]);

			else{
			cpcn[i][j]=cpcx[i][j]+DifI*(cpcx[i+1][j]+cpcx[i-1][j]+cpcx[i][j+1]+cpcx[i][j-1]-4*cpcx[i][j]);
			}


			cpcDif[i][j]=H(cpcn[i][j]);

		}//
	}//

}//¿

difundegl3egl3(){)
	for(i=0;i<cel;i++){
		for(j=0;j<cel;j++){


				if(i==0){
					if(j==0){
						gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i+1][j]+gl3egl3x[i][j+1]+gl3egl3x[i][cel-1]-3*gl3egl3x[i][j]);
					}
					else if(j==(cel-1)){
						gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i][j-1]+gl3egl3x[i+1][j]+gl3egl3x[i][0]-3*gl3egl3x[i][j]);
					}
					else{
						gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i+1][j]+gl3egl3x[i][j+1]+gl3egl3x[i][j-1]-3*gl3egl3x[i][j]);
					}
				}

				else if(i==(cel-1)){
					if(j==0){
						gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i-1][j]+gl3egl3x[i][j+1]+gl3egl3x[i][cel-1]-3*gl3egl3x[i][j]);
					}
					else if(j==(cel-1)){
						gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i-1][j]+gl3egl3x[i][j-1]+gl3egl3x[i][0]-3*gl3egl3x[i][j]);
					}
					else{
						gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i-1][j]+gl3egl3x[i][j+1]+gl3egl3x[i][j-1]-3*gl3egl3x[i][j]);
					}
				}

				else if(j==0&&i!=(cel-1)&&i!=0){
					gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i-1][j]+gl3egl3x[i][j+1]+gl3egl3x[i+1][j]+gl3egl3x[i][cel-1]-4*gl3egl3x[i][j]);
				}

				else if(j==(cel-1)&&i!=(cel-1)&&i!=0){
					gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i-1][j]+gl3egl3x[i+1][j]+gl3egl3x[i][j-1]+gl3egl3x[i][0]-4*gl3egl3x[i][j]);
				}

				else{
					gl3egl3n[i][j]=gl3egl3x[i][j]+DifA*(gl3egl3x[i+1][j]+gl3egl3x[i-1][j]+gl3egl3x[i][j+1]+gl3egl3x[i][j-1]-4*gl3egl3x[i][j]);
				}


			gl3egl3Dif[i][j]=H(gl3egl3n[i][j]);

		}//
	}//

}//

int Itera() { 

	//Para cada c�lula...
	for(i=0;i<cel;i++){
		for(j=0;j<cel;j++){

			//IC
			if(cpcx[i][j]==0&&tryx[i][j]==0&&etc1x[i][j]==0) IC[i][j]=0;
			if(cpcx[i][j]==0&&tryx[i][j]==0&&etc1x[i][j]==1) IC[i][j]=0;
			if(cpcx[i][j]==0&&tryx[i][j]==1&&etc1x[i][j]==0) IC[i][j]=0;
			if(cpcx[i][j]==0&&tryx[i][j]==1&&etc1x[i][j]==1) IC[i][j]=1;
			if(cpcx[i][j]==0&&tryx[i][j]==2&&etc1x[i][j]==0) IC[i][j]=1;
			if(cpcx[i][j]==0&&tryx[i][j]==2&&etc1x[i][j]==1) IC[i][j]=2;
			if(cpcx[i][j]==1&&tryx[i][j]==0&&etc1x[i][j]==0) IC[i][j]=1;
			if(cpcx[i][j]==1&&tryx[i][j]==0&&etc1x[i][j]==1) IC[i][j]=1;
			if(cpcx[i][j]==1&&tryx[i][j]==1&&etc1x[i][j]==0) IC[i][j]=1;
			if(cpcx[i][j]==1&&tryx[i][j]==1&&etc1x[i][j]==1) IC[i][j]=1;
			if(cpcx[i][j]==1&&tryx[i][j]==2&&etc1x[i][j]==0) IC[i][j]=2;
			if(cpcx[i][j]==1&&tryx[i][j]==2&&etc1x[i][j]==1) IC[i][j]=2;
			if(cpcx[i][j]==2&&tryx[i][j]==0&&etc1x[i][j]==0) IC[i][j]=1;
			if(cpcx[i][j]==2&&tryx[i][j]==0&&etc1x[i][j]==1) IC[i][j]=2;
			if(cpcx[i][j]==2&&tryx[i][j]==1&&etc1x[i][j]==0) IC[i][j]=2;
			if(cpcx[i][j]==2&&tryx[i][j]==1&&etc1x[i][j]==1) IC[i][j]=2;
			if(cpcx[i][j]==2&&tryx[i][j]==2&&etc1x[i][j]==0) IC[i][j]=2;
			if(cpcx[i][j]==2&&tryx[i][j]==2&&etc1x[i][j]==1) IC[i][j]=2;


      myb=werx[i][j]+myb23x[i][j];

			//MBW

			if(gl3egl3x[i][j]==0)AC[i][j]=0;
			if(myb==0)AC[i][j]=0;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==0)AC[i][j]=0;

			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==1)AC[i][j]=0;

			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==0&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==0&&IC[i][j]==1)AC[i][j]=0;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==0&&IC[i][j]==2)AC[i][j]=0;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==1&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==1&&IC[i][j]==1)AC[i][j]=0;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==1&&IC[i][j]==2)AC[i][j]=0;

			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==0&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==0&&IC[i][j]==1)AC[i][j]=0;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==0&&IC[i][j]==2)AC[i][j]=0;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==1&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==1&&IC[i][j]==1)AC[i][j]=1;
			if(gl3egl3x[i][j]==1&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==1&&IC[i][j]==2)AC[i][j]=0;

			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==1&&ttg2x[i][j]==0&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==1&&ttg2x[i][j]==0&&IC[i][j]==1)AC[i][j]=0;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==1&&ttg2x[i][j]==0&&IC[i][j]==2)AC[i][j]=0;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==1&&ttg2x[i][j]==1&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==1&&ttg2x[i][j]==1&&IC[i][j]==1)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==1&&ttg2x[i][j]==1&&IC[i][j]==2)AC[i][j]=0;

			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==0&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==0&&IC[i][j]==1)AC[i][j]=0;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==0&&IC[i][j]==2)AC[i][j]=0;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==1&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==1&&IC[i][j]==1)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==2&&ttg2x[i][j]==1&&IC[i][j]==2)AC[i][j]=1;

			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==0&&IC[i][j]==0)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==0&&IC[i][j]==1)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==0&&IC[i][j]==2)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==1&&IC[i][j]==0)AC[i][j]=2;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==1&&IC[i][j]==1)AC[i][j]=1;
			if(gl3egl3x[i][j]==2&&ttg1x[i][j]==1&&myb==3&&ttg2x[i][j]==1&&IC[i][j]==2)AC[i][j]=1;

			//WER.
			if (scmx[i][j]==0)wery[i][j]=2;
			else {
				if(scmx[i][j]==1&&salt==1)wery[i][j]=1;
				else wery[i][j]=0;
			}

			//TTG2
			if (AC[i][j]==0)ttg2y[i][j]=0;
			if (AC[i][j]==1)ttg2y[i][j]=1;
			if (AC[i][j]==2)ttg2y[i][j]=1;

      //MYB23
      if (ite<10)myb23y[i][j]=0;
			if (AC[i][j]==0&&(!(ite<25)))myb23y[i][j]=0;
			if (AC[i][j]==1&&(!(ite<25)))myb23y[i][j]=1;
			if (AC[i][j]==2&&(!(ite<25)))myb23y[i][j]=1;

			//ETC1
			if (AC[i][j]==0&&auxx==0)etc1y[i][j]=0;
			if (AC[i][j]==0&&auxx==1)etc1y[i][j]=0;
			if (AC[i][j]==0&&auxx==2)etc1y[i][j]=1;
			if (AC[i][j]==1&&auxx==0)etc1y[i][j]=0;
			if (AC[i][j]==1&&auxx==1)etc1y[i][j]=1;
			if (AC[i][j]==1&&auxx==2)etc1y[i][j]=1;
			if (AC[i][j]==2&&auxx==0)etc1y[i][j]=1;
			if (AC[i][j]==2&&auxx==1)etc1y[i][j]=1;
			if (AC[i][j]==2&&auxx==2)etc1y[i][j]=1;

			//GL3-EGL3
			if (AC[i][j]==0&&cpcx[i][j]==0)gl3egl3Local[i][j]=1;
			if (AC[i][j]==0&&cpcx[i][j]==1)gl3egl3Local[i][j]=2;
			if (AC[i][j]==0&&cpcx[i][j]==2)gl3egl3Local[i][j]=2;
			if (AC[i][j]==1&&cpcx[i][j]==0)gl3egl3Local[i][j]=0;
			if (AC[i][j]==1&&cpcx[i][j]==1)gl3egl3Local[i][j]=1;
			if (AC[i][j]==1&&cpcx[i][j]==2)gl3egl3Local[i][j]=1;
			if (AC[i][j]==2)gl3egl3Local[i][j]=0;

			//CPC
			if (AC[i][j]==2||(AC[i][j]!=0&&(auxx==2||(auxx!=0&&(wrky75x==0||zfp5x[i][j]==2))))) cpcLocal[i][j]=2;
			else {
				if ((AC[i][j]==0&&auxx==2)||(AC[i][j]==1&&(auxx==1&&zfp5x[i][j]<2&&wrky75x<2))) cpcLocal[i][j]=1;
				else cpcLocal[i][j]=0;
			}


			//ZFP5
			if (etx==0&&citoquininas==0)zfp5y[i][j]=0;
			if (etx==0&&citoquininas==1)zfp5y[i][j]=1;
			if (etx==0&&citoquininas==2)zfp5y[i][j]=2;
			if (etx==1&&citoquininas==0)zfp5y[i][j]=1;
			if (etx==1&&citoquininas==1)zfp5y[i][j]=2;
			if (etx==1&&citoquininas==2)zfp5y[i][j]=2;

			//TRY
			if ((AC[i][j]==2&&auxx!=0)||(AC[i][j]!=0&&auxx!=0&&wrky75x==0)) tryy[i][j]=2;
			else {
				if (AC[i][j]==1&&auxx==1&&wrky75x==1) tryy[i][j]=1;
				else tryy[i][j]=0;
			}


			//WRKY75
			if (phos==0)wrky75y=2;
			if (phos==1)wrky75y=1;
			if (phos==2)wrky75y=0;

			//AUX
			if (CO2==0&&phos==0&&wrky75x==0)auxy=1;
			if (CO2==0&&phos==0&&wrky75x==1)auxy=1;
			if (CO2==0&&phos==0&&wrky75x==2)auxy=2;
			if (CO2==0&&phos==1)auxy=1;
			if (CO2==0&&phos==2)auxy=0;
			if (CO2==1&&phos==0)auxy=2;
			if (CO2==1&&phos==1)auxy=2;
			if (CO2==1&&phos==2)auxy=1;



			//Ethylene
			if (phos == 0) ety=1;
			else ety=CO2;
			// reglas jkd
			if (posicion[i][j]==1) jkdy[i][j]=1;
			else jkdy[i][j]=0;

			//SCM
			if (jkdy[i][j]==1 && IC[i][j]>0 && AC[i][j]<2) scmy[i][j]=1;
			else scmy[i][j]=0;

			//GL2
			if(AC[i][j]==0)gl2y[i][j]=0;
			if(AC[i][j]>0)gl2y[i][j]=1;

			// Tricoblasto = 1, Atricoblasto = 0
			if(gl2y[i][j]==0)cell[i][j]=1; 
			if(gl2y[i][j]==1)cell[i][j]=0; 


			coin=rand()%10;
			if(coin<deggl3egl3)rvar=1;
			else rvar=0;
			gl3egl3y[i][j]=(gl3egl3Local[i][j]+gl3egl3Dif[i][j]-rvar);
			if(gl3egl3y[i][j]>2)gl3egl3y[i][j]=2;
			if(gl3egl3y[i][j]<0)gl3egl3y[i][j]=0;

			coin=rand()%10;
			if(coin<degcpc)rvar=1;
			else rvar=0;
			cpcy[i][j]=(cpcLocal[i][j]+cpcDif[i][j]-rvar);
			if(cpcy[i][j]>2)cpcy[i][j]=2;
			if(cpcy[i][j]<0)cpcy[i][j]=0;

		}//cierra j
	} //cierra i

	//Primero aplicas las reglas a toda la matriz y defines las c�lulas con base en eso. Despu�s cambias x por y
	for(i=0;i<cel;i++){//aqui el vector actualizado se convierte en la semilla de la sig. iteracion
		for(j=0;j<cel;j++){
			werx[i][j]=wery[i][j];
			myb23x[i][j]=myb23y[i][j];
			gl3egl3x[i][j]=gl3egl3y[i][j];
			cpcx[i][j]=cpcy[i][j];
			tryx[i][j]=tryy[i][j];
			etc1x[i][j]=etc1y[i][j];
			ttg2x[i][j]=ttg2y[i][j];
			gl2x[i][j]=gl2y[i][j];
			scmx[i][j]=scmy[i][j];
			wrky75x=wrky75y;
			auxx=auxy;
                        etx=ety;
			zfp5x[i][j]=zfp5y[i][j];
		}
	}

	if(ite==(veces-1)){
	Hectopic=0;
		for(i=0;i<cel;i++){
			for(j=0;j<cel;j++){
				if((j%(3)!=0)&&(cell[i][j]==1)) Hectopic++;
			}
		}
	}

	if(ite==(veces-1)){
	NHectopic=0;
		for(i=0;i<cel;i++){
			for(j=0;j<cel;j++){
				if((j%(3)==0)&&(cell[i][j]==0)) NHectopic++;
			}
		}
	}
} //Cierra itera

void imprimeC(){
	for (k=0;k<cel;k++){
		for(l=0;l<cel;l++){
			printf("%d ", cell[k][l]);
			}puts("");
	}puts("");

}


void imprimedata(){
FILE *dat=fopen("WT_patron.dat","w");
if (dat)
{
    	for (k=0;k<cel;k++){
                for(l=0;l<cel;l++){
                fprintf(dat,"%d ", cell[k][l]);
                }
                fprintf(dat,"\n");
    	}
    	fprintf(dat,"\n");
}
}

void imprimeVector(){
	for (k=0;k<cel;k++){
		for(l=0;l<cel;l++){
			printf("%d %d %d %d %d %d %d %d %d ||", AC[k][l],IC[k][l],werx[k][l],myb23x[k][l],gl2x[k][l],gl3egl3x[k][l],ttg1x[k][l],cpcx[k][l],tryx[k][l],etc1x[k][l],ttg2x[k][l],scmx[k][l],zfp5x[k][l]);
		}puts("");
	}puts("");

}



main() {
srand(time(NULL)); 

	for(i=0;i<cel;i++){ 
		for(j=0;j<cel;j++){
			if(j%(3)==0) posicion[i][j]=1;
			else posicion[i][j]=0;
		}
	}
	 system("color 70");

	
	citoquininas=1;
	salt=0;
	CO2=0;
	num=0,suma=0,promedio=0,sigma=0,desv_standar=0;i=0;


    printf("\n Enter the parameters for the simulation of the continuous variables of the network coupling.\n");
    printf("\n The diffusion parameter for CPC is: ");
    scanf("%f",&DifI); 
    printf("\n The diffusion parameter for EGL3/GL3 is: ");
    scanf("%f",&DifA);
    printf("\n The degradation parameter for CPC is:");
    scanf("%d",&degcpc); 
    printf("\n The degradation parameter for EGL3/GL3 is:");
    scanf("%d",&deggl3egl3);
    printf("\n Enter the phosphate value: 0)Low 1)normal, 2)High \n");
    scanf("%d", &phos); 
    printf("\n Enter the total number of simulations \n ");
	scanf("%d", &Numsim);




for (contador=1; contador<=Numsim; contador++){

	Iniciales();

		for (ite=0;ite<veces;ite++){
            
			Itera();

		
			difundecpc();
			difundegl3egl3();

			//while(i<Numsim){
        //for (i=0;i<Numsim;i++) vector[i]=Hectopic;}
		}

        Ptotal=Hectopic+(140-NHectopic);
        NPtotal= 360-NHectopic;
        printf("\n The total hair number in N-position is: ");
        printf("%d", Hectopic); 
		printf(" \n The total hair number is: ");
		printf("%d\n",Ptotal);
        sumaH=sumaH+Ptotal;
        sumaEH=sumaEH+Hectopic;
		 printf("\n The total hair number in H-position is:  ");
		printf("%f \n",sumaH);
		 printf("\n The total hair number in N-position is:  ");
		printf("%f \n",sumaEH);



}
PromPelosTotal=sumaH/Numsim;
PromEH=sumaEH/Numsim;

        printf("\n The mean total hairs is: \n");
        printf("%.1f",PromPelosTotal);
        printf("\n The mean of hairs in hair-position is: \n");
        printf("%.1f",PromEH);
        printf("\n This is the spatial pattern for the parameters defined above, values ​​of 1 are hairs and 0 are non-hairs. \n");
        imprimeC(); 
        //printf("\n The attractors are: \n");
        //imprimeVector();
        puts("");
        imprimedata();
        printf(vector);

        //imprimeVector();
}//main
