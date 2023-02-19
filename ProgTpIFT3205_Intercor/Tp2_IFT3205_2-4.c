/*------------------------------------------------------*/
/* Prog    : Tp2_IFT3205-2-4.c                          */
/* Auteur  :  Loïc Daudé Mondet (20243814)  Adel Abdeladim (20127626) */
/* Emails  :  loic.daude.mondet@umontreal.ca   adel.abdeladim@umontreal.ca  */
/* Date    : --/--/2010                                 */
/* version :                                            */
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo2.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define NAME_VISUALISER "display "
#define NAME_IMG_IN1  "UdM_1"
#define NAME_IMG_IN2  "UdM_2"
#define NAME_IMG_OUT "image-TpIFT3205-2-4"

/*------------------------------------------------*/
/* PROTOTYPE DE FONCTIONS  -----------------------*/
/*------------------------------------------------*/
void CenterImg2(float*** MatriceImgM,int length, int width){
    float ** CenteredMatrix = fmatrix_allocate_2d(length, width);
    int center_length = length%2==0 ? length/2-1 : length/2;
    int center_width = width%2==0 ? width/2-1 : width/2;
    for (int i=0; i<width; i++) {
        for (int j=0; j<length; j++) {
            if (i<=center_width) {
                if (j<=center_length) {
                    CenteredMatrix[i+center_width+1][j+center_length+1]=(*MatriceImgM)[i][j];
                } else {
                    CenteredMatrix[i+center_width+1][j-center_length-1]=(*MatriceImgM)[i][j];
                }
            } else {
                if (j<=center_length) {
                    CenteredMatrix[i-center_width-1][j+center_length+1]=(*MatriceImgM)[i][j];
                } else {
                    CenteredMatrix[i-center_width-1][j-center_length-1]=(*MatriceImgM)[i][j];
                }
            }
        }
    }

    free_fmatrix_2d(*MatriceImgM);
    *MatriceImgM=CenteredMatrix;
}

void applyLog(float ** MatriceImgM, int length, int width) {
    for (int i=0; i<width; i++) {
        for (int j=0; j<length; j++) {
            MatriceImgM[i][j]=logf(1+MatriceImgM[i][j]);
        }
    }
}

void rotateImg(float ** MatriceImg, float ** MatriceImgRotate, int length, int width, float Theta) {
    int i,j,k;
    int center_length = length%2==0 ? length/2-1 : length/2;
    int center_width = width%2==0 ? width/2-1 : width/2;
    for(i=0;i<length;i++)
        for(j=0;j<width;j++)
        {

            int x_r = roundf((i-center_length) * cosf(Theta) + (j-center_width) * sinf(Theta) + center_length);
            int y_r = roundf(-(i-center_length) * sinf(Theta) + (j-center_width) * cosf(Theta) + center_width);
            if (y_r<0 || y_r>=length || x_r<0 || x_r>=length) MatriceImgRotate[i][j]=0;
            else MatriceImgRotate[i][j]=MatriceImg[x_r][y_r];
        }

}

void rotateImgInterp(float ** MatriceImg, float ** MatriceImgRotate, int length, int width, float Theta) {
    int i,j,k;
    int center_length = length%2==0 ? length/2-1 : length/2;
    int center_width = width%2==0 ? width/2-1 : width/2;
    for(i=0;i<length;i++)
        for(j=0;j<width;j++)
        {
            float x_p = (i-center_length) * cosf(Theta) + (j-center_width) * sinf(Theta) + center_length;
            float y_p = -(i-center_length) * sinf(Theta) + (j-center_width) * cosf(Theta) + center_width;
            int x = roundf(x_p);
            x = x>x_p ? x-1 : x;
            int y = roundf(y_p);
            y = y>y_p ? y-1 : y;
            float f_xy = x<0 || y<0 || x>=length || y>= width ? 0 : MatriceImg[x][y];
            float f_x1y = x+1<0 || y<0 || x+1>=length || y>= width ? 0 : MatriceImg[x+1][y];
            float f_xy1 = x<0 || y+1<0 || x>=length || y+1>= width ? 0 : MatriceImg[x][y+1];
            float f_x1y1 = x+1<0 || y+1<0 || x+1>=length || y+1>= width ? 0 : MatriceImg[x+1][y+1];
            float f_xpy = f_xy + (x_p-x)*(f_x1y-f_xy);
            float f_xpy1 = f_xy1 + (x_p-x)*(f_x1y1-f_xy1);
            float f_xpyp = f_xpy + (y_p-y)*(f_xpy1-f_xpy);
            MatriceImgRotate[i][j]=f_xpyp;
        }

}

void copyMatrice(float *** dest, float ** src, int length, int width) {
    for (int i=0; i<length; i++) {
        for (int j=0; j<width; j++) {
            (*dest)[i][j]=src[i][j];
        }
    }
}

void printMatrice(float ** mat, int length, int width) {
    for (int i=0; i<length; i++) {
        for (int j=0; j<width; j++) {
            printf("%f", mat[i][j]);
        }
        printf("\n");
    }
}

float findTheta(float ** MatriceImgF, float ** MatriceImgG, int length, int width) {
    float** MatriceImgIF=fmatrix_allocate_2d(length, width);
    float** MatriceImgRF=fmatrix_allocate_2d(length, width);
    float** MatriceImgMF=fmatrix_allocate_2d(length, width);

    float** MatriceImgCG=fmatrix_allocate_2d(length, width);

    float** MatriceImgIGRotate=fmatrix_allocate_2d(length, width);
    float** MatriceImgRGRotate=fmatrix_allocate_2d(length, width);
    float** MatriceImgMGRotate=fmatrix_allocate_2d(length, width);


    for(int i=0;i<length;i++)
        for(int j=0;j<width;j++)
        {
            MatriceImgIF[i][j]=0.0;
            MatriceImgRF[i][j]=0.0;
            MatriceImgMF[i][j]=0.0;

            MatriceImgCG[i][j]=0.0;


            MatriceImgIGRotate[i][j]=0.0;
            MatriceImgRGRotate[i][j]=0.0;
            MatriceImgMGRotate[i][j]=0.0;
        }


    copyMatrice(&MatriceImgRF, MatriceImgF, length, width);
    copyMatrice(&MatriceImgCG, MatriceImgG, length, width);



    FFTDD(MatriceImgRF, MatriceImgIF, length, width);
    Mod(MatriceImgMF, MatriceImgRF, MatriceImgIF, length, width);
    applyLog(MatriceImgMF, length, width);
    Recal(MatriceImgMF,length,width);


    float start = -(PI/16);
    float stop = PI/16;
    float step = 0.005;
    double minError=0;
    float minTheta=start;
    double modulesDiff;

    rotateImgInterp(MatriceImgCG, MatriceImgRGRotate, length, width, start);
    FFTDD(MatriceImgRGRotate, MatriceImgIGRotate, length, width);
    Mod(MatriceImgMGRotate, MatriceImgRGRotate, MatriceImgIGRotate, length, width);
    applyLog(MatriceImgMGRotate, length, width);
    Recal(MatriceImgMGRotate,length,width);

    for (int i=0; i<length; i++) {
        for (int j=0; j<width; j++) {
            modulesDiff= MatriceImgMGRotate[i][j] - MatriceImgMF[i][j];
            modulesDiff = modulesDiff<0 ? -modulesDiff : modulesDiff;
            minError = minError +  modulesDiff;
        }
    }
    printf("[%f::%f]",start, minError);
    start = start + step;
    double currentError;
    while (start<stop) {

        for(int i=0;i<length;i++)
            for(int j=0;j<width;j++)
            {
                MatriceImgIGRotate[i][j]=0.0;
                MatriceImgRGRotate[i][j]=0.0;
                MatriceImgMGRotate[i][j]=0.0;
            }


        currentError=0;
        rotateImgInterp(MatriceImgCG, MatriceImgRGRotate, length, width, start);
        FFTDD(MatriceImgRGRotate, MatriceImgIGRotate, length, width);
        Mod(MatriceImgMGRotate, MatriceImgRGRotate, MatriceImgIGRotate, length, width);
        applyLog(MatriceImgMGRotate, length, width);
        Recal(MatriceImgMGRotate,length,width);
        for (int i=0; i<length; i++) {
            for (int j=0; j<width; j++) {
                modulesDiff = MatriceImgMGRotate[i][j] - MatriceImgMF[i][j];
                modulesDiff = modulesDiff < 0 ? -modulesDiff : modulesDiff;
                currentError = currentError +  modulesDiff;
            }
        }
        if (currentError<minError) {
            minError=currentError;
            minTheta=start;
        }
        printf(">[%f::%f]",start, currentError);
        start = start + step;
    }
    printf("\n");

    return minTheta;
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/
/*------------------------------------------------*/
int main(int argc,char **argv)
{
    int i,j,k;
    int length,width;
    float Theta0;
    int x0,y0;
    char BufSystVisuImg[100];

    //Constante
    length=512;
    width=512;

    //Allocation Memoire
    float** MatriceImgI1=fmatrix_allocate_2d(length,width);
    float** MatriceImgM1=fmatrix_allocate_2d(length,width);
    float** MatriceImgR1=fmatrix_allocate_2d(length,width);

    float** MatriceImgI2=fmatrix_allocate_2d(length,width);
    float** MatriceImgM2=fmatrix_allocate_2d(length,width);
    float** MatriceImgR2=fmatrix_allocate_2d(length,width);

    float** MatriceImgI3=fmatrix_allocate_2d(length,width);
    float** MatriceImgM3=fmatrix_allocate_2d(length,width);
    float** MatriceImgR3=fmatrix_allocate_2d(length,width);

    float** MatriceImg3=fmatrix_allocate_2d(length,width);

    //Lecture Image
    float** MatriceImg1=LoadImagePgm(NAME_IMG_IN1,&length,&width);
    float** MatriceImg2=LoadImagePgm(NAME_IMG_IN2,&length,&width);


    // .... .... .... .... .... .... ....
    float theta = findTheta(MatriceImg2, MatriceImg1, length, width);
    printf("Angle >>> [%f]\n", theta);
    float** MatriceImgGRotate=fmatrix_allocate_2d(length, width);
    rotateImgInterp(MatriceImg2, MatriceImgGRotate, length, width, -theta);

    //Sauvegarde
    SaveImagePgm(NAME_IMG_OUT,MatriceImgGRotate,length,width);

    //Commande systeme: VISU
    strcpy(BufSystVisuImg,NAME_VISUALISER);
    strcat(BufSystVisuImg,NAME_IMG_OUT);
    strcat(BufSystVisuImg,".pgm&");
    printf(" %s",BufSystVisuImg);
    system(BufSystVisuImg);
    strcpy(BufSystVisuImg,NAME_VISUALISER);
    strcat(BufSystVisuImg,NAME_IMG_IN1);
    strcat(BufSystVisuImg,".pgm&");
    printf(" %s",BufSystVisuImg);
    system(BufSystVisuImg);

    strcpy(BufSystVisuImg,NAME_VISUALISER);
    strcat(BufSystVisuImg,NAME_IMG_IN2);
    strcat(BufSystVisuImg,".pgm&");
    printf(" %s",BufSystVisuImg);
    system(BufSystVisuImg);


    //==End=========================================================

    //Liberation memoire
    free_fmatrix_2d(MatriceImgR1);
    free_fmatrix_2d(MatriceImgI1);
    free_fmatrix_2d(MatriceImgM1);
    free_fmatrix_2d(MatriceImgR2);
    free_fmatrix_2d(MatriceImgI2);
    free_fmatrix_2d(MatriceImgM2);
    free_fmatrix_2d(MatriceImgR3);
    free_fmatrix_2d(MatriceImgI3);
    free_fmatrix_2d(MatriceImgM3);
    free_fmatrix_2d(MatriceImg1);
    //free_fmatrix_2d(MatriceImg2);
    free_fmatrix_2d(MatriceImg3);

    //retour sans probleme
    printf("\n C'est fini ... \n\n");
    return 0;
}


