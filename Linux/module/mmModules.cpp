//14 Октября 2010. программа для поиска групп сайтов в окне определяемых моделью, использует найденные сайты mm
//основана на MatrixCatchTools,   все максимально вынесено в .h
//Доработана для Продорик (MaxMatrLength=35->45)

//#define Win32
#define Linux

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>

#ifdef Win32
// For Windows
#include <io.h>
#include <sys\stat.h>
#else
// For linux
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#endif


//from mmtools.cpp
#define outusr stdout
#define outerr stderr
#define MaxLineLength 200       // in any file
#define MaxSeqLength 2000000000   // Max sequence lenth
#define MaxMatrLength 45 // Max length of matrix
#define MaxMatrixNum 2000 // Max number of matrices
#define MaxMatrixNameLen 25 //
#define MaxNumofSeq 50000 // максимальное кол-во сиквенсов в файле
#define inittime {T0=time(NULL);}
#define printtime {T1=time(NULL);printf ("Time - %li ",T1-T0);T0=T1;}
#define printtimenL {T1=time(NULL);printf ("Time - %li\n",T1-T0);T0=T1;}
#define printtimeOutErrL {T1=time(NULL);fprintf (outerr,"Time - %li\n",T1-T0);T0=T1;}

//from cefinder.cpp
#define MaxEVNum 1500
#define MaxNumCompElem 50000
#define MaxFNum 6000
#define MaxSiteLength 200 // длина одного из сайтов КЭ, c краями +-30!
#define SeqBorders 10

//for Thresholds and ...
# define PwmThrNum 21  // количество значений порогов. Сами значения должны быть равномерно распределены, т.к. попадание к класс вычисляется, а не сравнивается
# define LenThrNum 1  //  вариация по длине
# define CSThrNum 1      //    по композ скору


//#define MCatchOutPut

#define PrintFirstOnly

#include "MMtypes.h"

// перменные для типов, сиквенсов, матриц...
enum formats {EMBLf,FASTAf,PlainSeqf,UnDef} format=UnDef;
int SeqLength,xl; // если unsigned  то когда Seq[-5] то это уже не -5!
int *SeqList[MaxNumofSeq],*Seq; // макс число сикв ~ 2*MaxNumCompElem
PentaP* PentaSeq, P;

char *YesSitesF,*NoSitesF;
MatrixElemP MatrixList [MaxMatrixNum], M;
int ***SitesDistrY, ***SitesDistrN,*MatrixIndex;

ModuleElemP CEMList[MaxNumCompElem],CE,CEMRList[MaxNumCompElem],CER;

//технические переменные, для данного метода
int NumberOfMatrix;
int NumOfSites=0,fast=0,all=0,TotalSitesLen=0,NumofSeq=0;
int TotalNumOfSitesY=0, TotalNumOfSitesN=0, TotalSeqLengthY=0,TotalSeqLengthN=0,CurrentNumOfSites=0;
InfBlockP InfoY,InfoN;
SiteElemP **SitesListY,**SitesListN, *SitesList, MLI,Site,MLIBegin,*Sites;
int MinSitesY,MaxSitesN;
int CEMNum,NumofCE,CEprfNum=0,TotalNumofCE=0,CEnum;
CETableElemP *CEDistrY,*CEDistrN;
float BestCoverAbs, BestCoverRel;

//технические переменные, без отношения к данному методу
FILE *inf1,*outf,*compel;
char outfileb[50],mmdatfile[100],*infile,*outfile,*OutStr, *libfile;
int outfb,MCatchoutb;
struct stat fst;
time_t T0,T1;
char SeqChars[5]={'a','c','g','t','n'}, CoreChars[5]={'A','C','G','T','N'};
char Name[MaxLineLength];
short int buffer[5], *BigBuffer;
char *SeqT;
float y,f1,f2,f3,FreqPer1Kb,ProfileAllValues=0,ProfileTuneBy=0;
int GContent[5];
int i,j,l,k,m,n,sum,f,SN,N,Strand,i1,i2,i3,i4,printed,BestI1,BestI2,BestI3,BestI4;
char s[MaxLineLength];
int needed[10];

// параметры метода
int Verbose=0, MakeProfile=0, FullOutput=0, WindowLength=200; // длиннна окна
//int PwmThrNum=20, LenThrNum=10, CSThrNum=20; // defined above
float MinCoverY=0.75,MaxCoverN=0.5,FreqDiff=1.5,CoverDiff=1.5;
float PwmThrRlx=0.4, LenThrRlx=3, CSThrRlx=0.8; // Maximum relaxation values. PwmThrMin=1-PwmThrRlx
//float PwmThrMin=0.6, LenThrMax=3, CSThrMin=1.2; // PwmThrMin=1-PwmThrRlx
//float Thr[20]={0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60};

//ProfileListP Profile[MaxMatrixNum],Prf;
//CE_ProfileElemP CEProfile[MaxNumCompElem],CEPrf;
//int BestFreqL,BestFreqR,Fr,Pos,BestPosL,BestPosR,BestMatrL,BestMatrR,Strand,BestMatrLStrand,BestMatrRStrand;
//float BestWeightR,BestWeightL,Weight,val;

#include "MMfunctions.h"

int FindMatrixbyName(char * MatrixName)
{
for (l=0;l<=NumberOfMatrix-1;l++)
 if (strcmp(MatrixList[l]->Name,MatrixName)==0) {return (l);}
printf ("Matrix not found: %s\n",MatrixName);
printf ("Matrix names in .dat:\n "); for (l=0;l<=NumberOfMatrix-1;l++) printf ("%s,",MatrixList[l]->Name);

exit(0);
}


 //*************************************************************************************************************************************************//
//*************************************************************************************************************************************************//
int main(int argc, char** argv)
{
inittime
printf ("Module Statistics v.1.0 - Okt 14 2010 Igor Deyneko\n");

if (argc < 4)
         {fprintf(outerr,"Error in parameters\n");
	       fprintf(outerr,"Usage: mmModules ModulesDatFile YesSites.bin NoSites.bin  [FreqDiff float MaxCoverN float MinCoverY float] [/v 0..3 - Verbose level] [/window bp - window length]/FullOutput - print found sites\n"); exit(0);}

libfile=argv[1];
YesSitesF=argv[2];
NoSitesF =argv[3];
//infile=argv[2];
//outfile=argv[3]; // ?????
ReadOptionalParameters (argc,argv);
printf ("Parameters: %s, %s, %s      ",libfile,YesSitesF,NoSitesF);
printf ("MinCoverY %3.2f, MaxCoverN %3.2f, FreqDiff %3.2f WindowLength %d\n",MinCoverY,MaxCoverN,FreqDiff,WindowLength);

//strcpy(mmdatfile,"G:\\igor\\programs\\MM\\Dat\\405595b08c7f9_matrixTFP81.lib.dat");// mm.dat");
strcpy(mmdatfile,"mm.dat"); printf ("Use local copy of mm.dat\n");
NumberOfMatrix=LoadAllMatrix (MatrixList, mmdatfile);
fprintf (outerr,"Number of matr=%d\n",NumberOfMatrix);

// прочитаем сайты
//YES
InfoY=(InfBlockP)malloc(sizeof(InfBlock));
#ifdef Win32
f=open(YesSitesF, O_RDONLY|O_BINARY);
#else
f=open(YesSitesF, O_RDONLY);
#endif
read (f,InfoY,sizeof(InfBlock)); // информационный блок

SitesListY=(SiteElemP**)malloc(sizeof(SiteElemP*)*InfoY->NumofSeq);
TotalNumOfSitesY=0;

//printf ("Start reading sites... ");
for (SN=0;SN<=InfoY->NumofSeq-1;SN++) // Большой цикл по кол-ву сиквенсов
{
NumOfSites=InfoY->NumofSites[SN];
SeqLength=InfoY->SeqLength [SN];
if ((NumOfSites)&&(!(SitesListY[SN]=(SiteElemP*)malloc(sizeof(SiteElemP)*NumOfSites)))) printf ("Memory SitesList! NumOfSites=%d\n",NumOfSites);
if (!(MLIBegin=MLI=(SiteElemP)(malloc(sizeof(SiteElem)*(NumOfSites+1))))) printf ("Memory MLIBegin!\n");
//printf ("Reading sites... ");
for (i=0; i<=NumOfSites-1;i++)
{Site=SitesListY[SN][i]=MLI++;
read (f,&buffer,10);         //!here 10 bit per site!
if ( buffer[0] > 0) {Site->ClassID=buffer[0]-1; Site->Strand=0;} else {Site->ClassID=-1-buffer[0];Site->Strand=1;} //memcpy (buffer+1,&Site->Begin,4);
memcpy (&Site->Begin,buffer+1,4); memcpy (&Site->Score,buffer+3,4); // Here it is absolute score!
} //all sites this seq read
SortSitesbyBegin (SitesListY[SN],NumOfSites,SeqLength);
TotalNumOfSitesY+=NumOfSites;TotalSeqLengthY+=SeqLength;
}
close (f);
printf ("NumofSeqY=%d, TotalNumOfSitesY=%d ",InfoY->NumofSeq,TotalNumOfSitesY);
printtimenL

//NO
InfoN=(InfBlockP)malloc(sizeof(InfBlock));
#ifdef Win32
f=open(NoSitesF, O_RDONLY|O_BINARY);
#else
f=open(NoSitesF, O_RDONLY);
#endif
read (f,InfoN,sizeof(InfBlock)); // информационный блок
SitesListN=(SiteElemP**)malloc(sizeof(SiteElemP*)*InfoN->NumofSeq);
TotalNumOfSitesN=0;

for (SN=0;SN<=InfoN->NumofSeq-1;SN++) // Большой цикл по кол-ву сиквенсов
{
NumOfSites=InfoN->NumofSites[SN]; //printf ("NumOfSites= %d\n",NumOfSites);
SeqLength=InfoN->SeqLength [SN];
if ((NumOfSites)&&(!(SitesListN[SN]=(SiteElemP*)malloc(sizeof(SiteElemP)*NumOfSites)))) printf ("Memory SitesListN! NumOfSites=%d\n",NumOfSites);
if (!(MLIBegin=MLI=(SiteElemP)(malloc(sizeof(SiteElem)*(NumOfSites+1))))) printf ("Memory MLIBegin, asked for %d bytes!\n",sizeof(SiteElem)*(NumOfSites+1));
//printf ("Reading sites... ");
for (i=0; i<=NumOfSites-1;i++)
   {Site=SitesListN[SN][i]=MLI++;
	 read (f,&buffer,10);
	 if ( buffer[0] > 0) {Site->ClassID=buffer[0]-1; Site->Strand=0;} else {Site->ClassID=-1-buffer[0];Site->Strand=1;} //memcpy (buffer+1,&Site->Begin,4);
	 memcpy (&Site->Begin,buffer+1,4); memcpy (&Site->Score,buffer+3,4); // Here it is absolute score!
   }
SortSitesbyBegin (SitesListN[SN],NumOfSites,SeqLength);
TotalNumOfSitesN+=NumOfSites;TotalSeqLengthN+=SeqLength;
}
close (f);
printf ("NumofSeqN=%d, TotalNumOfSitesN=%d ",InfoN->NumofSeq, TotalNumOfSitesN);
printtimenL
// ***********  Сайты готовы!  *************

// прочитаем модели  - лучше читать текстовик :(
compel=fopen(libfile,"r"); //
fgets (s,MaxLineLength,compel);
fgets (s,MaxLineLength,compel);
sscanf (s,"Number of Module models = %d\n",&CEMNum);
for (j=0;j<=CEMNum-1;j++)
    {CE=CEMList[j]=new ModuleElem;
     fgets (s,MaxLineLength,compel);
     if (sscanf (s,"M %d = %s %s\n",&CE->AC,CE->MatrName1, CE->MatrName2)!=3)   printf ("Bad model: %s\n",s);
//relax parameters
     CE->MatrID1=FindMatrixbyName(CE->MatrName1);   CE->MatrID2=FindMatrixbyName(CE->MatrName2);
     CE->Score1=MatrixList[CE->MatrID1]->Max;    CE->Score2=MatrixList[CE->MatrID2]->Max;
     CE->Min=CE->Max=CE->MaxRelaxed=CE->MinRelaxed=WindowLength;

     CE->Score1Relaxed=(int)(((float)CE->Score1)*((float)1-PwmThrRlx));
     CE->Score2Relaxed=(int)(((float)CE->Score2)*((float)1-PwmThrRlx));

     CE->TotalNum=0; //
    }
fclose (compel);
printf ("Number of Modules = %d\n",CEMNum);

//таблица рапределения КЭ CEDistr [Thr1][Thr2][L][CS]
CEDistrY=(CETableElemP*)malloc(sizeof(CETableElemP)*(InfoY->NumofSeq+3)); // +0 там будем хранить сумму сайтов по сиквенсам,+1 годность по ограничениям +2 и т.д.
for (m=0;m<=InfoY->NumofSeq+3-1;m++) {CEDistrY[m] = (CETableElemP)malloc(sizeof(CETableElem));}

CEDistrN=(CETableElemP*)malloc(sizeof(CETableElemP)*(InfoN->NumofSeq+3)); // +0 там будем хранить сумму сайтов по сиквенсам,+1 годность по ограничениям +2 и т.д.
for (m=0;m<=InfoN->NumofSeq+3-1;m++) {CEDistrN[m] = (CETableElemP)malloc(sizeof(CETableElem));}


// Всё готово можно начинать подсчет
// Цикл по всем моделям КЭ

for (CEnum=0;CEnum<=CEMNum-1;CEnum++)
  {//CEnum=42;
   CE=CEMList[CEnum];CER=CEMRList[CEnum];
   //NumofCE=0;
   fprintf (outerr,"%3d Module %3d ",CEnum,CE->AC); printtimeOutErrL

//обнулим таблицу распределения КЭ в сиквенсах
for (SN=0;SN<=InfoY->NumofSeq+3-1;SN++) //по всем сиквенсам
{ for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) CEDistrY[SN]->CETable[i1][i2][i3][i4]=0; //занулить что там есть
     //потом посмотрим на ск лучше memcpy (CEDistrN[m],CEDistrY[m],sizeof(CETableElem));
}
for (SN=0;SN<=InfoN->NumofSeq+3-1;SN++) //по всем сиквенсам
{ for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) CEDistrN[SN]->CETable[i1][i2][i3][i4]=0; //занулить что там есть
}

//заполним таблицу рапределения КЭ
//YES
//поскольку сразу ищем 1й или 2й сайты (т.е. прямой и обратный модуль), то надо учитывать порядок при присваивании i1 i2 !!
for (SN=0;SN<=InfoY->NumofSeq-1;SN++) //по всем сиквенсам
 { SitesList = SitesListY[SN];
   NumOfSites=InfoY->NumofSites[SN];
   //for (m=0;m<=InfoY->NumofSites[SN]-1;m++) // по всем сайтам на этом сикв

j=0;
needed[0]=1; needed[1]=1; needed[2]=1;
if (NumOfSites<=1) goto NextSeqY; //if there is sites at all
do
{
   if ((((SitesList[j]->ClassID == CE->MatrID1)&&(SitesList[j]->Score >= CE->Score1Relaxed))&&(m=1))||  // Нашли первый и пометили какой нашли
       (((SitesList[j]->ClassID == CE->MatrID2)&&(SitesList[j]->Score >= CE->Score2Relaxed))&&(m=2)))
     {needed[m-1]=0;
      k=j+1;
      while ((k<=NumOfSites-1)&&(SitesList[k]->Begin-SitesList[j]->Begin <= CE->MaxRelaxed))   // второй недалеко
         {if ((((needed[0])&&(SitesList[k]->ClassID == CE->MatrID1)&&(SitesList[k]->Score >= CE->Score1Relaxed)))||  // Нашли первый и пометили какой нашли
              (((needed[1])&&(SitesList[k]->ClassID == CE->MatrID2)&&(SitesList[k]->Score >= CE->Score2Relaxed))))
             { // нашли!
//             fprintf (outusr1,"CE00%03d %s %4.3f/%d %3d %s %4.3f/%d %d+ %4.3f %e ",CE->AC, MatrixList[CE->MatrIDL]->Name,SitesList[j]->Score/(float)MatrixList[CE->MatrIDL]->Max,CE->StrandL,
 //			                      SitesList[l]->Begin-SitesList[j]->Begin,MatrixList[CE->MatrIDR]->Name,SitesList[l]->Score/(float)MatrixList[CE->MatrIDR]->Max,CE->StrandR,SitesList[j]->Begin+1,x,val); // нашли второй

             // put this element in the table CEDistrY[SN]->CETable[][][][]
//   MAX OF 2 SCORES          i1 = (int)floor( (1-max(SitesList[j]->Score/(float)MatrixList[SitesList[j]->ClassID]->Max,SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max))*(PwmThrNum-1)/PwmThrRlx); // 1st coordinate, score of the 1st site
  //               if ((i1<0)||(i1>=PwmThrNum)) {printf ("Error PwmThrNum1 %d: %d %d, %d %d\n",i1,SitesList[j]->Score,MatrixList[SitesList[j]->ClassID]->Max,SitesList[k]->Score,MatrixList[SitesList[k]->ClassID]->Max);exit(0);}

/* РАЗОБРАТЬСЯ ->         ?? неправильно !!! - все порешал, PwmThrNum-1 везде, тк. отрезков -1, а точек PwmThrNum */
             i1 = (int)ceil( (1-SitesList[j]->Score/(float)MatrixList[SitesList[j]->ClassID]->Max)*(PwmThrNum-1)/PwmThrRlx); // 1st coordinate, score of the 1st site
                 if (i1 == PwmThrNum ) {i1=PwmThrNum-1;} // Иногда случается, когда при округлении целый скор становиться чуть меньше и тогда относ скор меньше порога: 0.9->0.9*12->10.8=10->10/12=0.83 12=MaxScore
                 if ((i1<0)||(i1>=PwmThrNum)) {printf ("Y: Error PwmThrNum1 %d: %d %d, %d %d\n",i1,SitesList[j]->Score,MatrixList[SitesList[j]->ClassID]->Max,SitesList[k]->Score,MatrixList[SitesList[k]->ClassID]->Max); exit(0);}

             i2 = (int)ceil( (1-SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max)*(PwmThrNum-1)/PwmThrRlx); // 1st coordinate, score of the 1st site
//             double y1=((double)SitesList[k]->Score)/(double)MatrixList[SitesList[k]->ClassID]->Max;
  //           printf ("%d %d = %f\n",SitesList[k]->Score,MatrixList[SitesList[k]->ClassID]->Max,((double)SitesList[k]->Score)/(double)MatrixList[SitesList[k]->ClassID]->Max);
                 if (i2 == PwmThrNum ) {i2=PwmThrNum-1;} // Иногда случается, когда при округлении целый скор становиться чуть меньше и тогда относ скор меньше порога: 0.9->0.9*12->10.8=10->10/12=0.83 12=MaxScore
    //             if ((i2<0)||(i2>=PwmThrNum)) {printf ("Error PwmThrNum1 %d: %d %d, %d %d\n",i2,SitesList[j]->Score,MatrixList[SitesList[j]->ClassID]->Max,SitesList[k]->Score,MatrixList[SitesList[k]->ClassID]->Max);exit(0);}

             if (m==2) {l=i1;i1=i2;i2=l;} // т.к. первым был найден второй элемент модуля!!!!
             CEDistrY[SN]->CETable[i1][i2][0][0]++;
             //NumofCE++;CE->TotalNum++;
             }
          k++;
         }
     needed[m-1]=1; //restore needed[]
     } //

j++;
}while (j<=NumOfSites-2);
NextSeqY:;
} //next seq

//NO
for (SN=0;SN<=InfoN->NumofSeq-1;SN++) //по всем сиквенсам
 { SitesList = SitesListN[SN];
   NumOfSites=InfoN->NumofSites[SN];
   //for (m=0;m<=InfoY->NumofSites[SN]-1;m++) // по всем сайтам на этом сикв

j=0;
needed[0]=1; needed[1]=1; needed[2]=1;
if (NumOfSites<=1) goto NextSeqN; //if there is sites at all
do
{
   if ((((SitesList[j]->ClassID == CE->MatrID1)&&(SitesList[j]->Score >= CE->Score1Relaxed))&&(m=1))||  // Нашли первый и пометили какой нашли
       (((SitesList[j]->ClassID == CE->MatrID2)&&(SitesList[j]->Score >= CE->Score2Relaxed))&&(m=2)))
     {needed[m-1]=0;
      k=j+1;
      while ((k<=NumOfSites-1)&&(SitesList[k]->Begin-SitesList[j]->Begin <= CE->MaxRelaxed))   // второй недалеко
         {if ((((needed[0])&&(SitesList[k]->ClassID == CE->MatrID1)&&(SitesList[k]->Score >= CE->Score1Relaxed)))||  // Нашли первый и пометили какой нашли
              (((needed[1])&&(SitesList[k]->ClassID == CE->MatrID2)&&(SitesList[k]->Score >= CE->Score2Relaxed))))
             { // нашли!
//             fprintf (outusr1,"CE00%03d %s %4.3f/%d %3d %s %4.3f/%d %d+ %4.3f %e ",CE->AC, MatrixList[CE->MatrIDL]->Name,SitesList[j]->Score/(float)MatrixList[CE->MatrIDL]->Max,CE->StrandL,
 //			                      SitesList[l]->Begin-SitesList[j]->Begin,MatrixList[CE->MatrIDR]->Name,SitesList[l]->Score/(float)MatrixList[CE->MatrIDR]->Max,CE->StrandR,SitesList[j]->Begin+1,x,val); // нашли второй

             // put this element in the table CEDistrY[SN]->CETable[][][][]
//  max of ..           i1 = (int)floor((1-max(SitesList[j]->Score/(float)MatrixList[SitesList[j]->ClassID]->Max,SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max))*(PwmThrNum-1)/PwmThrRlx); // 1st coordinate, score of the 1st site
  //               if ((i1<0)||(i1>=PwmThrNum)) {printf ("Error PwmThrNumN1 %d\n",i1);exit(0);}
             i1 = (int)ceil((1-SitesList[j]->Score/(float)MatrixList[SitesList[j]->ClassID]->Max)*(PwmThrNum-1)/PwmThrRlx); // 1st coordinate, score of the 1st site
                 if (i1 == PwmThrNum ) {i1=PwmThrNum-1;} // Иногда случается, когда при округлении целый скор становиться чуть меньше и тогда относ скор меньше порога: 0.9->0.9*12->10.8=10->10/12=0.83 12=MaxScore
                 if ((i1<0)||(i1>=PwmThrNum)) {printf ("NO: Error PwmThrNumN1 %d\n",i1);exit(0);}
             i2 = (int)ceil((1-SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max)*(PwmThrNum-1)/PwmThrRlx); // 1st coordinate, score of the 1st site
                 if (i2 == PwmThrNum ) {i2=PwmThrNum-1;} // Иногда случается, когда при округлении целый скор становиться чуть меньше и тогда относ скор меньше порога: 0.9->0.9*12->10.8=10->10/12=0.83 12=MaxScore
                 if ((i2<0)||(i2>=PwmThrNum)) {printf ("NO: Error PwmThrNumN1 %d\n",i2);exit(0);}

             if (m==2) {l=i1;i1=i2;i2=l;} // т.к. первым был найден второй элемент модуля!!!!
             CEDistrN[SN]->CETable[i1][i2][0][0]++;
             //NumofCE++;CE->TotalNum++;
             }
          k++;
         }
     needed[m-1]=1; //restore needed[]
     } //

j++;
}while (j<=NumOfSites-2);

NextSeqN:;
} //next seq

//printf("Distribution done "); printtime

//Ну а теперь проверки и поиск
//суммируем по порогам
for (SN=0;SN<=InfoY->NumofSeq-1;SN++) //по всем сиквенсам
{ for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=1;i4<=CSThrNum -1;i4++) CEDistrY[SN]->CETable[i1][i2][i3][i4]+=CEDistrY[SN]->CETable[i1][i2][i3][i4-1];
  for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i4=0;i4<=CSThrNum -1;i4++) for (i3=1;i3<=LenThrNum-1;i3++) CEDistrY[SN]->CETable[i1][i2][i3][i4]+=CEDistrY[SN]->CETable[i1][i2][i3-1][i4];
  for (i1=0;i1<=PwmThrNum-1;i1++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) for (i2=1;i2<=PwmThrNum-1;i2++) CEDistrY[SN]->CETable[i1][i2][i3][i4]+=CEDistrY[SN]->CETable[i1][i2-1][i3][i4];
  for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) for (i1=1;i1<=PwmThrNum-1;i1++) CEDistrY[SN]->CETable[i1][i2][i3][i4]+=CEDistrY[SN]->CETable[i1-1][i2][i3][i4];
} // потом переставить порядок суммирования - непомогло :(

for (SN=0;SN<=InfoN->NumofSeq-1;SN++) //по всем сиквенсам
{ for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=1;i4<=CSThrNum -1;i4++) CEDistrN[SN]->CETable[i1][i2][i3][i4]+=CEDistrN[SN]->CETable[i1][i2][i3][i4-1];
  for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i4=0;i4<=CSThrNum -1;i4++) for (i3=1;i3<=LenThrNum-1;i3++) CEDistrN[SN]->CETable[i1][i2][i3][i4]+=CEDistrN[SN]->CETable[i1][i2][i3-1][i4];
  for (i1=0;i1<=PwmThrNum-1;i1++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) for (i2=1;i2<=PwmThrNum-1;i2++) CEDistrN[SN]->CETable[i1][i2][i3][i4]+=CEDistrN[SN]->CETable[i1][i2-1][i3][i4];
  for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) for (i1=1;i1<=PwmThrNum-1;i1++) CEDistrN[SN]->CETable[i1][i2][i3][i4]+=CEDistrN[SN]->CETable[i1-1][i2][i3][i4];
} // потом переставить порядок суммирования

//суммируем по сиквенсам. Кол-во сайтов в выборке при данных порогах, результат кладем в InfoY->NumofSeq+0
for (SN=0;SN<=InfoY->NumofSeq-1;SN++) //по всем сиквенсам
{ for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) CEDistrY[InfoY->NumofSeq]->CETable[i1][i2][i3][i4]+=CEDistrY[SN]->CETable[i1][i2][i3][i4]; }
for (SN=0;SN<=InfoN->NumofSeq-1;SN++) //по всем сиквенсам
{ for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++) CEDistrN[InfoN->NumofSeq]->CETable[i1][i2][i3][i4]+=CEDistrN[SN]->CETable[i1][i2][i3][i4]; }

//printf("Sum done "); printtime
//проверим проходит ли по ограничениям,
// 1) Покрытие выборки, результат кладем в +1: n-Да, 0-Нет.
// MinSites - мин чилсло сайтов на 1секв. MinCover - сколько таких сикв в выборке
MinSitesY=1;
//MinCoverY=0.75;
for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++)
{n=0; for (SN=0;SN<=InfoY->NumofSeq-1;SN++) //по всем сиквенсам
           if (CEDistrY[SN]->CETable[i1][i2][i3][i4]>=MinSitesY) n++;
      if (n/(float)InfoY->NumofSeq >= MinCoverY) CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4]=n;
}

// No покрывается не более чем на MaxCoverN сайтами
MaxSitesN=1;
//MaxCoverN=0.5;
for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++)
{n=0; for (SN=0;SN<=InfoN->NumofSeq-1;SN++) //по всем сиквенсам
           if (CEDistrN[SN]->CETable[i1][i2][i3][i4]>=MaxSitesN) n++;
      if (n/(float)InfoN->NumofSeq <= MaxCoverN) if (n==0) CEDistrN[InfoN->NumofSeq+1]->CETable[i1][i2][i3][i4]=-1; else CEDistrN[InfoN->NumofSeq+1]->CETable[i1][i2][i3][i4]=n;
}

// 2) Разница по частоте сайтов в выборке, результат кладем в +2: 1-Да, 0-Нет.
//FreqDiff=2.0; // Должны отличаться по частотам мин в 2 раза

for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++)
{
if (CEDistrY[InfoY->NumofSeq]->CETable[i1][i2][i3][i4] == 0) CEDistrY[InfoY->NumofSeq+2]->CETable[i1][i2][i3][i4] = 0; else
if (CEDistrN[InfoN->NumofSeq]->CETable[i1][i2][i3][i4] == 0) CEDistrY[InfoY->NumofSeq+2]->CETable[i1][i2][i3][i4] = 1; else
   {f1=CEDistrY[InfoY->NumofSeq]->CETable[i1][i2][i3][i4]/(float)TotalSeqLengthY;
    f2=CEDistrN[InfoN->NumofSeq]->CETable[i1][i2][i3][i4]/(float)TotalSeqLengthN;
    if ((f1/f2 >= FreqDiff)||(f2/f1 >= FreqDiff)) CEDistrY[InfoY->NumofSeq+2]->CETable[i1][i2][i3][i4] = 1; }
}

//printf("Restrictions done "); printtime
//выводим результат
//printf ("Parameters: MinSitesY %d, MinCoverY %3.2f, MaxCoverN %3.2f, FreqDiff %3.2f \n",MinSitesY,MinCoverY,MaxCoverN,FreqDiff);
printed=0;
BestCoverAbs=0; BestCoverRel=0;
for (i1=0;i1<=PwmThrNum-1;i1++) for (i2=0;i2<=PwmThrNum-1;i2++) for (i3=0;i3<=LenThrNum-1;i3++) for (i4=0;i4<=CSThrNum -1;i4++)

     if ((CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4]) &&(CEDistrN[InfoN->NumofSeq+1]->CETable[i1][i2][i3][i4])&&(CEDistrY[InfoY->NumofSeq+2]->CETable[i1][i2][i3][i4]))
//     if ((CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4]) )

     { // print this one
       f1=CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4]/(float)InfoY->NumofSeq;
       f2=CEDistrN[InfoN->NumofSeq+1]->CETable[i1][i2][i3][i4]/(float)InfoN->NumofSeq;
       if (f2 < 0)
       { if ( (BestCoverAbs/f1-1) < 0.001) // better case
         {if (fabs(BestCoverAbs/f1-1) < 0.001 ) // same case, only to check if i1,i2,.. are better.
             {  if (UniformPoints(i1,i2,i3,i4,BestI1,BestI2,BestI3,BestI4)) {BestI1=i1;BestI2=i2;BestI3=i3;BestI4=i4;} else goto NextPoint1;
             }
             else // really better point
             { BestI1=i1;BestI2=i2;BestI3=i3;BestI4=i4;
               BestCoverAbs=f1;}
//          printf ("CE %3d|Thr=%d,%d,%d,%d|%4.3f/0.000 (%d/0)|AbsThr=%4.3f, %4.3f, %4.3f, %4.3f\n",CE->AC,i1,i2,i3,i4,f1,CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4],(i1*PwmThrRlx)/(float)PwmThrNum,(i2*PwmThrRlx)/(float)PwmThrNum,(i3*LenThrRlx)/(float)LenThrNum,(i4*CSThrRlx)/(float)CSThrNum);
          printf ("ModuleID %3d|Thr=%2d,%2d|Cover=%4.3f/0.000 (%d/0)|AbsThr=%4.3f, %4.3f\n",CE->AC,i1,i2,f1,CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4],                                                           (i1*PwmThrRlx)/(float)(PwmThrNum-1),(i2*PwmThrRlx)/(float)(PwmThrNum-1));
//Выведем сиквенсы в которых есть такие элементы
             if (FullOutput>=2)
             for (i=0;i<=InfoY->NumofSeq-1;i++) //по всем сиквенсам
              if (CEDistrY[i]->CETable[i1][i2][i3][i4]>=1)
                  {printf ("Seq%04d\n",i);
                  if (FullOutput>=3)   // выведем сайты которые подходят под этот модуль
                  { SitesList = SitesListY[i];
                  	NumOfSites=InfoY->NumofSites[i];
                  	for (k=0;k<=NumOfSites-1;k++)
                   	{ if ((SitesList[k]->ClassID == CE->MatrID1)||
                         (SitesList[k]->ClassID == CE->MatrID2))  //&&(SitesList[k]->Score >= (i1*PwmThrRlx)/(float)PwmThrNum
                     	{  printf ("f0 %s %d\n", MatrixList[SitesList[k]->ClassID]->Name, SitesList[k]->Begin);
	                     }

   	               }
                  }
                  }
 //       printed =1;
 //       if (printed) goto NextCE;
          //BestCoverRel=2; // теперь ищем только те где f2<0
          //  и будем улучшать f1
NextPoint1:;
         }
       }
       else
       { f3=f1/f2;  // f1/f2 should first go in mem, so that it is rounded and can be compared!!!
        if (( (BestCoverRel/f3-1) < 0.001)&&(BestCoverAbs == 0))
         {if (fabs(BestCoverRel/f3-1) < 0.001 ) // same case, only to check if i1,i2,.. are better.
             {  if (UniformPoints(i1,i2,i3,i4,BestI1,BestI2,BestI3,BestI4)) {BestI1=i1;BestI2=i2;BestI3=i3;BestI4=i4;} else goto NextPoint2;
             }
             else // really better point
             { BestI1=i1;BestI2=i2;BestI3=i3;BestI4=i4;
               BestCoverRel=f3;}

//             printf ("ModuleID %3d|Thr=%d,%d,%d,%d|%4.3f/%4.3f (%d/%d)|AbsThr=%4.3f, %4.3f, %4.3f, %4.3f\n",CE->AC,i1,i2,i3,i4,f1,f2,CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4],CEDistrN[InfoN->NumofSeq+1]->CETable[i1][i2][i3][i4],(i1*PwmThrRlx)/(float)PwmThrNum,(i2*PwmThrRlx)/(float)PwmThrNum,(i3*LenThrRlx)/(float)LenThrNum,(i4*CSThrRlx)/(float)CSThrNum);
             printf ("ModuleID %3d|Thr=%2d,%2d|Cover=%4.3f/%4.3f (%d/%d)|AbsThr=%4.3f, %4.3f\n",CE->AC,i1,i2,f1,f2,CEDistrY[InfoY->NumofSeq+1]->CETable[i1][i2][i3][i4],CEDistrN[InfoN->NumofSeq+1]->CETable[i1][i2][i3][i4],(i1*PwmThrRlx)/(float)(PwmThrNum-1),(i2*PwmThrRlx)/(float)(PwmThrNum-1));
//Выведем сиквенсы в которых есть такие элементы
             if (FullOutput>=2)
             for (i=0;i<=InfoY->NumofSeq-1;i++) //по всем сиквенсам
              if (CEDistrY[i]->CETable[i1][i2][i3][i4]>=1)
                  { printf ("Seq%04d\n",i);
                  if (FullOutput)   // выведем сайты которые подходят под этот модуль
                    {SitesList = SitesListY[i];
                  	NumOfSites=InfoY->NumofSites[i];
                  	for (k=0;k<=NumOfSites-1;k++)
                   		{    // print all sites for testing only
	                     //if ((SitesList[k]->ClassID == CE->MatrID1)||(SitesList[k]->ClassID == CE->MatrID2)) printf ("%20s %3d-%3d %4.3f\n", MatrixList[SitesList[k]->ClassID]->Name, SitesList[k]->Begin+1,SitesList[k]->Begin+MatrixList[SitesList[k]->ClassID]->OrigLength,(1-SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max) );

      	               if (((SitesList[k]->ClassID == CE->MatrID1)&&((1-SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max) <= (i1*PwmThrRlx)/(float)(PwmThrNum-1)))||
                         ((SitesList[k]->ClassID == CE->MatrID2)&&((1-SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max) <= (i2*PwmThrRlx)/(float)(PwmThrNum-1))))
         		            {  printf ("%20s %3d-%3d %4.3f\n", MatrixList[SitesList[k]->ClassID]->Name, SitesList[k]->Begin+1,SitesList[k]->Begin+MatrixList[SitesList[k]->ClassID]->OrigLength,(1-SitesList[k]->Score/(float)MatrixList[SitesList[k]->ClassID]->Max) );
               		      }

                   		}
                    }
                  }
 //            printed =1;
//        if (printed) goto NextCE;
//             BestCoverRel=f3;
NextPoint2:;
          }
       }
     }



NextCE:;
//printf("Output done "); printtime
//printf("\n");
} //Цикл по всем моделям КЭ

//printtime
//fprintf (outerr,"All done\n");
printf ("All done\n");

exit(0);

}






