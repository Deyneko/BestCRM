
#define Win32
//#define Linux

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

#define outusr stdout
#define outerr stdout //stderr
#define MaxLineLength 200       // in any file
#define MaxSeqLength 2000000000   // Max sequence lenth
#define MinSeqLength 50 // Min
#define MaxMatrLength 45 // Max length of matrix
#define MaxMatrixNum 2000 // Max number of matrices
#define MaxMatrixNameLen 25 //
#define MaxNumofSeq 50000 // максимальное кол-во сиквенсов в файле
#define inittime {T0=time(NULL);}
#define printtime //{T1=time(NULL);printf ("Time - %li\n",T1-T0);T0=T1;}

// Here define output format
//#define MCatchOutPut  - используется только в MCatch.  MCatchtools уже читают простой бинари формат.
//#define Text
#define Binary
//#define EncodeSeq
//#define MCatchFreqTable - больше не использую. Теперь просто ещем сайты, потом mmtools делаем таблицу.
//#define CGreport
//#define Triplets
//#define CountEntropyProfile

//#define UseDefaultDir  - сохранять .дат файлы в спец директории.

typedef struct Penta_
        {int Pos; Penta_* Next;} Penta, *PentaP;
typedef struct MatrHash_
        {int AccN; int Val;} MatrHash, *MatrHashP;
typedef struct MatrixElem_
	{ char Name[MaxMatrixNameLen];
     int AC,ID;
     int OrigLength,Length;
     int Max,Min,CoreMin,CoreMax,CoreThr,Thr;
     float Thresh;
     int Val[5][MaxMatrLength], ValR[5][MaxMatrLength];
     MatrHash H[3125], HR[3125];
     int CoreStart, CoreStartR; // сколько отступать от кора до краев матрицы
     int LeftSide, RightSide; // кол-во троек по сторонам
     int Shift1, Shift2; // дополнительных мест в расширенной матрице слева, справа
     int ValT[125][(int)MaxMatrLength/3+1],ValTR[125][(int)MaxMatrLength/3+1]; // Значения троек
//     float StatThresh [BestSitesNum];
//     float Sigma[BEstSitesNum];
//     float Mu[BestSitesNum];
   } MatrixElem, * MatrixElemP;
typedef struct ProfileList_
   { float CoreThr,Thr;
     char Name [MaxMatrixNameLen] ,Num[10];} ProfileList, *ProfileListP;

typedef struct InfBlock_
	{int NumofSeq;
    int SeqLength [MaxNumofSeq];
    int NumofSites [MaxNumofSeq]; } InfBlock, *InfBlockP;

enum formats {EMBLf,FASTAf,PlainSeqf,UnDef} format=UnDef;
int SeqLength,xl; // если unsigned  то когда Seq[-5] то это уже не -5!
char Name[MaxLineLength],*Seq,*infile1,*outfile,*libfile;
FILE* inf1,*outf;
int x,i,j,l,sum,outbin;
PentaP* PentaSeq;
MatrixElemP MatrixList [MaxMatrixNum];
int NumberOfMatrix;
MatrixElemP M;
PentaP P;
time_t T0,T1;
int NumOfSites=0,fast=0,all=0,TotalSitesLen=0,NumofSeq=0;
char outfileb[100], mmdatfile[100];
char *SeqT;
int TotalNumOfSites=0, TotalSeqLength=0,CurrentNumOfSites=0,TotalSiteStat[MaxMatrixNum];
int OnSeqOnly=0,Verbose=0;
char OutStr[MaxLineLength];
int outfb,MCatchoutb,Strand;;
char SeqChars[5]={'a','c','g','t','n'};
char CoreChars[5]={'A','C','G','T','N'};
float y,ProfileTuneBy=0,ProfileAllValues=0;
InfBlockP Info;
int GContent[5], DiContent[25];
struct stat fst, fst1;
short int buffer[16];
short int *BigBuffer;
int BigBufferSize=1024*1024/2+10, BigBufferEnd=0; //BufferSize=500000 in items (=*2 bytes) 524288*2+20=1048576+20 = 1 Mb +20, 20 это защитный участок, пишем ровно по 1мег!
const char *DefaultDir="G:\\igor\\programs\\MM\\Dat";
char **SeqSet;
int *SeqLen,GCWindow;
float **GCDistr;
float  GCAvrg[5], GCStDev[5],GCThr,ATThr;
char *SeqNames[MaxNumofSeq];


unsigned int LoadSeq (FILE* inf, char* seq) // return value - number of char read, 0 - error
{ unsigned int i=0;
int j,k;
char ch,s[MaxLineLength],*s1;

if (format==UnDef) // in what format sequence ?
   {fgets (s,MaxLineLength,inf);
    if (s[0]=='I') format=EMBLf;
    else if (s[0]=='>') format=FASTAf;
    else format=PlainSeqf;
//    printf ("Seq format = %d\n",format);
    rewind (inf);
    for (j=-1;j>=-50;j--) seq[j]=4; // в первых 50 четверки и в конце тоже будут четверки.
   }

i=0;
if (format==PlainSeqf)
{
strcpy (Name,"No name");
while ((ch=getc(inf))&&(!feof(inf)))
	{
	  switch (ch) {
       case 'a' : seq[i++]=0;break;
       case 'A' : seq[i++]=0;break;
       case 'c' : seq[i++]=1;break;
       case 'C' : seq[i++]=1;break;
       case 'g' : seq[i++]=2;break;
       case 'G' : seq[i++]=2;break;
       case 't' : seq[i++]=3;break;
       case 'T' : seq[i++]=3;break;
       case 'n' : seq[i++]=4;break;
       case 'N' : seq[i++]=4;break;
       case ' ' : break;
       case 10  : break;
       case 13  : break;
       default  : {fprintf (outusr,"Your sequence is wrong: %c - position %u, skiped\n",ch,i); /*return (0);*/} }
	}
}
else if (format==EMBLf)
{
while ((fgets(s,MaxLineLength,inf)!=NULL)&&(strstr(s,"//")==NULL)&&(!feof(inf)))
{
if ((s[0]=='I')&&(s[1]=='D')) { k=0; while ((Name[k]=s[k+3])&&(s[k+3]!=13)&&(s[k+3]!=10)) k++; Name[k]=0;} //strcpy (Name,s+3); if (Name[strlen(Name)-1]==13) Name[strlen(Name)-1]=0; if (Name[strlen(Name)-1]==10) Name[strlen(Name)-1]=0; }
if ((s[0]=='S')&&(s[1]=='Q'))
	{while (((s1=fgets(s,MaxLineLength,inf))!=NULL)&&(strstr(s,"//")==NULL))
   {while (*s1)
	  {switch (*s1) {
       case 'a' : seq[i++]=0;break;
       case 'A' : seq[i++]=0;break;
       case 'c' : seq[i++]=1;break;
       case 'C' : seq[i++]=1;break;
       case 'g' : seq[i++]=2;break;
       case 'G' : seq[i++]=2;break;
       case 't' : seq[i++]=3;break;
       case 'T' : seq[i++]=3;break;
       case 'n' : seq[i++]=4;break;
       case 'N' : seq[i++]=4;break;
       case ' ' : break;
       case 10  : break;
       case 13  : break;
       default  : if (((*s1)<'0')||(*s1>'9')) {fprintf (outusr,"Your sequence is wrong: %c - position %u, skiped\n",*s1,i+1);/* return (0);*/}}
     s1++;
     }
   }
   goto SeqReady; //  когда много сиквенсов - то она их склеивает  - надо выходить
   }
}
}
else if (format==FASTAf)
{
fgets(s,MaxLineLength,inf);
k=0; while ((Name[k]=s[k+1])&&(s[k+1]!=13)&&(s[k+1]!=10)) k++; Name[k]=0; //strcpy (Name,s+1);
while (((ch=getc(inf))!='>')&&(!feof(inf)))
	  {switch (ch) {
       case 'a' : seq[i++]=0;break;
       case 'A' : seq[i++]=0;break;
       case 'c' : seq[i++]=1;break;
       case 'C' : seq[i++]=1;break;
       case 'g' : seq[i++]=2;break;
       case 'G' : seq[i++]=2;break;
       case 't' : seq[i++]=3;break;
       case 'T' : seq[i++]=3;break;
       case 'n' : seq[i++]=4;break;
       case 'N' : seq[i++]=4;break;
       case ' ' : break;
       case 10  : break;
       case 13  : break;
       default  : fprintf (outusr,"Your sequence is wrong: %c - position %u, skiped\n",ch,i); /*return (0);*/}
     }
if (ch=='>') ungetc('>',inf);
}

SeqReady:
for (j=0;j<=49;j++) seq[i+j]=4; // вот и концовку тоже оформили.
return (i);
};

PentaP* MakePentaSeq (char * seq, int N) // возвращ указатель на список указ на пятерки
{
PentaP Last,List,*PS,*EndsList;
int i,addr;
int j;
PS=(PentaP*)malloc(3125*sizeof(PentaP));
List=(PentaP)malloc(N*sizeof(Penta));
EndsList=(PentaP*)malloc(3125*sizeof(PentaP));
for (i=0;i<=3125-1;i++) {PS[i]=NULL;EndsList[i]=NULL;}

//for (j=0;j<=99;j++) // Для проверки
{
Last=List;
for (i=0;i<=N-1-4;i++)
  {//addr=seq[i]+seq[i+1]*5+seq[i+2]*25+seq[i+3]*125+seq[i+4]*625;
   addr=seq[i]+5*(seq[i+1]+5*(seq[i+2]+5*(seq[i+3]+5*seq[i+4])));   // так быстрее - 216/226 ;)
   if (PS[addr]==NULL) // первый раз
      {(EndsList[addr]=PS[addr]=Last++)->Pos=i;}
   else
      {(EndsList[addr]=EndsList[addr]->Next=Last++)->Pos=i;}
  }
}


for (i=0;i<=3125-1;i++) if (EndsList[i]) EndsList[i]->Next=NULL; // оформим концы

// проверим структуру PentaSeq
/*addr=0;
for (i=0;i<=SeqLength-5;i++)
     addr+=List[i].Pos;
printf ("%d\n",addr);
*/
free (EndsList);
return (PS);
};

int LoadMatrix (MatrixElemP* ML, char* libfile, char* profile)
{
FILE* inf1,inf2;
char  s[MaxLineLength];
int i=0,j,j1,k,i1,i2,i3,i4,i5,f;
MatrixElemP matr;
int TotalNumofMatr=0;
float x1,x2,x3,x4;
//struct ftime ft,ft1;
ProfileList PL[MaxMatrixNum];
int MultFactor=10000;

//struct stat fst; //
off_t FileSize; // long can be used (instead of off_t)
time_t FileDate;

if (!(strcmp(libfile,"mm.dat"))) goto FileReady; //use default lib

//check if HahedMatrix.dat exist and have proper size and date
stat(libfile,&fst);

#ifdef UseDefaultDir
//generate name for mmdatfile
if (stat(DefaultDir,&fst1) == 0) {  //check if path exists
 //Add .lib name to mm.dat
strcpy(s,libfile); //s[strlen(s)-4]=0; //And check if it was relative path
i=strlen (libfile); while (--i >= 0) {if ((libfile[i]=='/')||(libfile[i]=='\\')) {strcpy(s,libfile+i+1);/*s[strlen(s)-4]=0;*/ i=-1;} }
sprintf (mmdatfile,"%s\\%x%x_%s.dat",DefaultDir,fst.st_mtime,fst.st_size,s);
//printf ("%08x %08x ",fst.st_mtime,fst.st_size);
//printf ("sizeof(FileDate)=%u\n",sizeof(fst.st_mtime));
//strcpy(mmdatfilename,=%s%s\n",DefaultDir,s);
printf ("mmdatfilename=%s\n",mmdatfile);
}
#endif

#ifdef Win32
if ((f=open(mmdatfile, O_RDONLY|O_BINARY))==-1) {close (f);goto buildfile;}
#else
if ((f=open(mmdatfile, O_RDONLY))==-1) {close (f); goto buildfile;}//linux
#endif
read (f,&FileDate,sizeof(time_t)); if (FileDate!=fst.st_mtime) {close (f); goto buildfile;}
read (f,&FileSize,sizeof(off_t)); if (FileSize!=fst.st_size) {close (f); goto buildfile;}  //  date, length
close(f);
goto FileReady;

buildfile:
// build hashed matrix file
i=0;
fprintf (outerr,"Building matrixfile\n");
if (!(inf1=fopen(libfile,"r"))) {fprintf (outerr,"Error opening matrix file\n");exit(0);}
//matr = new MatrixElem;
matr = (MatrixElemP)malloc(sizeof(MatrixElem)); //linux
while (!(feof(inf1)))
   {fgets (s,MaxLineLength,inf1);
    s[strlen(s)-1]=0; // delete new line char
    if (s[strlen(s)-1]==0x0D) s[strlen(s)-1]=0; // for case str ends 0D 0A (we think 0A)
    if (strstr(s,"AC ")!=NULL) {matr->AC=atoi(s+4);}
    if (strstr(s,"ID ")!=NULL) {if (strlen(s+4) >= MaxMatrixNameLen) {printf ("Too long matrix name - %s\n",s+4); exit(0);} else strcpy(matr->Name,s+4);}
//    if (matr->name[0]!='F') goto NextMatr;
    if (strstr(s,"MATR_LENGTH")!=NULL) {matr->Length=atoi(s+12);matr->OrigLength=matr->Length;
                                        if (matr->Length>=(MaxMatrLength-4)) {fprintf (outerr,"Error: too long matrix - %s\n",matr->Name); exit(0);} }
    if (strstr(s,"MAXIMAL")!=NULL) {matr->Max=(int)(MultFactor*(float)atof(s+7));
                                    if (MultFactor*atof(s+7)>2E9) {printf ("Error: too big value of MultFactor %f\n",(float)atof(s+7)); exit(0);}}
    if (strstr(s,"MINIMAL")!=NULL) {matr->Min=(int)(MultFactor*(float)atof(s+7));}
    if (strstr(s,"THRESHOLD")!=NULL) {matr->Thresh=(float)atof(s+9);}
    if (strstr(s,"CORE_START")!=NULL) {matr->CoreStart=atoi(s+10)-1;} // стартуем от нуля
    if (strstr(s,"CORE_LENGTH")!=NULL) {if(atoi(s+11)!=5) {fprintf (outerr,"Different CoreLength - %d\n",atoi(s+11));exit(0);}}
    if (strstr(s,"WEIGHTS")!=NULL)
    	  for (j=0;j<=matr->Length-1;j++)
        		{fgets (s,MaxLineLength,inf1);
             if (sscanf (s,"%i A:%f C:%f G:%f T:%f\n",&k,&x1,&x2,&x3,&x4)!=5) fprintf (outerr,"Bad matrix %s\n",matr->Name);
             *(matr->Val[0]+j)=(int)(MultFactor*x1);
             *(matr->Val[1]+j)=(int)(MultFactor*x2);
             *(matr->Val[2]+j)=(int)(MultFactor*x3);
             *(matr->Val[3]+j)=(int)(MultFactor*x4);
             *(matr->ValR[3]+matr->Length-1-j)=(int)(MultFactor*x1); // обратная матрица
             *(matr->ValR[2]+matr->Length-1-j)=(int)(MultFactor*x2);
             *(matr->ValR[1]+matr->Length-1-j)=(int)(MultFactor*x3);
             *(matr->ValR[0]+matr->Length-1-j)=(int)(MultFactor*x4);

             if (x1>x2) x1=x2; if (x1>x3) x1=x3; if (x1>x4) x1=x4; // в N  лежит минимум
             *(matr->Val[4]+j)=(int)(MultFactor*x1);
             *(matr->ValR[4]+matr->Length-1-j)=(int)(MultFactor*x1);
             }
    if (strstr(s,"//")!=NULL)
        {//matr->CoreLeft=matr->CoreStart;
         matr->ID=i; // глобальный номер матрицы
         matr->CoreStartR = matr->Length - matr->CoreStart - 5; //  старт кора в обратной послед (абсолютные коорд)
         for (i5=0;i5<=4;i5++)
         for (i4=0;i4<=4;i4++)
         for (i3=0;i3<=4;i3++)
         for (i2=0;i2<=4;i2++)
         for (i1=0;i1<=4;i1++)
             {j=i1+i2*5+i3*25+i4*125+i5*625;
              matr->H[j].AccN=j;
              matr->H[j].Val=matr->Val[i1][matr->CoreStart+0]+matr->Val[i2][matr->CoreStart+1]+
                              matr->Val[i3][matr->CoreStart+2]+matr->Val[i4][matr->CoreStart+3]+
                              matr->Val[i5][matr->CoreStart+4];
              matr->HR[j].AccN=j;
              matr->HR[j].Val=matr->ValR[i1][matr->CoreStartR+0]+matr->ValR[i2][matr->CoreStartR+1]+
                               matr->ValR[i3][matr->CoreStartR+2]+matr->ValR[i4][matr->CoreStartR+3]+
                               matr->ValR[i5][matr->CoreStartR+4];
             }
         do   // отсортируем
         { i1=0;
           for (j=0;j<=3123;j++)
             if (matr->H[j].Val < matr->H[j+1].Val)
                {i2=matr->H[j].AccN; i3=matr->H[j].Val;
                    matr->H[j].AccN=matr->H[j+1].AccN; matr->H[j].Val=matr->H[j+1].Val;
                    matr->H[j+1].AccN=i2; matr->H[j+1].Val=i3;
                 i1=1;}
         }
         while (i1);
         do   //  и для второй цепи
         { i1=0;
           for (j=0;j<=3123;j++)
             if (matr->HR[j].Val < matr->HR[j+1].Val)
                {i2=matr->HR[j].AccN; i3=matr->HR[j].Val;
                    matr->HR[j].AccN=matr->HR[j+1].AccN; matr->HR[j].Val=matr->HR[j+1].Val;
                    matr->HR[j+1].AccN=i2; matr->HR[j+1].Val=i3;
                 i1=1;}
         }
         while (i1);
         matr->CoreMax=matr->H[0].Val;    // min max  для кора
         matr->CoreMin=matr->H[3124].Val;

         // переделаем Макс и Мин, (из-за всяких округлений надо)
         j=0; i3=0;
         for (i1=0;i1<=matr->Length-1;i1++)
             {i2=matr->Val[0][i1];
              if (i2 < matr->Val[1][i1]) i2 = matr->Val[1][i1];
              if (i2 < matr->Val[2][i1]) i2 = matr->Val[2][i1];
              if (i2 < matr->Val[3][i1]) i2 = matr->Val[3][i1];
              j+=i2;
              i3+=matr->Val[4][i1];
             }
         matr->Max=j; // точные макс и мин
         matr->Min=i3;

// расширим матрицу
             matr->Shift1 = (3 - matr->CoreStart%3)%3; // столько надо добавить вначале
             matr->Shift2 = (3 - (matr->Length-matr->CoreStart-5)%3)%3; //  в конце

             for (i1=matr->Length-1;i1>=0;i1--) for (i2=0;i2<=4;i2++) matr->Val[i2][i1+matr->Shift1]=matr->Val[i2][i1];
             for (i1=0;i1<=matr->Shift1-1;i1++) for (i2=0;i2<=4;i2++) matr->Val[i2][i1]=0;
             for (i1=0;i1<=matr->Shift2-1;i1++) for (i2=0;i2<=4;i2++) matr->Val[i2][i1+matr->Length+matr->Shift1]=0;

             for (i1=matr->Length-1;i1>=0;i1--) for (i2=0;i2<=4;i2++) matr->ValR[i2][i1+matr->Shift2]=matr->ValR[i2][i1];
             for (i1=0;i1<=matr->Shift2-1;i1++) for (i2=0;i2<=4;i2++) matr->ValR[i2][i1]=0;
             for (i1=0;i1<=matr->Shift1-1;i1++) for (i2=0;i2<=4;i2++) matr->ValR[i2][i1+matr->Length+matr->Shift2]=0;

             matr->Length+=(matr->Shift1+matr->Shift2);
             matr->CoreStart+=matr->Shift1;
             matr->CoreStartR+=matr->Shift2;

             matr->LeftSide = matr->CoreStart/3;
             matr->RightSide = (matr->Length-matr->CoreStart-4)/3;

             for (j=0;j<=matr->LeftSide-1;j++)  // сделаем тройки
             for (i3=0;i3<=4;i3++)
             for (i2=0;i2<=4;i2++)
             for (i1=0;i1<=4;i1++)
             matr->ValT[i1+5*i2+25*i3][j]=matr->Val[i1][j*3]+matr->Val[i2][j*3+1]+matr->Val[i3][j*3+2];
             for (j=matr->LeftSide;j<=matr->LeftSide+matr->RightSide-1;j++)
             for (i3=0;i3<=4;i3++)
             for (i2=0;i2<=4;i2++)
             for (i1=0;i1<=4;i1++)
             matr->ValT[i1+5*i2+25*i3][j]=matr->Val[i1][5+j*3]+matr->Val[i2][5+j*3+1]+matr->Val[i3][5+j*3+2];

             for (j=0;j<=matr->RightSide-1;j++)
             for (i3=0;i3<=4;i3++)
             for (i2=0;i2<=4;i2++)
             for (i1=0;i1<=4;i1++)
             matr->ValTR[i1+5*i2+25*i3][j]=matr->ValR[i1][j*3]+matr->ValR[i2][j*3+1]+matr->ValR[i3][j*3+2];
             for (j=matr->RightSide;j<=matr->LeftSide+matr->RightSide-1;j++)
             for (i3=0;i3<=4;i3++)
             for (i2=0;i2<=4;i2++)
             for (i1=0;i1<=4;i1++)
             matr->ValTR[i1+5*i2+25*i3][j]=matr->ValR[i1][5+j*3]+matr->ValR[i2][5+j*3+1]+matr->ValR[i3][5+j*3+2];

// все матрица стала побольше

         ML [i++]=matr;
        	//matr = new MatrixElem;
         matr = (MatrixElemP)malloc(sizeof(MatrixElem)); //linux
         fprintf (outerr,"Matricies processed - %d\n",i);
         }
    NextMatr: ;
   }
TotalNumofMatr=i;
fclose (inf1);
// Теперь надо все записать
/*f=open (libfile,O_RDONLY|O_TEXT); // возмем информацию о файле матриц
   i1=filelength(f);
   getftime (f,&ft);
close(f);*/
stat (libfile,&fst);
#ifdef Win32
if ((f=open(mmdatfile, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY,S_IREAD|S_IWRITE))==-1) printf ("smth wrong - can't open mm.dat for writing\n");
#else
if ((f=open(mmdatfile, O_WRONLY|O_CREAT|O_TRUNC|S_IREAD|S_IWRITE))==-1) printf ("smth wrong - can't open mm.dat for writing\n");
#endif
//j=sizeof(ftime);write (f,&j,4);write (f,&i1,4);write (f,&ft,sizeof(ftime)); // порядок: размер стуктуры, длина, дата
write (f,&(fst.st_mtime),sizeof(time_t));write (f,&(fst.st_size),sizeof(off_t)); // time, size

write (f,&TotalNumofMatr,4);
j=sizeof(MatrixElem);
printf ("MitrixElement size is %i\n",j);
write  (f,&j, sizeof j);
for (i=0;i<=TotalNumofMatr-1;i++)
{ write (f,(void*)ML[i],j);
  //delete(ML[i]);
  free (ML[i]); //linux
}
#ifdef Win32
#else
chmod (mmdatfile,438); // mode -rw-rw-rw linux
#endif
close (f);


FileReady://здесь читаем файл HashedMatrix.dat
fprintf (outerr,"Reading matrixfile\n");

// Для начала изучим что в профайле
if (!ProfileAllValues){
if (!(inf1=fopen(profile,"r"))) {fprintf (outerr,"Error opening profile\n");exit(0);}
fgets (s,MaxLineLength,inf1); fgets (s,MaxLineLength,inf1);fgets (s,MaxLineLength,inf1);fgets (s,MaxLineLength,inf1); // skip 1st 4 lines
i=0;
while (!(feof(inf1))&&(fgets (s,MaxLineLength,inf1))&&(strstr(s,"//")==NULL))
   {
    if (sscanf (s,"%f %f %f %s %s\n",&x1,&PL[i].CoreThr,&PL[i].Thr,PL[i].Num,PL[i].Name)!=5) {fprintf (outerr,"Bad profile line %s\n",s); exit(0);}
//printf ("%f %f %s %s\n",PL[i].CoreThr,PL[i].Thr,PL[i].Num,PL[i].Name);
    PL[i].Thr += ProfileTuneBy;
    i++;}
fclose(inf1);
}

#ifdef Win32
f=open(mmdatfile, O_RDONLY|O_BINARY);
#else
f=open(mmdatfile, O_RDONLY);
#endif
//read (f,&j,4);read (f,&j,4); read (f,&ft1,sizeof(ftime));
read (f,&FileDate,sizeof (time_t));read (f,&FileSize,sizeof(off_t));// linux
read (f,&TotalNumofMatr,4); read (f,&i1,4); if (i1!=sizeof(MatrixElem)) {fprintf (outerr,"MatrixElem size differ: here %d in file %d- rebuild mm.dat\n",sizeof(MatrixElem),i1);exit(0);}
i2=0; matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
if (!ProfileAllValues)
{for (i3=0;i3<=TotalNumofMatr-1;i3++)
    {read (f,matr,i1);
     for (i4=0;i4<=i-1;i4++)
         if (strcmp(PL[i4].Name,matr->Name)==0)
            {matr->CoreThr=(int)(matr->CoreMax*PL[i4].CoreThr);  // вот она нужная матрица
             matr->Thr=(int)(matr->Max*PL[i4].Thr);
             ML[i2++]=matr;
             matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
//             i4=i;
             goto ReadNextMatr;
            }
     ReadNextMatr: ;
    }
}
else // Читаем все матрицы и для всех одно значение порогов. Коровый порог = 0.75 от главного
{
for (i3=0;i3<=TotalNumofMatr-1;i3++)
    {read (f,matr,i1);
             matr->CoreThr=(int)(matr->CoreMax*ProfileAllValues*0.75);  // вот она нужная матрица
             matr->Thr=(int)(matr->Max*ProfileAllValues);
             ML[i2++]=matr;
             matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
    }
}
close(f);
NumberOfMatrix=i2;

// проверка матрицы  - выведем

/*matr=ML[0];
printf ("Debug : %d\n", matr->Length);
printf ("Shifts : %d<->%d\n", matr->Shift1,matr->Shift2);
printf ("LRsides : %d<->%d\n", matr->LeftSide,matr->RightSide);

for (i=0;i <= matr->Length-1;i++) printf ("%d-%d\t%d\t%d\t%d\t%d\n",i,matr->Val[0][i],matr->Val[1][i],matr->Val[2][i],matr->Val[3][i],matr->Val[4][i]);
for (i=0;i <= matr->LeftSide+matr->RightSide-1;i++)
     { for (i1=0;i1<=124;i1++) printf ("%d\t",matr->ValT[i1][i]);
       printf ( "\n");
     }
*/

return (i2);
}


int PrintSeq (char* St,char *End,char *CoreSt)
{
for (;St<=CoreSt-1;St++)
    switch (*St)
    {case 0 : fprintf (outf,"a");break;
     case 1 : fprintf (outf,"c");break;
     case 2 : fprintf (outf,"g");break;
     case 3 : fprintf (outf,"t");break;
     case 4 : fprintf (outf,"n");break;}

for (;St<=CoreSt+4;St++)
    switch (*St)
    {case 0 : fprintf (outf,"A");break;
     case 1 : fprintf (outf,"C");break;
     case 2 : fprintf (outf,"G");break;
     case 3 : fprintf (outf,"T");break;
     case 4 : fprintf (outf,"N");break;}

for (;St<=End;St++)
    switch (*St)
    {case 0 : fprintf (outf,"a");break;
     case 1 : fprintf (outf,"c");break;
     case 2 : fprintf (outf,"g");break;
     case 3 : fprintf (outf,"t");break;
     case 4 : fprintf (outf,"n");break;}
}


int AddSeq (char* St,char *End,char *CoreSt)
{
for (;St<=CoreSt-1;St++)
     OutStr[l++]=SeqChars[*St];

//St+=5;l+=5;
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];

for (;St<=End;St++)
     OutStr[l++]=SeqChars[*St];
}

int AddCoreSeq (char* St)
{
l=60;
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
OutStr[l++]=CoreChars[*(St++)];
}

int ReadOptionalParameters (int n, char **str)
{
int i=0;
Verbose=0; // by default
i=0;
while (++i <= n-1)
{       if (strstr(str[i],"/v")==str[i]) {if ((i+1 <= n-1)&&(str[i+1])) {Verbose=atoi(str[++i]);} else {fprintf(outerr,"Error in parameters - Verbose must have value. Exmpl: /v 3\n"); exit (0); } }
   else if (strstr(str[i],"/mdprf")==str[i]) {ProfileTuneBy=atof(str[++i]);}
   else if (strstr(str[i],"/prf")==str[i]) {ProfileAllValues=atof(str[++i]);}
   else if (strstr(str[i],"+")==str[i]) OnSeqOnly=0; else OnSeqOnly=1;
}

}


#include "OtherFunctions.h"

//*********************************************************************************************************//
//*********************************************************************************************************//
int main(int argc, char** argv)
{
inittime
#ifdef CGreport
printf ("GC content in window (based on mm) - 4 Feb 2011  Igor Deyneko (c)\nUse: mm SeqFile WinLength\n");
goto Param;
#endif
printf ("FastMatch v.2.5 - 13 Jul 2010  Igor Deyneko (c)\n");
#ifdef MCatchOutPut
printf ("***************Special edition for MCatch******************************\n");
#endif
#ifdef MCatchFreqTable
printf ("***************Special edition for MCatch******************************\n");
#endif
printf ("Output: ");
#ifdef MCatchOutPut
printf ("BinaryMcatch ");
#endif
#ifdef Binary
printf ("Binary ");
#endif
#ifdef Text
printf ("Text ");
#endif
printf ("\n");
fprintf (outerr,"Info: matrixElem size: %d\n",sizeof(MatrixElem) );

Param:
#ifdef CGreport
goto BendCalculations;

if (argc < 2) {fprintf(outerr,"Error in parameters\n"); exit(0);}
infile1=argv[1];
GCWindow=atoi(argv[2]);
GC();

fprintf (outerr,"All done\n");
exit(1);

BendCalculations:
if (argc < 2) {fprintf(outerr,"Error in parameters\n"); exit(0);}
infile1=argv[1];
GCWindow=atoi(argv[2]);
printf ("BendCalculations, GCWindow=%d\n",GCWindow);

Bend();
fprintf (outerr,"All done\n");
exit(1);

#endif

#ifdef Triplets
printf ("Countung triplets in DNA seqs\n Use: .exe seqfile\n");
if (argc < 2) {fprintf(outerr,"Error in parameters\n"); exit(0);}
infile1=argv[1];
CountTriplets();

fprintf (outerr,"All done\n");
exit(1);
#endif

#ifdef CountEntropyProfile
if (argc < 2) {fprintf(outerr,"Error in parameters\n"); exit(0);}
infile1=argv[1]; GCWindow=atoi(argv[2]);
EntropyProfile();
fprintf (outerr,"All done\n"); exit(1);
#endif

if (argc < 5)
//  if ((argc>=3)&&(strstr(argv[2],"RemakeDat")!=NULL))
  //   {MakeDatFile (MatrixList,argv[1]); exit(0);}
    // else
         {fprintf(outerr,"Error in parameters\n");
	       fprintf(outerr,"Usage: mm matrixlibfile seqfile outfile profile [/prf float] [/mdprf float] [+] [/v 0..3]\n");
	       fprintf(outerr," /prf - Use all PWM with these core and main threshods (profile in not required)\n /mdprf Val - modify profile = profile main thr + Val\n /v 0..3 - Verbose level\n + - output all sites, including those on sequence edges\n OR \n mm matrixlibfile RemakeDat - to rebuild mm.dat file - NOT yet implemented\n");
          exit(0);}
libfile=argv[1];
infile1=argv[2];
outfile=argv[3];

ReadOptionalParameters (argc,argv);
// check if relative path
mmdatfile[0]=0;
i=strlen (libfile); while (i-->=0) {if ((libfile[i]=='/')||(libfile[i]=='\\')) {strcpy (mmdatfile,libfile);mmdatfile[i+1]=0; i=-1;} }
strcpy (mmdatfile+strlen(mmdatfile),"mm.dat");
fprintf (outerr,"Path for mm.dat=%s\n",mmdatfile);

//i=open (infile1,O_RDONLY|O_TEXT); xl=filelength(i); close(i);
stat (infile1,&fst);
xl=fst.st_size;

if (!(inf1=fopen(infile1,"r"))) { fprintf (outerr,"Error opening file: %s\n",infile1); exit(0);}

Seq=(char*)malloc(xl+100)+50; // 100 for ends
if (Seq==NULL) {printf( "Memory problem\n");}

//for seq encoding
/*SeqLength=LoadSeq (inf1,Seq);
for (i=0;i<=1E7-1;)
{for (j=1;j<=80;j++) {printf ("%d ",Seq[i]); i++;}
//printf("\n");
}
exit (0);*/
// end


NumberOfMatrix=LoadMatrix (MatrixList,argv[1],argv[4]);
fprintf (outerr,"Number of matr=%d\n",NumberOfMatrix);
printtime

// Generate correspondence MatrixID - MatrixName
/*for (i=0;i<=NumberOfMatrix-1;i++)
{ printf ("%d %d %s %s\n",i,MatrixList[i]->ID,MatrixList[i]->Name,MatrixList[i]->Name+2);
}
exit (0);
  */

#ifdef EncodeSeq
goto ReversedSeq;
// прилада такая - выдает сиквенс как сигнал, а сигнар равен скору, отстой короче
SeqLength=LoadSeq (inf1,Seq);
M=MatrixList[0]; printf ("%s\n",M->Name);
for (i=0;i<=SeqLength-1;i++)
{ sum=0;
  for (j=M->Shift1;j<=M->Length-M->Shift2-1;j++)
                       sum+=M->Val[Seq[i+j]][j];
  printf ("%f\n",sum/((float)M->Max));
}
exit (0);

//Reversed sequence ...  - пишем обратную последовательность(и) и всё.
ReversedSeq:
if (!(outf=fopen(outfile,"w"))) { fprintf (outerr,"Error opening file: %s\n",outfile); exit(0);}
while (!feof(inf1))
{
// загружаем очередной сиквенс и вперед
if (!(SeqLength=LoadSeq (inf1,Seq))) {continue;} //return (0); // все, скорей всего конец файла
fprintf (outf,">%s\n",Name);
for (i=SeqLength-1;i>=0;i--)
  { if (Seq[i]==0) fprintf (outf,"T"); else {
    if (Seq[i]==1) fprintf (outf,"G"); else {
    if (Seq[i]==2) fprintf (outf,"C"); else {
    if (Seq[i]==3) fprintf (outf,"A"); else {fprintf (outf,"N");}}}}
    if (((SeqLength-i)%80 ==0)&&(SeqLength-i)) fprintf (outf,"\n");
  }
fprintf (outf,"\n");
}
fclose (outf);
exit (0);
//end
#endif

// бинарный файл для результата
#ifdef Binary
Info=(InfBlockP)malloc(sizeof(InfBlock));
strcpy (outfileb, outfile); strcpy (outfileb+strlen(outfileb),".bin");
BigBuffer=(short int*)malloc(BigBufferSize*sizeof(short int));
BigBufferEnd=0;
#ifdef Win32
if ((outbin=open(outfileb, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY,S_IREAD|S_IWRITE))==-1) printf ("smth wrong - can't open for writing OutFileBinary\n");
#else
if ((outbin=open(outfileb, O_WRONLY|O_CREAT|O_TRUNC,S_IREAD|S_IWRITE))==-1) printf ("smth wrong - can't open for writing OutFileBinary\n");
#endif
write (outbin,Info,sizeof(InfBlock)); // информационный блок
//write (outb,&i,4); // оставим место для количества сайтов
//запишем длину сиквенса
//write (outb,&i,4); // тоже оставим место - потом после поиска запишем
#endif
#ifdef MCatchOutPut
#ifdef Win32
MCatchoutb=open("sites.bin", O_WRONLY|O_CREAT|O_TRUNC|O_BINARY,S_IREAD|S_IWRITE);
#else
MCatchoutb=open("sites.bin", O_WRONLY|O_CREAT|O_TRUNC,S_IREAD|S_IWRITE);
#endif
//write (MCatchoutb,&i,4);write (MCatchoutb,&i,4); // и длину
Info=(InfBlockP)malloc(sizeof(InfBlock));
write (MCatchoutb,Info,sizeof(InfBlock)); // информационный блок
OnSeqOnly=1; // нам надо только НА послед
#endif

//write (outb,&NumberOfMatrix,4);
//for (i=0;i<=NumberOfMatrix-1;i++)
//{  write (outb,MatrixList[i]->Name,15);write (outb,&i,4);}


// проверка
/*printf ("%s %d %d Thr=%d/%d \n",MatrixList[0]->Name,MatrixList[0]->CoreMax,MatrixList[0]->Max,MatrixList[0]->CoreThr,MatrixList[0]->Thr);
printf ("%s %d %d Thr=%d/%d \n",MatrixList[1]->Name,MatrixList[1]->CoreMax,MatrixList[1]->Max,MatrixList[1]->CoreThr,MatrixList[1]->Thr);
printf ("%s %d %d Thr=%d/%d \n",MatrixList[2]->Name,MatrixList[2]->CoreMax,MatrixList[2]->Max,MatrixList[2]->CoreThr,MatrixList[2]->Thr);
M=MatrixList[0];
for (i=0; i<=M->Length-1;i++)
{printf ("%d\t%d\t%d\t%d\t%d\n",M->ValR[0][i],M->ValR[1][i],M->ValR[2][i],M->ValR[3][i],M->ValR[4][i]);
}
for (i=0; i<=3124;i++)
{printf ("%d - %d\n",M->HR[i].AccN,M->HR[i].Val);
} */

// ну и собственно поиск
//  теперь новый !!!
/*if (!(outf=fopen(outfile,"w"))) { fprintf (outerr,"Error opening file: %s\n",outfile); exit(0);}
fprintf (outf,"Search for sites by WeightMatrix library: %s\n",argv[1]);
fprintf (outf,"Sequence file: %s\n",argv[2]);
fprintf (outf,"Site selection profile: %s\n\n",argv[4]);
*/
#ifdef Text
// пишем теперь по новому - текстовый файл
#ifdef Win32
if ((outfb=open(outfile, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY,S_IREAD|S_IWRITE))==-1) { fprintf (outerr,"Error opening file: %s\n",outfile); exit(0);}
#else
if ((outfb=open(outfile, O_WRONLY|O_CREAT|O_TRUNC,S_IREAD|S_IWRITE))==-1) { fprintf (outerr,"Error opening file: %s\n",outfile); exit(0);}
#endif
strcpy(OutStr,"Search for sites by WeightMatrix library: "); strcpy (OutStr+42,argv[1]); write (outfb,OutStr,strlen(OutStr));
OutStr[0]=13;OutStr[1]=10;strcpy(OutStr+2,"Sequence file: "); strcpy (OutStr+17,argv[2]); write (outfb,OutStr,strlen(OutStr));
OutStr[0]=13;OutStr[1]=10;strcpy(OutStr+2,"Site selection profile: "); strcpy (OutStr+26,argv[4]); write (outfb,OutStr,strlen(OutStr));
write (outfb,OutStr,2);
//OutStr[2]=13;OutStr[3]=10; write (outfb,OutStr,4);
#endif

for (i=0;i<=NumberOfMatrix-1;i++) {TotalSiteStat[i]=0;} // глобальная статистика
for (j=0;j<=4;j++) GContent[j]=0;
for (j=0;j<=24;j++) DiContent[j]=0;

while (!feof(inf1))
{
// загружаем очередной сиквенс и вперед
if (!(SeqLength=LoadSeq (inf1,Seq))) {continue;} //return (0); // все, скорей всего конец файла
//if (!(SeqLength)) continue;
if (SeqLength+100 > MaxSeqLength ) {fprintf (outusr,"Your sequence too long\n"); return (0);}
if (SeqLength <= MinSeqLength-1) {fprintf (outusr,"your sequence too short"); return (0);}
//printf ("%d/%d\n",SeqLength,xl);
// проверим сиквенс
//   for (xl=-49;xl<=SeqLength+46;xl++) {printf ("%c",Seq[xl]+97);}
// сделать seqTriple
SeqT=(char*)malloc(SeqLength+100)+50; // 100 for ends
for (xl=-49;xl<=SeqLength+46;xl++)
    {SeqT[xl]=Seq[xl]+5*Seq[xl+1]+25*Seq[xl+2];}

if (Verbose>=2)
{ for (xl=0;xl<=SeqLength-1;xl++) GContent[Seq[xl]]++;
  for (xl=0;xl<=SeqLength-2;xl++) DiContent[Seq[xl]*5+Seq[xl+1]]++;
}

// проверим эту штуку
/*for (xl=0;xl<=100;xl++) printf ("%d ",(int)Seq[xl]); printf ("\n");
for (xl=0;xl<=100;xl++) printf ("%d ",(int)SeqT[xl]);
*/
printtime

PentaSeq=MakePentaSeq (Seq,SeqLength);
printtime

NumOfSites=0;
NumofSeq++;

#ifdef Text
//fprintf (outf,"Inspecting sequence ID   %s\n",Name);
memcpy(OutStr+2,"Inspecting sequence ID   ",25); strcpy (OutStr+27,Name); write (outfb,OutStr,strlen(OutStr)); write (outfb,OutStr,2);
strcpy(OutStr+2,"                       |          (+) |        |        | "); // такой вот шаблон, первые 2 байта - перенос строки!
#endif
for (i=0;i<=NumberOfMatrix-1;i++)
{  M=MatrixList[i];
   l=0; while (OutStr[l+++3]=M->Name[l]); for (l+=2;l<=24;l++) OutStr[l]=32;
   OutStr[37]='+';
   Strand=0;
   j=0;
   CurrentNumOfSites=0;  // сайты текущего типа
   while ((M->H[j].Val >= M->CoreThr)&&(M->H[j].Val+M->Max-M->CoreMax >= M->Thr)&&(j<=3124)) // если кор скор сильно маленький
   {   P=PentaSeq[M->H[j].AccN];
       sprintf (OutStr+43,"%5.3f",((float)M->H[j].Val)/M->CoreMax); OutStr[46+2]=' ';
//       l=M->CoreStart+60; AddCoreSeq (Seq+P->Pos);    // типа хотим кор сразу вывести а потом только края //  смотреть чтоб P!=Null
       while (P)
            {//sum=0;
             //for (l=0;l<=M->Length-1;l++) sum+=M->Val[Seq[P->Pos-M->CoreStart+l]][l];
             sum=M->H[j].Val;
//             if (sum+M->Min-M->CoreMin > M->Thr) fast++; else all++;
             xl=P->Pos-M->CoreStart-3;
             for (l=0;l<=M->LeftSide-1;l++) sum+=M->ValT[SeqT[(xl+=3)]][l];
             xl+=5;
             for (;l<=M->LeftSide+M->RightSide-1;l++) sum+=M->ValT[SeqT[(xl+=3)]][l];
             if (sum > M->Thr)
             if ((!OnSeqOnly)||((P->Pos-M->CoreStart+1+M->Shift1 > 0)&&(P->Pos-M->CoreStart+M->Length-M->Shift2 <= SeqLength)))
                {output:
                 CurrentNumOfSites++; //NumOfSites++; //TotalSitesLen+=M->Length;
// /*Text1*/       fprintf (outf,"%-10s (+) %d %5.4f %5.4f \n",M->Name,P->Pos-M->CoreStart+1+M->Shift1,((float)M->H[j].Val)/M->CoreMax,((float)sum)/M->Max);  // вот ОН!
// /*Text1*/       fprintf (outf," %-22s|%9d (+) |  %5.3f |  %5.3f | ",M->Name,P->Pos-M->CoreStart+1+M->Shift1,((float)M->H[j].Val)/M->CoreMax,((float)sum)/M->Max);  // формат матча
// /*Text1*/       PrintSeq (Seq+(P->Pos-M->CoreStart+M->Shift1),Seq+(P->Pos-M->CoreStart+M->Length-M->Shift2-1),Seq+P->Pos);
// /*Text1*/       fprintf (outf,"\n");
                 x=P->Pos-M->CoreStart+M->Shift1; // здесь выводим начало с нуля
#ifdef Binary
/*Binary out*/    //write(outbin,&M->ID,4);write(outbin,&x,4);write(outbin,&Strand,4);write(outbin,&sum,4); //write(outb,&P->Pos,4); // глобальный номер матрицы, позиция, начало кора (абс коорд).
/*New fast and compact)*/
                 //buffer[0]=M->ID; memcpy (buffer+1,&x,4); memcpy (buffer+3,&sum,4);
                 //write(outbin,buffer,10);
/*Noch schneller*/
                  BigBuffer[BigBufferEnd]=(M->ID+1);memcpy (BigBuffer+BigBufferEnd+1,&x,4); memcpy (BigBuffer+BigBufferEnd+3,&sum,4);BigBufferEnd+=5;  // M->ID+1 потому что обр послед записывается с минусом
                  if (BigBufferEnd >= BigBufferSize-20) {write(outbin,BigBuffer,BigBufferEnd*2); BigBufferEnd=0;}
#endif
     #ifdef MCatchOutPut
                 write(MCatchoutb,&M->ID,4);write(MCatchoutb,&x,4);write(MCatchoutb,&Strand,4);write(MCatchoutb,&sum,4); // глобальный номер матрицы, позиция, цепь, вес.
	  #endif
//                   fprintf (outf,"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n");
//                 sprintf (OutStr+52,"%5.3f",(((float)sum)/M->Max));
#ifdef Text
/*Text2*/        if (sum==M->Max) memcpy (OutStr+52,"1.000",5); else {OutStr[52]='0';OutStr[53]='.'; sprintf (OutStr+54,"%3d",(int)(1000*((sum)/(float)M->Max)));OutStr[57]=' ';} // !!sprintf copies 0 as end of string!!!
/*Text2*/        sprintf(OutStr+26,"%9d",P->Pos-M->CoreStart+1+M->Shift1); OutStr[35]=' ';
/*Text2*/        l=60;AddSeq(Seq+(P->Pos-M->CoreStart+M->Shift1),Seq+(P->Pos-M->CoreStart+M->Length-M->Shift2-1),Seq+P->Pos);
       sprintf (OutStr+43,"%5d",P->Pos-M->CoreStart+1+M->Shift1+M->OrigLength-1); OutStr[46+2]=' '; // до
/*Text2*/        write (outfb,OutStr,l);
#endif
#ifdef MCatchFreqTable
                  y=sum/(float)M->Max; write(outbin,&M->ID,4);write(outbin,&y,4); // глобальный номер матрицы, удельный скор.
#endif
                }
             P=P->Next;
            }
       j++;
   }
   j=0;
   OutStr[37]='-';
   Strand=1;
   while ((M->HR[j].Val >= M->CoreThr)&&(M->HR[j].Val+M->Max-M->CoreMax >= M->Thr)&&(j<=3124))
   {   P=PentaSeq[M->HR[j].AccN];
       sprintf (OutStr+43,"%5.3f",((float)M->HR[j].Val)/M->CoreMax);  OutStr[46+2]=' ';
       while (P)
            {//sum=0;
             //for (l=0;l<=M->Length-1;l++) sum+=M->ValR[Seq[P->Pos-M->CoreStartR+l]][l];
             sum=M->HR[j].Val;
             xl=P->Pos-M->CoreStartR-3;
             for (l=0;l<=M->RightSide-1;l++) sum+=M->ValTR[SeqT[(xl+=3)]][l];
             xl+=5;
             for (;l<=M->RightSide+M->LeftSide-1;l++) sum+=M->ValTR[SeqT[(xl+=3)]][l];
             if (sum > M->Thr)
             if (!(OnSeqOnly)||((P->Pos-M->CoreStartR+1+M->Shift2)>0)&&(P->Pos-M->CoreStartR+M->Length-M->Shift1 <= SeqLength))
                {
                 CurrentNumOfSites++;//NumOfSites++; //TotalSitesLen+=M->Length;
// /*Text1*/       fprintf (outf,"%-10s (-) %d %5.4f %5.4f \n",M->Name,P->Pos-M->CoreStartR+1+M->Shift2,((float)M->HR[j].Val)/M->CoreMax,((float)sum)/M->Max);  // вот ОН!
// /*Text1*/       fprintf (outf," %-22s|%9d (-) |  %5.3f |  %5.3f | ",M->Name,P->Pos-M->CoreStartR+1+M->Shift2,((float)M->HR[j].Val)/M->CoreMax,((float)sum)/M->Max);  // формат матча
// /*Text1*/       PrintSeqR (Seq+(P->Pos-M->CoreStart+M->Length-M->Shift2-1),Seq+(P->Pos-M->CoreStart+M->Shift1),Seq+(P->Pos-2*M->CoreStart+M->Shift1+M->Length-M->Shift2-1)); // обратная последовательность
// /*Text1*/       PrintSeq (Seq+(P->Pos-M->CoreStartR+M->Shift2),Seq+(P->Pos-M->CoreStartR+M->Length-M->Shift1-1),Seq+P->Pos);
// /*Text1*/       fprintf (outf,"\n");

                 x=P->Pos-M->CoreStartR+M->Shift2;
#ifdef Binary
/*Binary out*/   //write (outbin,&M->ID,4);write (outbin,&x,4);write(outbin,&Strand,4);write(outbin,&sum,4); // write(outb,&P->Pos,4);
/*New fast and compact)*/
                 //buffer[0]=-1*(M->ID); memcpy (buffer+1,&x,4); memcpy (buffer+3,&sum,4);
                 //write(outbin,buffer,10);
/*Noch schneller*/
                  BigBuffer[BigBufferEnd]=-1*(M->ID+1);memcpy (BigBuffer+BigBufferEnd+1,&x,4); memcpy (BigBuffer+BigBufferEnd+3,&sum,4);BigBufferEnd+=5; // class position score
                  if (BigBufferEnd >= BigBufferSize-20) {write(outbin,BigBuffer,BigBufferEnd*2); BigBufferEnd=0;}
#endif
       #ifdef MCatchOutPut
                 write(MCatchoutb,&M->ID,4);write(MCatchoutb,&x,4);write(MCatchoutb,&Strand,4);write(MCatchoutb,&sum,4); // глобальный номер матрицы, позиция, цепь, вес.
		 #endif
//                   fprintf (outf,"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n");
//                   write (outb,Seq,80);
//                 sprintf (OutStr+52,"%5.3f",(((float)sum)/M->Max));
#ifdef Text
/*Text2*/        if (sum==M->Max) memcpy (OutStr+52,"1.000",5); else {OutStr[52]='0';OutStr[53]='.'; sprintf (OutStr+54,"%3d",(int)(1000*((sum)/(float)M->Max))); OutStr[57]=' ';}
/*Text2*/        sprintf(OutStr+26,"%9d",P->Pos-M->CoreStartR+1+M->Shift2);  OutStr[35]=' ';
/*Text2*/        l=60;AddSeq(Seq+(P->Pos-M->CoreStartR+M->Shift2),Seq+(P->Pos-M->CoreStartR+M->Length-M->Shift1-1),Seq+P->Pos);
       sprintf (OutStr+43,"%5d",P->Pos-M->CoreStartR+1+M->Shift2+M->OrigLength-1); OutStr[46+2]=' '; // до
/*Text2*/        write (outfb,OutStr,l);
#endif
#ifdef MCatchFreqTable
                  y=sum/(float)M->Max; write(outbin,&M->ID,4);write(outbin,&y,4); // глобальный номер матрицы, скор.
#endif
                }
             P=P->Next;
            }
       j++;
   } // след пятерка
if (Verbose>=3) printf ("%d %.30s\t- %d\n",i,M->Name,CurrentNumOfSites);
NumOfSites+=CurrentNumOfSites;
TotalSiteStat[i]+=CurrentNumOfSites;
}  //след матрица
printtime


TotalSeqLength+=SeqLength;
TotalNumOfSites+=NumOfSites;

if (Verbose>=2) {
fprintf (outerr,"NumOfSites=%-10d SeqLength=%d\n",NumOfSites,SeqLength);
//fprintf (outerr,"SeqLength=%d\n",SeqLength);
}

#if defined(Binary) || defined (MCatchOutPut)
Info->SeqLength[NumofSeq-1]=SeqLength;
Info->NumofSites[NumofSeq-1]=NumOfSites;
#endif


free (SeqT-50);
free (PentaSeq[Seq[0]+5*(Seq[1]+5*(Seq[2]+5*(Seq[3]+5*Seq[4])))]);
free (PentaSeq);
//printf ("Next seq\n");
// Здесь можно переходить к след послед.
}

if (Verbose>=1)
  { printf ("Global statistics (sites per 1Kb):\n");
    for (i=0;i<=NumberOfMatrix-1;i++) printf ("%3d %-20s\t- %d\t%f\n",i,MatrixList[i]->Name,TotalSiteStat[i],1000*TotalSiteStat[i]/(float)TotalSeqLength);
  }

if (Verbose>=2)
{ xl=GContent[0]+GContent[1]+GContent[2]+GContent[3];
  printf ("Nucleotide freq:\nA: %4.3f\nC: %4.3f\nG: %4.3f\nT: %4.3f\n",GContent[0]/(float)xl,GContent[1]/(float)xl,GContent[2]/(float)xl,GContent[3]/(float)xl);
//DiNucleotide content
  xl=DiContent[0]+DiContent[1]+DiContent[2]+DiContent[3]+DiContent[5]+DiContent[6]+DiContent[7]+DiContent[8]+DiContent[10]+DiContent[11]+DiContent[12]+DiContent[13]+DiContent[15]+DiContent[16]+DiContent[17]+DiContent[18];
  printf ("DiNucleotide freq:\nAA: %4.3f\tAC: %4.3f\tAG: %4.3f\tAT: %4.3f\n",DiContent[0]/(float)xl,DiContent[1]/(float)xl,DiContent[2]/(float)xl,DiContent[3]/(float)xl);
  printf ("CA: %4.3f\tCC: %4.3f\tCG: %4.3f\tCT: %4.3f\n",DiContent[5]/(float)xl,DiContent[6]/(float)xl,DiContent[7]/(float)xl,DiContent[8]/(float)xl);
  printf ("GA: %4.3f\tGC: %4.3f\tGG: %4.3f\tGT: %4.3f\n",DiContent[10]/(float)xl,DiContent[11]/(float)xl,DiContent[12]/(float)xl,DiContent[13]/(float)xl);
  printf ("TA: %4.3f\tTC: %4.3f\tTG: %4.3f\tTT: %4.3f\n",DiContent[15]/(float)xl,DiContent[16]/(float)xl,DiContent[17]/(float)xl,DiContent[18]/(float)xl);
}

#ifdef Binary
Info->NumofSeq=NumofSeq;
//Fluch BigBuffer
if (BigBufferEnd > 0) {write(outbin,BigBuffer,BigBufferEnd*2); BigBufferEnd=0;}
close(outbin);
#endif
#ifdef MCatchOutPut
Info->NumofSeq=NumofSeq;
close (MCatchoutb);
#endif

#ifdef Binary
// а теперь запишем наш блок     //количество сайтов в файл
#ifdef Win32
if ((outbin=open(outfileb, O_WRONLY|O_BINARY,S_IREAD|S_IWRITE))==-1) printf ("smth wrong - can't open for writing OutFileBinary, InfoBlock\n");
#else
if ((outbin=open(outfileb, O_WRONLY,S_IREAD|S_IWRITE))==-1) printf ("smth wrong -  - can't open for writing OutFileBinary, InfoBlock\n");
chmod (outfileb,438); // mode -rw-rw-rw linux
#endif
write (outbin,Info,sizeof (InfBlock));
//write (outb,&NumOfSites,4);
//write (outb,&TotalSeqLength,4); // и длину
close (outbin);
#endif
#ifdef MCatchOutPut
#ifdef Win32
MCatchoutb=open("sites.bin", O_WRONLY|O_BINARY,S_IREAD|S_IWRITE);
#else
MCatchoutb=open("sites.bin", O_WRONLY,S_IREAD|S_IWRITE);
#endif
//write (MCatchoutb,&NumOfSites,4); write (MCatchoutb,&TotalSeqLength,4); // и длину
write (MCatchoutb,Info,sizeof (InfBlock));
close (MCatchoutb);
#endif

#ifdef Text
//printf ("%d/%d(%f)\n",fast,fast+all,((float)fast)/(fast+all));
//fprintf (outf,"\n Total sequences length=%d\n\n Total number of found sites=%d\n\n Frequency of sites per nucleotide=%f\n",TotalSeqLength,TotalNumOfSites,TotalNumOfSites/((float)TotalSeqLength));
write (outfb,OutStr,2);
sprintf (OutStr+2," Total sequences length=%d",TotalSeqLength); write (outfb,OutStr,strlen(OutStr));write (outfb,OutStr,2);
sprintf (OutStr+2," Total number of found sites=%d",TotalNumOfSites);write (outfb,OutStr,strlen(OutStr));write (outfb,OutStr,2);
sprintf (OutStr+2," Frequency of sites per nucleotide=%f",TotalNumOfSites/((float)TotalSeqLength));write (outfb,OutStr,strlen(OutStr));write (outfb,OutStr,2); write (outfb,OutStr,2);
close (outfb);
#endif

fprintf (outerr,"TotalNumOfSites=%d\n",TotalNumOfSites);
fprintf (outerr,"TotalSeqLength=%d\n",TotalSeqLength);
fprintf (outerr,"TotalNumofSeq=%d\n",NumofSeq);
//fclose(outf);

fprintf (outerr,"All done\n");
exit(1);
}

