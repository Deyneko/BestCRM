//функции дл€ MCtools, потом расширю на все программы.


int ReadOptionalParameters (int n, char **str)
{
int i=0;
Verbose=0; // by default
i=0;
//CoverDiff=2.0;
//FreqDiff=2.0;
//MaxCoverN=0.5;
//MinCoverY=0.75;
FullOutput=0;

while (++i <= n-1)
{       if (strstr(str[i],"/v")==str[i]) {Verbose=atoi(str[++i]);}
   else if (strstr(str[i],"MinCoverY")==str[i]) {MinCoverY=atof(str[++i]);}
   else if (strstr(str[i],"MaxCoverN")==str[i]) {MaxCoverN=atof(str[++i]);}
   else if (strstr(str[i],"FreqDiff")==str[i]) {FreqDiff=atof(str[++i]);}
   else if (strstr(str[i],"/mkprf")==str[i]) {MakeProfile=1;FreqPer1Kb=atof(str[++i]);}
   else if (strstr(str[i],"/FullOutput")==str[i]) {FullOutput=1;i++;}
   else if (strstr(str[i],"/window")==str[i]) {WindowLength=atoi(str[++i]);}

   }

}


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

PentaP* MakePentaSeq (char * seq, int N) // возвращ указатель на список указ на п€терки
{
PentaP Last,List,*PS,*EndsList;
int i,addr;
int j;
PS=(PentaP*)malloc(3125*sizeof(PentaP));
List=(PentaP)malloc(N*sizeof(Penta));
EndsList=(PentaP*)malloc(3125*sizeof(PentaP));
for (i=0;i<=3125-1;i++) {PS[i]=NULL;EndsList[i]=NULL;}

//for (j=0;j<=99;j++) // ƒл€ проверки
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

int LoadAllMatrix (MatrixElemP* ML, char* mmdatfile)
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

struct stat fst; //
 off_t FileSize;
 time_t FileDate;

//check if HahedMatrix.dat exist and have proper size and date
fprintf (outerr,"Did not checked if .lib == mm.dat, use mm.dat from %s\n",mmdatfile);
//stat(mmdatfile,&fst);

//if ((f=open(mmdatfile, O_RDONLY|O_BINARY))==-1) {fprintf (outerr,"open error\n"); goto buildfile;}
//read (f,&FileDate,sizeof(time_t)); if (FileDate!=fst.st_mtime) {fprintf (outerr,"FileDate error\n"); close (f); exit (0);}
//read (f,&FileSize,sizeof(off_t));  if (FileSize!=fst.st_size)  {fprintf (outerr,"FileSize error\n"); close (f); exit (0);}  //  date, length
//close(f);

FileReady://здесь читаем файл HashedMatrix.dat
fprintf (outerr,"Reading matrixfile\n");

#ifdef Win32
if ((f=open(mmdatfile, O_RDONLY|O_BINARY))==-1) {fprintf (outerr,"open error, mmdatfile=%s\n",mmdatfile);exit(1);}
#else
if ((f=open(mmdatfile, O_RDONLY))==-1) {fprintf (outerr,"open error, mmdatfile=%s\n",mmdatfile);exit(1);}
#endif
read (f,&FileDate,sizeof (time_t));read (f,&FileSize,sizeof(off_t));
read (f,&TotalNumofMatr,4); read (f,&i1,4); if (i1!=sizeof(MatrixElem)) {fprintf (outerr,"MatrixElem size differ: here %d in file %d - rebuild mm.dat\n",sizeof(MatrixElem),i1);exit(0);}
i2=0; matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
for (i3=0;i3<=TotalNumofMatr-1;i3++)
    {read (f,matr,i1);
//     for (i4=0;i4<=i-1;i4++)
//         if (strcmp(PL[i4].Name,matr->Name)==0)
//            {matr->CoreThr=(int)(matr->CoreMax*PL[i4].CoreThr);  // вот она нужна€ матрица
//             matr->Thr=(int)(matr->Max*PL[i4].Thr);
             ML[i2++]=matr;
             matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
//             i4=i;
//             goto ReadNextMatr;
//            }
//     ReadNextMatr:
    }

close(f);
NumberOfMatrix=i2;

return (i2);
}

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

struct stat fst; //
 off_t FileSize;
 time_t FileDate;

//check if HahedMatrix.dat exist and have proper size and date
stat(libfile,&fst);

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
    if (strstr(s,"CORE_START")!=NULL) {matr->CoreStart=atoi(s+10)-1;} // стартуем от нул€
    if (strstr(s,"CORE_LENGTH")!=NULL) {if(atoi(s+11)!=5) {fprintf (outerr,"Different CoreLength - %d\n",atoi(s+11));exit(0);}}
    if (strstr(s,"WEIGHTS")!=NULL)
    	  for (j=0;j<=matr->Length-1;j++)
        		{fgets (s,MaxLineLength,inf1);
             if (sscanf (s,"%i A:%f C:%f G:%f T:%f\n",&k,&x1,&x2,&x3,&x4)!=5) fprintf (outerr,"Bad matrix %s\n",matr->Name);
             *(matr->Val[0]+j)=(int)(MultFactor*x1);
             *(matr->Val[1]+j)=(int)(MultFactor*x2);
             *(matr->Val[2]+j)=(int)(MultFactor*x3);
             *(matr->Val[3]+j)=(int)(MultFactor*x4);
             *(matr->ValR[3]+matr->Length-1-j)=(int)(MultFactor*x1); // обратна€ матрица
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
         do   //  и дл€ второй цепи
         { i1=0;
           for (j=0;j<=3123;j++)
             if (matr->HR[j].Val < matr->HR[j+1].Val)
                {i2=matr->HR[j].AccN; i3=matr->HR[j].Val;
                    matr->HR[j].AccN=matr->HR[j+1].AccN; matr->HR[j].Val=matr->HR[j+1].Val;
                    matr->HR[j+1].AccN=i2; matr->HR[j+1].Val=i3;
                 i1=1;}
         }
         while (i1);
         matr->CoreMax=matr->H[0].Val;    // min max  дл€ кора
         matr->CoreMin=matr->H[3124].Val;

         // переделаем ћакс и ћин, (из-за вс€ких округлений надо)
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
// “еперь надо все записать
/*f=open (libfile,O_RDONLY|O_TEXT); // возмем информацию о файле матриц
   i1=filelength(f);
   getftime (f,&ft);
close(f);*/
stat (libfile,&fst);
#ifdef Win32
if ((f=open(mmdatfile, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY,S_IREAD|S_IWRITE))==-1) printf ("smth wrong - can't open mm.dat for writing\n");
#else
if ((f=open(mmdatfile, O_WRONLY|O_CREAT|O_TRUNC|S_IREAD|S_IWRITE,0666))==-1) printf ("smth wrong - can't open mm.dat for writing\n");
#endif
//j=sizeof(ftime);write (f,&j,4);write (f,&i1,4);write (f,&ft,sizeof(ftime)); // пор€док: размер стуктуры, длина, дата
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

// ƒл€ начала изучим что в профайле
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
read (f,&TotalNumofMatr,4); read (f,&i1,4); if (i1!=sizeof(MatrixElem)) {fprintf (outerr,"MatrixElem size differ: here %d in file %d - rebuild mm.dat\n",sizeof(MatrixElem),i1);exit(0);}
i2=0; matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
if (!ProfileAllValues)
{for (i3=0;i3<=TotalNumofMatr-1;i3++)
    {read (f,matr,i1);
     for (i4=0;i4<=i-1;i4++)
         if (strcmp(PL[i4].Name,matr->Name)==0)
            {matr->CoreThr=(int)(matr->CoreMax*PL[i4].CoreThr);  // вот она нужна€ матрица
             matr->Thr=(int)(matr->Max*PL[i4].Thr);
             ML[i2++]=matr;
             matr= (MatrixElemP)malloc(sizeof(MatrixElem)); //linux matr=new MatrixElem;
//             i4=i;
             goto ReadNextMatr;
            }
     ReadNextMatr: ;
    }
}
else // „итаем все матрицы и дл€ всех одно значение порогов.  оровый порог = 0.75 от главного
{
for (i3=0;i3<=TotalNumofMatr-1;i3++)
    {read (f,matr,i1);
             matr->CoreThr=(int)(matr->CoreMax*ProfileAllValues*0.75);  // вот она нужна€ матрица
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

SiteElemP* SortSitesbyBegin (SiteElemP* SL,int N,int L) // сортирует сайты по началу, возвр тот же список  но сортированный. пам€ть дл€ результата уже должна быть
{                                           // N - number of sites L - seq Length
SiteElemP SP;
int changed,i,j,ElongFactor=5;   //сортируем список сайтов по позизи€м  - круто сортируем!!!
//printf ("Sorting sites....");
SiteElemP* hashlist=(SiteElemP*)malloc(sizeof(SiteElemP)*(N+1000)*ElongFactor); for (i=0;i<=(ElongFactor*N+1000)-1;i++) hashlist[i]=NULL;
for (i=0;i<=N-1;i++)
    {j=floor(((float)SL[i]->Begin*(N*ElongFactor))/(float)L);
     while ((hashlist[j])) j++;
     hashlist[j]=SL[i]; }
//j=0; for (i=0;i<=N-1;)  if (hashlist[j]) SL[i++]=hashlist[j++]; else j++; // собираем все обратно
j=0; for (i=0;i<=N-1;)  if (hashlist[j]) SL[i++]=hashlist[j++]; else j++; // собираем все обратно
ElongFactor=0;
do                         // ну и надо довести дело до конца
{ changed=0; ElongFactor++;
  for (i=0;i<=N-2;i++)
       if (SL[i]->Begin > SL[i+1]->Begin)
           {SP=SL[i];
            SL[i]=SL[i+1];
            SL[i+1]=SP;
            changed=1;}
  for (i=N-2;i>=0;i--)
       if (SL[i]->Begin > SL[i+1]->Begin)
           {SP=SL[i];SL[i]=SL[i+1];SL[i+1]=SP;changed=1;}
}while (changed);
free (hashlist);
//printf ("  positions sorted, cycles=%d ",ElongFactor);
//printtime
return (SL);
}


int UniformPoints(int i1,int i2,int i3,int i4,int BestI1,int BestI2,int BestI3,int BestI4)
{
if ((i1>=BestI1)&&(i2>=BestI2)&&(i3>=BestI3)&&(i4>=BestI4)) return (0);
if (fmax(i1,fmax(i2,fmax(i3,i4))) > fmax(BestI1,fmax(BestI2,fmax(BestI3,BestI4)))) return (0);
//if (abs(i1-i2)+abs(i1-i3)+abs(i1-i4)+abs(i2-i3)+abs(i2-i4)+abs(i3-i4) <
  //  abs(BestI1-BestI2)+abs(BestI1-BestI3)+abs(BestI1-BestI4)+abs(BestI2-BestI3)+abs(BestI2-BestI4)+abs(BestI3-BestI4)) return (1); else return (0);
return (1);
}

