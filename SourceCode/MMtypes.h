//типы дл€ MCtools, потом расширю на все программы.

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
     int ValT[125][(int)MaxMatrLength/3+1],ValTR[125][(int)MaxMatrLength/3+1]; // «начени€ троек
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


typedef struct SiteElem_
	{int ClassID,Begin,Strand,Score; 
   } SiteElem, *SiteElemP;

typedef struct CEModelelem1_ // модель  Ё
	{ int AC;
//     int Beg, End,LBeg,LEnd,RBeg,REnd;
     int MatrIDL,MatrIDR; // номера матриц в  .dat
     int StrandL,StrandR; // ориентаци€ матриц
     int Min, Max; //  рассто€ние между началами ќригинальное
     int MinRelaxed, MaxRelaxed; //  рассто€ние между началами расслабленное
     int ScoreL,ScoreR; //
     int ScoreLRelaxed,ScoreRRelaxed; //
     int Gr;
     float ScoreC, Avrg;

//     char ID[50]; // тип
     int TotalNum;
   } CEModelelem1, * CEModelelem1P;

typedef struct CE_ProfileElem_
   { float LeftRelax, RightRelax, LengthRelax, CompositeRelax;
     int ID;
   } CE_ProfileElem, *CE_ProfileElemP;

typedef struct CETableElem_ {
 int CETable[PwmThrNum][PwmThrNum][LenThrNum][CSThrNum];
 } CETableElem, * CETableElemP;

typedef struct ModuleElem_ // модель  Ё
	{ int AC;
//     int Beg, End,LBeg,LEnd,RBeg,REnd;
     int MatrID1,MatrID2; // номера матриц в  .dat
     int Strand1,Strand2; // ориентаци€ матриц
     int Min,Max; //  рассто€ние между началами ќригинальное
     int MinRelaxed, MaxRelaxed; //  рассто€ние между началами расслабленное
     int Score1,Score2; //
     int Score1Relaxed,Score2Relaxed; //
//     int Gr;
//     float ScoreC, Avrg;
//     char ID[50]; // тип
     int TotalNum;
     char MatrName1[MaxMatrixNameLen], MatrName2[MaxMatrixNameLen];
   } ModuleElem, * ModuleElemP;
