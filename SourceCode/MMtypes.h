//���� ��� MCtools, ����� ������� �� ��� ���������.

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
     int CoreStart, CoreStartR; // ������� ��������� �� ���� �� ����� �������
     int LeftSide, RightSide; // ���-�� ����� �� ��������
     int Shift1, Shift2; // �������������� ���� � ����������� ������� �����, ������
     int ValT[125][(int)MaxMatrLength/3+1],ValTR[125][(int)MaxMatrLength/3+1]; // �������� �����
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

typedef struct CEModelelem1_ // ������ ��
	{ int AC;
//     int Beg, End,LBeg,LEnd,RBeg,REnd;
     int MatrIDL,MatrIDR; // ������ ������ �  .dat
     int StrandL,StrandR; // ���������� ������
     int Min, Max; //  ���������� ����� �������� ������������
     int MinRelaxed, MaxRelaxed; //  ���������� ����� �������� �������������
     int ScoreL,ScoreR; //
     int ScoreLRelaxed,ScoreRRelaxed; //
     int Gr;
     float ScoreC, Avrg;

//     char ID[50]; // ���
     int TotalNum;
   } CEModelelem1, * CEModelelem1P;

typedef struct CE_ProfileElem_
   { float LeftRelax, RightRelax, LengthRelax, CompositeRelax;
     int ID;
   } CE_ProfileElem, *CE_ProfileElemP;

typedef struct CETableElem_ {
 int CETable[PwmThrNum][PwmThrNum][LenThrNum][CSThrNum];
 } CETableElem, * CETableElemP;

typedef struct ModuleElem_ // ������ ��
	{ int AC;
//     int Beg, End,LBeg,LEnd,RBeg,REnd;
     int MatrID1,MatrID2; // ������ ������ �  .dat
     int Strand1,Strand2; // ���������� ������
     int Min,Max; //  ���������� ����� �������� ������������
     int MinRelaxed, MaxRelaxed; //  ���������� ����� �������� �������������
     int Score1,Score2; //
     int Score1Relaxed,Score2Relaxed; //
//     int Gr;
//     float ScoreC, Avrg;
//     char ID[50]; // ���
     int TotalNum;
     char MatrName1[MaxMatrixNameLen], MatrName2[MaxMatrixNameLen];
   } ModuleElem, * ModuleElemP;
