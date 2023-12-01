int GC()

{
int sum2;

SeqSet=(char**)malloc(MaxNumofSeq*sizeof(char*));
SeqLen=(int*)malloc(MaxNumofSeq*sizeof(int));
//float  GCAvrg[5], GCStDev[5], *GCDistr[5];

//read seq
NumofSeq=0;
stat (infile1,&fst);
xl=fst.st_size;
if (!(inf1=fopen(infile1,"r"))) { fprintf (outerr,"Error opening file: %s\n",infile1); exit(0);}
Seq=(char*)malloc(int(xl+100)+50);  if (Seq==NULL) {printf( "Memory problem\n");}
while (!feof(inf1))
{if (!(SeqLength=LoadSeq (inf1,Seq))) {continue;} //return (0); // ўбҐ, бЄ®аҐ© ўбҐЈ® Є®­Ґж д ©« 
 SeqSet[NumofSeq]=Seq;
 SeqLen[NumofSeq]=SeqLength;
 Seq+=SeqLength;
 SeqNames[NumofSeq]=(char*)malloc(MaxLineLength);
 strcpy(SeqNames[NumofSeq],Name);
 NumofSeq++;
}

GCDistr=(float**)malloc(5*sizeof(float*));
for (j=0;j<=4;j++) GCDistr[j]=(float*)malloc(NumofSeq*sizeof(float));

for (i=0;i<=NumofSeq-1;i++)
{ SeqLength=SeqLen[i];
  Seq=SeqSet[i];

  for (xl=0;xl<=SeqLength-1;xl++) GCDistr[Seq[xl]][i]+=1/(float)SeqLength; //for (xl=0;xl<=SeqLength-2;xl++) DiDistr[i][Seq[xl]*5+Seq[xl+1]]+=1;
  GCDistr[4][i]=GCDistr[1][i]+GCDistr[2][i]; //g+c
//  for (j=0;j<=4;j++) GCDistr[i][j]/=;
}
//for (i=0;i<=NumofSeq-1;i++) printf ("%f %f %f %f %f\n",GCDistr[0][i],GCDistr[1][i],GCDistr[2][i],GCDistr[3][i],GCDistr[4][i]); exit(0);

for (j=0;j<=4;j++) GCAvrg[j]=0; //for (j=0;j<=24;j++) DiDistr[i][j]=0;
for (i=0;i<=NumofSeq-1;i++) { for (j=0;j<=4;j++) GCAvrg[j]+=GCDistr[j][i]/(float)NumofSeq; }
for (i=0;i<=NumofSeq-1;i++) { for (j=0;j<=4;j++) GCStDev[j]+=(GCDistr[j][i]-GCAvrg[j])*(GCDistr[j][i]-GCAvrg[j])/(float)NumofSeq; }


//for (j=0;j<=4;j++)

printf ("Mean StDev:\n");
printf ("A:  %4.3f %f\n", GCAvrg[0], sqrt(GCStDev[0]));
printf ("C:  %4.3f %f\n", GCAvrg[1], sqrt(GCStDev[1]));
printf ("G:  %4.3f %f\n", GCAvrg[2], sqrt(GCStDev[2]));
printf ("T:  %4.3f %f\n", GCAvrg[3], sqrt(GCStDev[3]));
printf ("GC: %4.3f %f\n", GCAvrg[4], sqrt(GCStDev[4]));

// GC in Window
/*
for (i=0;i<=NumofSeq-1;i++)
{ SeqLength=SeqLen[i];
  Seq=SeqSet[i];

//printf ("%d\n",SeqLength-GCWindow);
for (xl=0;xl<=SeqLength-GCWindow-1;xl++)
{ for (j=0;j<=4;j++) GCAvrg[j]=0;
  for (l=0;l<=GCWindow-1;l++) GCAvrg[Seq[xl+l]]+=1/(float)GCWindow;
  printf ("%d  %4.3f\n",xl,  GCAvrg[1]+GCAvrg[2]);
}

} */

// GC in Window with Threshould
GCThr=0.65; // 10.488+4*0.028961; //0.5785; //0.65; //0.5785;
ATThr=0.30; // 0.488-4*0.028961; // 0.4001;  //0.65; //0.5785;
sum2=sum=0;
printf ("GCThr=%f, ATThr=%f\n",GCThr, ATThr);

for (i=0;i<=NumofSeq-1;i++)
{ SeqLength=SeqLen[i];
  Seq=SeqSet[i];

//printf ("%d\n",SeqLength-GCWindow);
for (xl=0;xl<=SeqLength-GCWindow-1;xl++)
{ for (j=0;j<=4;j++) GCAvrg[j]=0;
  for (l=0;l<=GCWindow-1;l++) GCAvrg[Seq[xl+l]]+=1/(float)GCWindow;
  if (GCAvrg[1]+GCAvrg[2] >=GCThr) Seq[xl]=6;
  if (GCAvrg[1]+GCAvrg[2] <=ATThr) Seq[xl]=7;
}

printf ("Inspecting sequence ID   %s\n\n",SeqNames[i]);
// A8                    |      510 (+) |    517 |  0.875 | AAAAAaac
j=0;
x=0;
for (xl=0;xl<=SeqLength-GCWindow-1;xl++)
 {if ((Seq[xl]==6)&&(j==0)) j=xl+1; //Begin
  if ((Seq[xl]!=6)&&(j>0)) {printf (" GCHigh                |   %6d (+) | %6d | gcgcgc\n", j,xl-1+1+GCWindow); j=0;xl+=GCWindow;x=1; } //end
 }
if ((j>0)) {printf (" GCHigh                |   %6d (+) | %6d | gcgcgc\n", j,xl-1+1+GCWindow); x=1; j=0;} //end
if (x==1) sum++;

//GCLow
j=0;
x=0;
for (xl=0;xl<=SeqLength-GCWindow-1;xl++)
 {if ((Seq[xl]==7)&&(j==0)) j=xl+1; //Begin
  if ((Seq[xl]!=7)&&(j>0)) {printf (" ATHigh                |   %6d (+) | %6d | atatat\n", j,xl-1+1+GCWindow); j=0;xl+=GCWindow;x=1; } //end
 }
if ((j>0)) {printf (" ATHigh                |   %6d (+) | %6d | atatat\n", j,xl-1+1+GCWindow); x=1; j=0;} //end
if (x==1) sum2++;

} //next seq

printf ("\n\nTotal seqs with GCHigh = %d,  %4.3f\n",sum, sum/(float)NumofSeq);
printf ("Total seqs with ATHigh = %d,  %4.3f\n",sum2, sum2/(float)NumofSeq);
}


int PrintNucleotide(int n)
{
    switch (n)
    {case 0 : printf ("a");break;
     case 1 : printf ("c");break;
     case 2 : printf ("g");break;
     case 3 : printf ("t");break;
     case 4 : printf ("n");break;}


}



int CountTriplets()
{
int TotalSeqTLength=0, TriplFreq[5*5*5],TriplTemp[5*5*5],TriplCover[5*5*5],DiFreq[5*5],MonoFreq[5],i1,i2,i3;
float Class1Weights[64]={1.265202559,0.990660649,1.000849317,1.152218386,1.156139106,0.899125899,0.889661319,1.046570973,1.064894141,0.826022811,1.044341722,1.083770855,1.218208445,1.114123837,0.898197183,1.270693385,1.006906864,1.029893493,0.893154559,1.220579132,0.869409414,0.874439834,0.866389928,0.796079075,0.844336175,0.940530279,0.987339415,0.952330195,0.785723924,1.12693464,0.856260528,0.989177489,0.87704003,0.857971687,1.052839876,1.113559501,0.978766647,0.777412713,0.908002559,0.868590081,0.909261261,0.885579687,0.882815084,0.882735292,1.092199961,0.966173362,1.058814267,1.155276413,1.348008712,1.046850867,1.072178236,0.991644512,1.149285963,0.930156546,1.14740958,1.055477774,1.084085903,0.828702208,0.705419287,1.494829872,1.205445545,1.114766082,1.1404322,1.195817905}; // aaa,aac, aag


stat (infile1,&fst);
xl=fst.st_size;
if (!(inf1=fopen(infile1,"r"))) { fprintf (outerr,"Error opening file: %s\n",infile1); exit(0);}

Seq=(char*)malloc(xl+100)+50; // 100 for ends
if (Seq==NULL) {printf( "Memory problem\n");}

for (i=0;i<=5*5*5-1;i++) {TriplFreq[i]=0; TriplCover[i]=0;}
for (i=0;i<=5*5-1;i++) {DiFreq[i]=0;}
for (i=0;i<=5-1;i++) {MonoFreq[i]=0;}
NumofSeq=0;

while (!feof(inf1))
{ // загружаем очередной сиквенс и вперед
if (!(SeqLength=LoadSeq (inf1,Seq))) {continue;}
SeqT=(char*)malloc(SeqLength+100)+50; // 100 for ends
for (xl=-49;xl<=SeqLength+46;xl++)
    {SeqT[xl]=Seq[xl]+5*Seq[xl+1]+25*Seq[xl+2];}
NumofSeq++;
TotalSeqTLength+=(SeqLength-2);
printf ("SeqLength=%d\n",SeqLength);

//count triplets
for (i=0;i<=5*5*5-1;i++) {TriplTemp[i]=0;}
for (i=0;i<=SeqLength-1-2;i++)
{ TriplFreq[SeqT[i]]++;
  TriplTemp[SeqT[i]]=1;
}
for (i=0;i<=5*5*5-1;i++) TriplCover[i]+=TriplTemp[i]; // coverage by tripets

//count duplets
for (i=0;i<=SeqLength-1-1;i++)
{ DiFreq[Seq[i]+5*Seq[i+1]]++; }

//count mono
for (i=0;i<=SeqLength-1;i++)
{ MonoFreq[Seq[i]]++; }

//count promoter weight
y=0;
for (i=0;i<=SeqLength-1-2;i++) {y+=Class1Weights[4*4*Seq[i]+4*Seq[i+1]+Seq[i+2]]; }
printf ("%-50s %f\n",Name,y/(SeqLength-2));


} //while inf1


//output
printf ("Frequencies\n");

printf ("Mono\n");
for (i1=0;i1<=5-1;i1++) {if (i1==4) continue;
PrintNucleotide(i1); printf ("  ");
printf (" %7d %f\n",MonoFreq[i1], MonoFreq[i1]/(float)TotalSeqTLength);
}

printf ("Di\n");
for (i1=0;i1<=5-1;i1++) {if (i1==4) continue;for (i2=0;i2<=5-1;i2++) {if (i2==4) continue;
PrintNucleotide(i1);PrintNucleotide(i2); printf (" ");
printf (" %7d %f\n",DiFreq[i1+5*i2],DiFreq[i1+5*i2]/(float)TotalSeqTLength); // coverage by tripets
}}

printf ("Tri\n");
for (i1=0;i1<=5-1;i1++) {if (i1==4) continue;for (i2=0;i2<=5-1;i2++) {if (i2==4) continue;for (i3=0;i3<=5-1;i3++) {if (i3==4) continue;
PrintNucleotide(i1);PrintNucleotide(i2);PrintNucleotide(i3);
printf (" %7d %f\n",TriplFreq[i1+5*i2+25*i3], TriplFreq[i1+5*i2+25*i3]/(float)TotalSeqTLength); // coverage by tripets
}}}

/*printf ("Coverage\n");
for (i1=0;i1<=5-1;i1++) {if (i1==4) continue;for (i2=0;i2<=5-1;i2++) {if (i2==4) continue;for (i3=0;i3<=5-1;i3++) {if (i3==4) continue;
PrintNucleotide(i1);PrintNucleotide(i2);PrintNucleotide(i3);
printf (" %f\n",TriplCover[i1+5*i2+25*i3]/(float)NumofSeq); // coverage by tripets
}}}
*/


printf ("TotalSeqTLength=%d\n",TotalSeqTLength);

}

//*******************************************************************************************************************//
int Bend()
{
int i2;
float Curvature, GTilt,GRoll,AccTwist[50];
float TwistDeg[16]={38.9,31.12,32.15,33.81,41.41,34.96,32.91,32.15,41.31,38.5,34.96,31.12,33.28,41.31,41.41,38.9};
float TwistRad[16]={0.678933079,0.543146463,0.561123355,0.59009582,0.722740843,0.610167106,0.574387857,0.561123355,0.720995514,0.671951762,0.610167106,0.543146463,0.580845575,0.720995514,0.722740843,0.678933079};
float Roll[16]={0.5,0.4,2.9,-0.6,1.1,6.5,6.6,2.9,-0.1,-7,6.5,0.4,2.6,-0.1,1.1,0.5};
float Tilt[16]={-0.4,-0.9,-2.6,0,0.6,-1.1,0,-2.6,-0.4,0,-1.1,-0.9,0,-0.4,0.6,-0.4};

stat (infile1,&fst);
xl=fst.st_size;
if (!(inf1=fopen(infile1,"r"))) { fprintf (outerr,"Error opening file: %s\n",infile1); exit(0);}

Seq=(char*)malloc(xl+100)+50; // 100 for ends
if (Seq==NULL) {printf( "Memory problem\n");}
NumofSeq=0;

while (!feof(inf1))
{ // загружаем очередной сиквенс и вперед
if (!(SeqLength=LoadSeq (inf1,Seq))) {continue;}
NumofSeq++;
printf ("SeqLength=%d\n",SeqLength);

//printf ("%d %d - %f\n",Seq[1],Seq[1+1],TwistDeg[4*Seq[1]+Seq[1+1]]); exit(0);

for (i=0;i<=SeqLength-1-1-GCWindow;i++)
{
//count accumulated twist
AccTwist[0]=TwistRad[4*Seq[i]+Seq[i+1]];
for (j=1;j<=GCWindow-1;j++) AccTwist[j]=AccTwist[j-1]+TwistRad[4*Seq[i+j]+Seq[i+j+1]];

//now global Roll and Tilt
GTilt=0;GRoll=0;
for (j=0;j<=GCWindow-1;j++) GTilt+= Tilt[4*Seq[i+j]+Seq[i+j+1]]*cos(AccTwist[j]) + Roll[4*Seq[i+j]+Seq[i+j+1]]*sin(AccTwist[j]);
for (j=0;j<=GCWindow-1;j++) GRoll+= Roll[4*Seq[i+j]+Seq[i+j+1]]*cos(AccTwist[j]) - Tilt[4*Seq[i+j]+Seq[i+j+1]]*sin(AccTwist[j]);
Curvature=sqrt(GTilt*GTilt+GRoll*GRoll);

printf ("%3d %4.3f\n",i,Curvature);

/*if (Curvature > 20) {
printf ("%3d %4.3f ",i,Curvature);
for (j=0;j<=GCWindow-1;j++) PrintNucleotide(Seq[i+j]);
i2=0;
for (j=0;j<=GCWindow-1;j++) if ((Seq[i+j]==0)||(Seq[i+j]==3)) i2++;
printf (" %3.2f",(float)i2/GCWindow);
printf ("\n");
} */




} //next i

} //next seq


}


//*******************************************************************************************************************//
int EntropyProfile()
{
//int TotalSeqTLength=0, TriplFreq[5*5*5],TriplTemp[5*5*5],TriplCover[5*5*5],DiFreq[5*5],MonoFreq[5],i1,i2,i3;
float FreeEnergyChange[16]={-1.2,-1.5,-1.5,-0.9,-1.7,-2.1,-2.8,-1.5,-1.5,-2.3,-2.1,-1.5,-0.9,-1.5,-1.7,-1.2};
float EntropyChange[16]={-21.9,-25.5,-16.4,-15.2,-21,-28.4,-29,-16.4,-23.5,-26.4,-28.4,-25.5,-18.4,-23.5,-21,-21.9};

stat (infile1,&fst);
xl=fst.st_size;
if (!(inf1=fopen(infile1,"r"))) { fprintf (outerr,"Error opening file: %s\n",infile1); exit(0);}


Seq=(char*)malloc(xl+100)+50; // 100 for ends
if (Seq==NULL) {printf( "Memory problem\n");}

//for (i=0;i<=5*5*5-1;i++) {TriplFreq[i]=0; TriplCover[i]=0;}

NumofSeq=0;

while (!feof(inf1))
{ // загружаем очередной сиквенс и вперед
if (!(SeqLength=LoadSeq (inf1,Seq))) {continue;}

NumofSeq++;
TotalSeqLength+=(SeqLength);
printf ("Seq: %s SeqLength=%d\n",Name,SeqLength);


for (i=0;i<=SeqLength-GCWindow;i++)
{
y=0;
for (j=0;j<=GCWindow-1-1;j++) y+=EntropyChange[4*Seq[i+j]+Seq[i+j+1]]; //FreeEnergyChange
//printf ("%-4d %f\n",i,y/GCWindow);
printf ("%f\n",y/GCWindow);
}


} //while inf1

printf ("TotalSeqLength=%d\n",TotalSeqLength);

}



