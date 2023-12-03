mm.exe matrixTFP81.lib Breast500.seq  Breast500.out all81.prf 
mm.exe matrixTFP81.lib BreastNot500.seq  BreastNot500.out  all81.prf 

mmmodules.exe  ModulesBestSitesAll.dat  Breast500.out.bin   BreastNot500.out.bin   MinCoverY 0.75 MaxCoverN 0.5 FreqDiff 1.5 >  Breast500.CRMout    
