mm.exe matrixTFP64.lib Breast500.seq  Breast500.out all64.prf 
mm.exe matrixTFP64.lib BreastNot500.seq  BreastNot500.out  all64.prf 

mmmodules.exe  ModulesBestSitesAll.dat  Breast500.out.bin   BreastNot500.out.bin   MinCoverY 0.75 MaxCoverN 0.5 FreqDiff 1.5 >  Breast500.CRMout    
