echo "Set 1" %TIME%
mmmodules.exe BacillusModules.dat bsubtilis_promoters_250bp.out.bin bsubtilis_background_250bp.out.bin /window 100 FreqDiff 1.0 MinCoverY 0.5 MaxCoverN 0.25 /FullOutput 3
echo "Set 1" %TIME%

echo "Set 2" %TIME%
mmmodules.exe EcoliModules.dat ecoli_promoters_250bp.out.bin ecoli_background_250bp.out.bin         /window 100 FreqDiff 1.0 MinCoverY 0.5 MaxCoverN 0.25 /FullOutput 3
echo "Set 2" %TIME%

echo "Set 3" %TIME%
mmmodules.exe PseudomonasModules.dat pao1_promoters_250bp.out.bin pao1_background_250bp.out.bin       /window 100 FreqDiff 1.0 MinCoverY 0.5 MaxCoverN 0.25 /FullOutput 3
echo "Set 3" %TIME%

echo "Set 1" %TIME%
mmmodules.exe BacillusModules.dat bsubtilis_promoters_250bp.out.bin bsubtilis_background_250bp.out.bin /window 100 FreqDiff 1.0 MinCoverY 0.25 MaxCoverN 0.1 /FullOutput 3
echo "Set 1" %TIME%

echo "Set 2" %TIME%
mmmodules.exe EcoliModules.dat  ecoli_promoters_250bp.out.bin ecoli_background_250bp.out.bin         /window 100 FreqDiff 1.0 MinCoverY 0.25 MaxCoverN 0.1 /FullOutput 3
echo "Set 2" %TIME%

echo "Set 3" %TIME%
mmmodules.exe PseudomonasModules.dat pao1_promoters_250bp.out.bin pao1_background_250bp.out.bin           /window 100 FreqDiff 1.0 MinCoverY 0.25 MaxCoverN 0.1 /FullOutput 3
echo "Set 3" %TIME%

