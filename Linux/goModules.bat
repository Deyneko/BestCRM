echo "Set 1" %TIME%
mmmodules.exe AllModules.dat bsubtilis_promoters_250bp.out.bin bsubtilis_background_250bp.out.bin /window 200 FreqDiff 1.0
echo "Set 1" %TIME%

echo "Set 2" %TIME%
mmmodules.exe AllModules.dat ecoli_promoters_250bp.out.bin ecoli_background_250bp.out.bin /window 200 FreqDiff 1.0
echo "Set 2" %TIME%

echo "Set 3" %TIME%
mmmodules.exe AllModules.dat pao1_promoters_250bp.out.bin pao1_background_250bp.out.bin /window 200 FreqDiff 1.0
echo "Set 3" %TIME%

