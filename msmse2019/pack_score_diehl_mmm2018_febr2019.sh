#!/bin/bash

##organize SCORE simultion result run18
##Markus Kuehbach, 2019/02/23, m.kuehbach at mpie.de


mkdir -p Rediscretized
mkdir -p Profiling
mkdir -p SingleGrainData
mkdir -p ThreadProfilingGrowth
mv *.Rediscretized.* Rediscretized/
mv *.ThreadProfilingGrowth.csv ThreadProfilingGrowth/
mv *.SingleGrainData.csv SingleGrainData/
mv *.ProfilingLog.csv Profiling/
mv *.STD* Profiling/

#rm -rf ang/
#rm -rf png/
#rm -rf prof/
#rm -rf res/

#simid=6403
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=6665
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=6666
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=64030
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=64031
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=64032
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=64033
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

simid=64034
mv *.$simid*scalebar.* SCORE_$simid/
cd SCORE_$simid/
tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
cd ..

simid=64035
mv *.$simid*scalebar.* SCORE_$simid/
cd SCORE_$simid/
tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
cd ..

#simid=66660
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=66661
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=66662
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

#simid=66663
#mv *.$simid*scalebar.* SCORE_$simid/
#cd SCORE_$simid/
#tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
#tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
#cd ..

simid=66664
mv *.$simid*scalebar.* SCORE_$simid/
cd SCORE_$simid/
tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
cd ..

simid=66665
mv *.$simid*scalebar.* SCORE_$simid/
cd SCORE_$simid/
tar -czvf SCORE_$simid.Snapshots2D.ang.tar.gz ang/
tar -czvf SCORE_$simid.Images2D.png.tar.gz png/
cd ..

