

rm *Stroma*

mkdir -p L9/P1 L9/P12L1 L9/P12L2 L9/P13 251B/ 201B/ 193B/P4L1 193B/P4L3 193B/P5 193B/P10L1 193B/P10L2 193B/P11L1 193B/P11L2

for file in *P11L1C1*; do cp $file ${file/P11L1C1/P11L2C1}; done

for file in *P10L3C3*; do mv $file ${file/P10L3C3/P10L2C3}; done

for file in *P10L1C4*; do cp $file ${file/P10L1C4/P10L2C5}; done

mv *L9_P1C* L9/P1/
mv *L9_P12L1* L9/P12L1/
mv *L9_P12L2* L9/P12L2/
mv *P13* L9/P13/
mv *251* 251B/
mv *201* 201B/
mv *P5* 193B/P5/
mv *P10L1* 193B/P10L1/
mv *P10L2* 193B/P10L2/
mv *P11L1* 193B/P11L1/
mv *P11L2* 193B/P11L2/
mv *P4L1* 193B/P4L1/
mv *P4L3* 193B/P4L3/



mv *L9_P1Bulk* L9/P1/ 
cp *P12Bulk* ./L9/P12L1/
mv *P12Bulk* ./L9/P12L2/
cp *P4Bulk* 193B/P4L1/
mv *P4Bulk* 193B/P4L3/
cp *P10Bulk* 193B/P10L1/
mv *P10Bulk* 193B/P10L2/
cp *P11Bulk* 193B/P11L1/
mv *P11Bulk* 193B/P11L2/