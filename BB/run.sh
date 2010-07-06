#@files=( "004DDf.dat",  "005DDM.dat", 006GRM.dat 008SHF.dat 013MBF.dat 014MBF.dat 020WHM.dat 021WHF.dat 022TFF.dat 024KRF.dat 026RLF.dat 027FCF.dat 028SHF.dat 031DDM.dat 036LCF.dat 039NPM.dat 047LMF.dat )

#for i in $files
for i in `ls data/`
do
        ifile=`(echo $i | grep "\.dat")`
	echo R Running $ifile
        #run the R run
	time R < supplement_A-2.R --save --slave -i data/$ifile > data/run.$ifile.txt
	#move our Rplots to data, renaming
	mv Rplots.pdf data/$ifile.Rplots.pdf
done