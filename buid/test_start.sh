N_MACHINE="3 5 10"
for f in $(find .dat ../resources)
do
	for i in $N_MACHINE
	do
		echo $f $i
	done
done
now=$(date +"%Y_%m_%d_%H_%M")
R_DIR="../results/${N_MACHINE// /_}"

echo $now
mkdir -p $R_DIR

echo "test" >> $R_DIR/result_$now
