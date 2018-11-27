
#!/bin/bash

N_FILES=100
FUN="npqw"
OFILE="Sim_${FUN}.R"

for((i=1;i<=$N_FILES;i++))
do 
	FNAME="Sim_${FUN}_${i}.R"
	FNAME2="Sim_${FUN}_${i}"
	echo $FNAME

	cp $OFILE $FNAME
	sed -i -e "s/i = 0/i = ${i}/g" $FNAME

	cp RunFile.bash RunFile_temp.bash
	sed -i -e "s|filename|${FNAME}|g" RunFile_temp.bash
	sed -i -e "s|out|${FNAME2}|g" RunFile_temp.bash
	bsub < RunFile_temp.bash
	rm RunFile_temp.bash
done
