
#!/bin/bash

N_FILES=9
OFILE=application_week.R
for((i=1;i<=$N_FILES;i++)); do
	FNAME="week_${i}.R"
	FNAME2="week_${i}"
	echo $FNAME
	cp $OFILE $FNAME
	sed -i -e "s/i = 0/i = ${i}/g" $FNAME
	cp RunFile.bash RunFile_temp.bash
	sed -i -e "s|filename|${FNAME}|g" RunFile_temp.bash
	sed -i -e "s|out|${FNAME2}|g" RunFile_temp.bash
	bsub < RunFile_temp.bash
	rm RunFile_temp.bash
done
