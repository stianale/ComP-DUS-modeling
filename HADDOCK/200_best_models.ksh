#!/bin/ksh93

rm -rf $1*pdb.gz

rm -rf $1_load_and_save.py

rm -rf $1_number_residues_ensemble.py

rm -rf Good_models_$1.txt

rm -rf copy_tasklist_$1.sh

rm -rf score_file_$1.txt

rm -rf best_models_$1.txt

caprimodels=$(ls ./{10_lowest,10_highest,10}_$1*[1-9]*_model*/run1-test/analysis/11_caprieval_analysis/capri_ss.tsv)

for capri in $caprimodels;
do
fivebest=$(cat ${capri})
modified_capri=$(echo "${capri}" | sed "s/\//_/g")
echo "${fivebest}" | sed "s/10_seletopclusts/${modified_capri}/g" | sed "s/.*\._10_/10_/g" | sed "s/_run1-test_analysis_11_caprieval_analysis_capri_ss.tsv/\/run1-test\/10_seletopclusts/g" | sed "s/.pdb/.pdb.gz/g" | awk '{print $1,$4,$5,$6,$7,$8,$9}' | sed "/model score/d" >> score_file_$1.txt
done

sort -V -k2 -r score_file_$1.txt > tmp.txt
mv tmp.txt score_file_$1.txt
#cat score_file_$1.txt

while IFS=$' ' read -r file score irmsd fnat lrmsd ilrmsd dockq; do
    newfile=$(echo "${file}" | sed "s/\/run1-test\/10_seletopclusts\//_/g")
    #echo -e "${file}\t${score}${irmsd}\t${lrmsd}\t${ilrsmd}\t${dockq}"
	#echo "${newfile}"
	echo "${file} ${newfile} ${score} ${irmsd} ${fnat} ${lrmsd} ${ilrmsd} ${dockq}" >> best_models_$1.txt
done < score_file_$1.txt

head -n 200 best_models_$1.txt > tmp.txt
mv tmp.txt best_models_$1.txt

while IFS=$' ' read -r file newfile score irmsd fnat lrmsd ilrmsd dockq;
do
#if (( $(echo "${fnat} >= 0.5" | bc -l) )) && (( $(echo "${irmsd} <= 1" | bc -l) )) || (( $(echo "${fnat} >= 0.5" | bc -l) )) && (( $(echo "${lrmsd} <= 1" | bc -l) )); 
#then
echo "cp ${file} ${newfile}" >> copy_tasklist_$1.sh
#fi
done < best_models_$1.txt

#head -n 40 copy_tasklist_$1.sh > tmp.txt
#mv tmp.txt copy_tasklist_$1.sh
sh copy_tasklist_$1.sh

awk '{print $3}' copy_tasklist_$1.sh > Good_models_$1.txt

cp load_pdb_files_and_save_ensemble.py $1_load_and_save.py

sed -i "s/N_subflava/$1/g" $1_load_and_save.py

sed -i 's/prefix +/"10_*" + prefix +/g' $1_load_and_save.py

sed -i "s/multimodel.pdb/200_best_multimodel_$1.pdb/g" $1_load_and_save.py

pymol -c $1_load_and_save.py

cp number_residues_ensemble.py $1_number_residues_ensemble.py

sed -i "s/multimodel.pdb/200_best_multimodel_$1.pdb/g" $1_number_residues_ensemble.py

rm -rf Good_models_$1.txt

rm -rf copy_tasklist_$1.sh

#rm -rf score_file_$1.txt

#pymol $1_number_residues_ensemble.py

rm -rf $1_load_and_save.py

rm -rf $1_number_residues_ensemble.py
