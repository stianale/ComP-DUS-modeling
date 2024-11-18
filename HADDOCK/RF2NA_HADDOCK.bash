#!/bin/bash

prefix_ComP_removed=$(echo "$1" | sed "s/_ComP/_RF2NA/g")

rm -rf sorted_dirs_list_$1.txt copy_tasklist_$1.sh

cp $1 /path/to/RoseTTAFold2NA/example/
cp $2 /path/to/RoseTTAFold2NA/example/

cd /path/to/RoseTTAFold2NA/example/

lisnum=$(ls -d *$1*$2*/ | wc -l)

if [[ ${lisnum} == 0 ]]; then

source ~/miniconda3/etc/profile.d/conda.sh

conda activate RF2NA

for num in {1..100};
do
../run_RF2NA.sh $1_$2_${num} $1.fa D:$2.fa
done

source ~/miniconda3/etc/profile.d/conda.sh

conda deactivate

else

find ./*$1*$2* -name "model_00.npz" > numpy_list_$1_$2.txt

sed "s/numpy_list.txt/numpy_list_$1_$2.txt/g" read_numpy_general.py > read_numpy_$1_$2.py

best_RF2NA_result=$(python read_numpy_$1_$2.py | sort -V -k2 | head -n 1 | awk '{print $1}' | sed "s/.npz/.pdb/g")

best_RF2NA_result_score=$(python read_numpy_$1_$2.py | sort -V -k2 | head -n 1 | awk '{print $2}')

if (( $(echo "${best_RF2NA_result_score} >= 10" | bc -l) ));

then

echo "Best RF2NA prediction has PAE score: ${best_RF2NA_result_score}"

source ~/miniconda3/etc/profile.d/conda.sh

conda activate RF2NA

tosum_1=$(echo "1")
tosum_2=$(echo "15")
lowernum=$(($lisnum + $tosum_1))
uppernum=$(($lisnum + $tosum_2))

for interval in $(seq "$lowernum" "$uppernum"); 
do
../run_RF2NA.sh "$1_$2_${interval}" "$1.fa" "D:$2.fa"
done

source ~/miniconda3/etc/profile.d/conda.sh

conda deactivate

fi

fi

find ./*$1*$2* -name "model_00.npz" > numpy_list_$1_$2.txt

sed "s/numpy_list.txt/numpy_list_$1_$2.txt/g" read_numpy_general.py > read_numpy_$1_$2.py

best_RF2NA_result=$(python read_numpy_$1_$2.py | sort -V -k2 | head -n 1 | awk '{print $1}' | sed "s/.npz/.pdb/g")

best_RF2NA_result_score=$(python read_numpy_$1_$2.py | sort -V -k2 | head -n 1 | awk '{print $2}')

if (( $(echo "${best_RF2NA_result_score} <= 10" | bc -l) ));

then

best_RF2NA_result_DNA_isolated=$(echo "$1" | sed "s/ComP/$2/g")

echo -e "load ${best_RF2NA_result}\nselect chain B chain C\nremove sele\nsave $1.pdb" > pymol_script_isolate_ComP_$1.pml

echo -e "load ${best_RF2NA_result}\nselect chain A\nremove sele\nsave ${best_RF2NA_result_DNA_isolated}.pdb" > pymol_script_isolate_DUS_$2.pml

pymol -c pymol_script_isolate_ComP_$1.pml

pymol -c pymol_script_isolate_DUS_$2.pml

rm -rf pymol_script_isolate_ComP_$1.pml pymol_script_isolate_DUS_$2.pml

rm -rf sorted_dirs_list_$1.txt

mkdir Pipeline_HADDOCK_${prefix_ComP_removed}
cd Pipeline_HADDOCK_${prefix_ComP_removed}
cp ../Find_normal_modes.R .
cp /path/to/RoseTTAFold2NA/example/$1.pdb .
cp /path/to/RoseTTAFold2NA/example/${best_RF2NA_result_DNA_isolated}.pdb .

sed -i "s/Pipeline_HADDOCK/Pipeline_HADDOCK_${prefix_ComP_removed}/g" Find_normal_modes.R

Rscript Find_normal_modes.R $1.pdb

Rscript Find_normal_modes.R ${best_RF2NA_result_DNA_isolated}.pdb

for mode in $(ls Normal_mode*.pdb | sed "s/.pdb//g" | sed "/model/d");
do
/home/stian/miniconda3/bin/pdb_selmodel -10 "${mode}".pdb > "${mode}"_model_10.pdb
/home/stian/miniconda3/bin/pdb_selmodel -27 "${mode}".pdb > "${mode}"_model_27.pdb
reduce -BUILD "${mode}"_model_10.pdb > protonated_"${mode}"_model_10.pdb
reduce -BUILD "${mode}"_model_27.pdb > protonated_"${mode}"_model_27.pdb
done

pdb_mkensemble protonated_Normal_mode_*ComP_7_model_10.pdb protonated_Normal_mode_*ComP_8_model_10.pdb protonated_Normal_mode_*ComP_9_model_10.pdb protonated_Normal_mode_*ComP_10_model_10.pdb protonated_Normal_mode_*ComP_11_model_10.pdb protonated_Normal_mode_*ComP_12_model_10.pdb protonated_Normal_mode_*ComP_13_model_10.pdb protonated_Normal_mode_*ComP_14_model_10.pdb protonated_Normal_mode_*ComP_15_model_10.pdb protonated_Normal_mode_*ComP_16_model_10.pdb | sed "s/ ATOM/\nATOM/g" | sed "s/ ENDMDL/\nENDMDL/g" > 10_lowest_$1_model_10.pdb

pdb_mkensemble protonated_Normal_mode_*ComP_7_model_27.pdb protonated_Normal_mode_*ComP_8_model_27.pdb protonated_Normal_mode_*ComP_9_model_27.pdb protonated_Normal_mode_*ComP_10_model_27.pdb protonated_Normal_mode_*ComP_11_model_27.pdb protonated_Normal_mode_*ComP_12_model_27.pdb protonated_Normal_mode_*ComP_13_model_27.pdb protonated_Normal_mode_*ComP_14_model_27.pdb protonated_Normal_mode_*ComP_15_model_27.pdb protonated_Normal_mode_*ComP_16_model_27.pdb | sed "s/ ATOM/\nATOM/g" | sed "s/ ENDMDL/\nENDMDL/g" > 10_lowest_$1_model_27.pdb

for somefile in $(ls protonated*pdb | sed "/ComP/d"); do pdbfixer "${somefile}" --replace-nonstandard --output pdbfixed_"${somefile}"; pdb_selchain -B,C pdbfixed_"${somefile}" | pdb_chain -B | pdb_reres -1 | pdb_tidy > chain_B_"${somefile}"; done


pdb_mkensemble chain_B_protonated_Normal_mode_*DUS_7_model_10.pdb chain_B_protonated_Normal_mode_*DUS_8_model_10.pdb chain_B_protonated_Normal_mode_*DUS_9_model_10.pdb chain_B_protonated_Normal_mode_*DUS_10_model_10.pdb chain_B_protonated_Normal_mode_*DUS_11_model_10.pdb | sed "s/ ATOM/\nATOM/g" | sed "s/ ENDMDL/\nENDMDL/g" > 10_lowest_"${prefix_ComP_removed}"_DUS_model_10.pdb

pdb_mkensemble chain_B_protonated_Normal_mode_*DUS_7_model_27.pdb chain_B_protonated_Normal_mode_*DUS_8_model_27.pdb chain_B_protonated_Normal_mode_*DUS_9_model_27.pdb chain_B_protonated_Normal_mode_*DUS_10_model_27.pdb chain_B_protonated_Normal_mode_*DUS_11_model_27.pdb | sed "s/ ATOM/\nATOM/g" | sed "s/ ENDMDL/\nENDMDL/g" > 10_lowest_"${prefix_ComP_removed}"_DUS_model_27.pdb

pdb_mkensemble chain_B_protonated_Normal_mode_*DUS_12_model_10.pdb chain_B_protonated_Normal_mode_*DUS_13_model_10.pdb chain_B_protonated_Normal_mode_*DUS_14_model_10.pdb chain_B_protonated_Normal_mode_*DUS_15_model_10.pdb chain_B_protonated_Normal_mode_*DUS_16_model_10.pdb | sed "s/ ATOM/\nATOM/g" | sed "s/ ENDMDL/\nENDMDL/g" > 10_highest_"${prefix_ComP_removed}"_DUS_model_10.pdb

pdb_mkensemble chain_B_protonated_Normal_mode_*DUS_12_model_27.pdb chain_B_protonated_Normal_mode_*DUS_13_model_27.pdb chain_B_protonated_Normal_mode_*DUS_14_model_27.pdb chain_B_protonated_Normal_mode_*DUS_15_model_27.pdb chain_B_protonated_Normal_mode_*DUS_16_model_27.pdb | sed "s/ ATOM/\nATOM/g" | sed "s/ ENDMDL/\nENDMDL/g" > 10_highest_"${prefix_ComP_removed}"_DUS_model_27.pdb


for file in $(ls 10_lowest*ComP*); 
do
for file2 in $(ls 10_lowest*DUS* 10_highest*DUS*);      
do
newdir=$(echo "${file}_${file2}" | sed "s/.pdb//g" | sed "s/ComP_${prefix_ComP_removed}.*DUS/ComP/g");
echo -e "${file}\t${file2}\t${newdir}" >> sorted_dirs_list_$1.txt;
done
done

sort -r -k3 -V sorted_dirs_list_$1.txt > tmp.txt
mv tmp.txt sorted_dirs_list_$1.txt

source ~/miniconda3/etc/profile.d/conda.sh

conda activate haddock3

while read -r pipeline_ComP pipeline_DUS sorted;
	do
	if [ ! -d "${sorted}" ]; then
    mkdir "${sorted}"
	mkdir "${sorted}"/data
	cp ${prefix_ComP_removed}/"${pipeline_ComP}" "${sorted}"/data
	cp ${prefix_ComP_removed}/"${pipeline_DUS}" "${sorted}"/data
	cp $3 "${sorted}"/data
	cp config_to_copy.cfg "${sorted}"
	cd "${sorted}"/data
	
	fileprefix=$(echo ${pipeline_ComP} | sed "s/.pdb//g")
	for newnum in {1..10}; do pdb_selmodel -"${newnum}" "${fileprefix}".pdb > model_"${newnum}"_"${fileprefix}".pdb; done
	freesasa model_"${newnum}"_"${fileprefix}".pdb --format=rsa > model_"${newnum}"_"${fileprefix}".rsa
	awk '{if (NF==13 && $5>40) print $0; if (NF==14 && $6>40) print $0}' model_"${newnum}"_"${fileprefix}".rsa | awk '{print $3,$2,$4}' | sed "s/ /\./g" | sed "s/^/R /g" > surface_accessible_model_"${newnum}"_"${fileprefix}".txt
	cat surface_accessible_model_"${newnum}"_"${fileprefix}".txt >> master_file_surface_accessible.rsa
	grep -f $3 master_file_surface_accessible.rsa > tmp_ambig_"${fileprefix}".txt
	awk -F'.' '{print $3}' tmp_ambig_"${fileprefix}".txt | tr '\n' ' ' | sed "s/$/\n /g" > active_ComP_"${fileprefix}".txt
	echo -e "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24\n " > active_DNA.txt
	python ../../active-passive-to-ambig.py active_ComP_"${fileprefix}".txt active_DNA.txt > final_ambig_"${fileprefix}".tbl
	rm -rf model_"${newnum}"_"${fileprefix}".pdb
	
	cd ..
	
	sed -i "s/freesasa_filtered_ambig.tbl/final_ambig_${fileprefix}.tbl/g" config_to_copy.cfg
	
	sed -i "s/sampling = 4000/sampling = 50000/g" config_to_copy.cfg
	
	sed -i "s/delenph = false/delenph = true/g" config_to_copy.cfg
	
	sed -i "s/randremoval = false/randremoval = true/g" config_to_copy.cfg
	
	sed -i "s/output_protein_final_2nba.pdb/${pipeline_ComP}/g" config_to_copy.cfg
	
	sed -i "s/output_dna_final_2nba.pdb/${pipeline_DUS}/g" config_to_copy.cfg
	
	haddock3 config_to_copy.cfg
	
	cd ..
	
	fi
	
	done < sorted_dirs_list_$1.txt