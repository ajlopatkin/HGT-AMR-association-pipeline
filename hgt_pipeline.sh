#!/bin/bash


CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

#Arguments

#Location of fasta files
in_dir=$1

#results directory
out_dir=$2

#mefinder coverage: Float, 0 to 1
cov_mef=$3

#Abricate coverage: Int, 0 to 100
cov_abr=$4

#Mob coverage: Int, 0 to 100
cov_mob=$5

#Abricate percent identity: Int, 0 to 100
pid_abr=$6

#Mob percent identity: Int, 0 to 100
pid_mob=$7

#Mefinder e-value: Int, 1eX, example -6 is 1e-6
e_mef=$8

#Mob evalue: Float, min allowed 0.00001
e_mob=$9

#Visuals: 0 = no, 1 = yes
visuals=$10



FILES="$in_dir"/*.fasta

mkdir "$out_dir"/raw_output


#Get ARGS with abricate on CARD

mkdir "$out_dir"/raw_output/card

echo Running Abricate db:card...

conda activate abricate_env

abricate --minid "$pid_abr" --mincov "$cov_abr" --db card "$in_dir"/*.fasta > "$out_dir"/raw_output/card/card_results.txt 

abricate --summary "$out_dir"/raw_output/card/card_results.txt > "$out_dir"/raw_output/card/card_summary.txt

conda deactivate 



#Abricate - plasmidfinder

mkdir "$out_dir"/raw_output/plasmidfinder

echo Running Abricate db: plasmidfinder...

conda activate abricate_env

abricate --minid "$pid_abr" --mincov "$cov_abr" --db plasmidfinder "$in_dir"/*.fasta > "$out_dir"/raw_output/plasmidfinder/plasmidfinder_results.txt

abricate --summary "$out_dir"/raw_output/plasmidfinder/plasmidfinder_results.txt > "$out_dir"/raw_output/plasmidfinder/plasmidfinder_summary.txt

conda deactivate



#Mobile Element Finder

mkdir "$out_dir"/raw_output/mefinder

conda activate mefinder_env

for f in $FILES
do
  f_name=$(basename $f)
  f_name=${f_name%.fasta}
  echo Running mefinder on $f_name
  mefinder find --contig "$f" --min-coverage "$cov_mef" --max-evalue "$e_mef" "$out_dir"/raw_output/mefinder/"$f_name"
done

conda deactivate



#Mobtyper
mkdir "$out_dir"/raw_output/mobtyper

conda activate updated_mob_env

for f in $FILES
do
  f_name=$(basename $f)
  f_name=${f_name%.fasta}
  echo Running mobtyper on $f_name
  mob_recon --infile "$f" --outdir "$out_dir"/raw_output/mobtyper/"$f_name" --min_con_ident "$pid_mob" --min_rep_ident "$pid_mob" --min_con_cov "$cov_mob" --min_rep_cov "$cov_mob" --min_con_evalue "$e_mob" --min_rep_evalue "$e_mob"
done



#Generate output
mkdir "$out_dir"/compiled_output

mkdir "$out_dir"/plasmid_diagrams

conda activate genomediagram_env

echo compiling results...

python hgt_analyze.py -d "$out_dir" -f "$in_dir" -v "$visuals"

conda deactivate 

