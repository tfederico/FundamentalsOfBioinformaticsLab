#!/bin/bash

python run_local_blast.py -ids ../data/uniprot_id_list.txt -q ../queries/ -db ../db/db.fasta -o ../results/blast_out.txt -evalue 0.5 -opng ../results/DistrubutionEvalue05.png

python run_local_blast.py -ids ../data/uniprot_id_list.txt -q ../queries/ -db ../db/db.fasta -o ../results/psiblast_out.txt -evalue 0.5 -opng ../results/DistrubutionEvalue05_psi.png -psi

python run_local_blast.py -ids ../data/uniprot_id_list.txt -q ../queries/ -db ../db/db.fasta -o ../results/blast_out.txt -evalue 1 -opng ../results/DistrubutionEvalue1.png

python run_local_blast.py -ids ../data/uniprot_id_list.txt -q ../queries/ -db ../db/db.fasta -o ../results/psiblast_out.txt -evalue 1 -opng ../results/DistrubutionEvalue1_psi.png -psi

cd ../results/
git add DistubutionEvalue*
git commit -m "Fix plot scales"
git push
