#~/bin/sh
mkdir -p out
t=`date "+%d%m%Y"`
./generate_taxonomy.py --output out/taxonomy_${t}.txt --output_red out/taxonomy_reduced_${t}.txt --pickle out/taxonomy_${t}.pkl | tee out/generate_taxonomy_${t}.txt
./repophlan_get_viruses.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir out/viruses_${t} --out_summary out/repophlan_viruses_${t}.txt | tee out/repophlan_viruses_${t}.log
./repophlan_get_microbes.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir out/microbes_${t} --nproc 15 --out_summary out/repophlan_microbes_${t}.txt | tee out/repophlan_microbes_${t}.log
cut -f 2 out/repophlan_microbes_${t}.txt  | cut -f 1-7 -d '|' | sort | uniq -c | sort -k 1 -g -r > out/repophlan_microbes_${t}.u.txt
./repophlan_get_euks.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir_fungi out/fungi_${t} --out_dir_protozoa out/protozoa_${t} --out_summary_fungi out/repophlan_fungi_${t}.txt --out_summary_protozoa out/repophlan_protozoa_${t}.txt | tee out/repophlan_euks_${t}.log
