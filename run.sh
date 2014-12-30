#~/bin/sh
mkdir -p out
t=`date "+%d%m%Y"`
./generate_taxonomy.py --output out/taxonomy_${t}.txt --output_red out/taxonomy_reduced_${t}.txt --pickle out/taxonomy_${t}.pkl | tee out/generate_taxonomy_${t}.txt
./repophlan_get_microbes.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir out/microbes_${t} --nproc 15 --out_summary out/repophlan_microbes_${t}.txt | tee out/repophlan_microbes_${t}.log
python screen.py --nproc 10 --in_summary out/repophlan_microbes_${t}.txt --out_summary out/repophlan_microbes_${t}_wscores.txt
./repophlan_get_viruses.py --taxonomy out/taxonomy_reduced_${t}.txt --out_dir out/viruses_${t} --out_summary out/repophlan_viruses_${t}.txt | tee out/repophlan_viruses_${t}.log
