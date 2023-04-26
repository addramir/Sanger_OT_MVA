cd ~/projects/MVA_output/02_PD_gwas/full_gwas
source activate ldsc

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats nallsEtAl2019_excluding23andMe.txt \
--out nallsEtAl2019_excluding23andMe \
--ignore eaf \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats GIP1_4_traits_ldsc.txt \
--out GIP1_4_traits \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--h2 GIP1_4_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out GIP1_4_traits_h2
#Total Observed scale h2: 0.0167 (0.0018)
#Lambda GC: 1.1207
#Mean Chi^2: 1.1658
#Intercept: 1.0251 (0.0069)
#Ratio: 0.1517 (0.0414)

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--h2 nallsEtAl2019_excluding23andMe.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out nallsEtAl2019_excluding23andMe_h2
#Total Observed scale h2: 0.0162 (0.0015)
#Lambda GC: 1.0926
#Mean Chi^2: 1.1378
#Intercept: 0.9857 (0.0067)
#Ratio < 0 (usually indicates GC correction).

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg nallsEtAl2019_excluding23andMe.sumstats.gz,GIP1_4_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg