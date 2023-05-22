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

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats GIP1_5_traits.txt \
--out GIP1_5_traits \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats GIP1_5_traits_neff.txt \
--out GIP1_5_traits_neff \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats 20230518_MA_over_14.txt \
--out 20230518_MA_over_14 \
--signed-sumstat b_ma,0 \
--N 350000 \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats 20230518_MA_over_all.txt \
--out 20230518_MA_over_all \
--signed-sumstat b_ma,0 \
--N 350000 \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats 20230518_MA_over_1-4.txt \
--out 20230518_MA_over_1-4 \
--signed-sumstat b_ma,0 \
--N 350000 \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist





~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats ~/projects/MVA_output/02_PD_gwas/generic-metal/METAANALYSIS1.TBL \
--out MA_F_N \
--chunksize 500000 \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist \
--snp MarkerName \
--frq Freq1 \
--a1 Allele1 \
--a2 Allele2 \
--p P-value \
--signed-sumstat Effect,0 \
--N-col TotalSampleSize



~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--h2 GIP1_5_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out GIP1_5_traits_h2
#Total Observed scale h2: 0.0173 (0.0019)
#Lambda GC: 1.1238
#Mean Chi^2: 1.1562
#Intercept: 1.03 (0.0067)
#Ratio: 0.1918 (0.0427)




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
--h2 GIP1_5_traits_neff.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out GIP1_5_traits_neff_h2


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--h2 MA_F_N.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out MA_F_N
#Total Observed scale h2: 0.0104 (0.001)
#Lambda GC: 1.1144
#Mean Chi^2: 1.1538
#Intercept: 1.0029 (0.0073)
#Ratio: 0.019 (0.0478)


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg nallsEtAl2019_excluding23andMe.sumstats.gz,GIP1_5_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg
#0.9319  0.0151  61.7506  0.

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg nallsEtAl2019_excluding23andMe.sumstats.gz,GIP1_4_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg
#0.9788  0.0073  133.8661  0.0

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg GIP1_4_traits.sumstats.gz,GIP1_5_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg
#0.8804  0.0175  50.1893  0.0


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg MA_F_N.sumstats.gz,GIP1_5_traits.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg MA_F_N.sumstats.gz,GIP1_5_traits_neff.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg GIP1_5_traits.sumstats.gz,GIP1_5_traits_neff.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg MA_F_N.sumstats.gz,20230518_MA_over_1-4.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg MA_F_N.sumstats.gz,20230518_MA_over_all.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg MA_F_N.sumstats.gz,20230518_MA_over_14.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg 20230518_MA_over_all.sumstats.gz,20230518_MA_over_14.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg

~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--rg 20230518_MA_over_all.sumstats.gz,20230518_MA_over_1-4.sumstats.gz \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out rg


