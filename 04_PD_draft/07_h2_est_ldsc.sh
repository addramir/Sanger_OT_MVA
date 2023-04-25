cd ~/projects/MVA_output/02_PD_gwas/full_gwas
source activate ldsc

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats nallsEtAl2019_excluding23andMe.txt \
--out allsEtAl2019_excluding23andMe \
--ignore eaf \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist

~/projects/MVA_output/02_PD_gwas/ldsc/munge_sumstats.py \
--sumstats GIP1_4_traits.txt \
--out GIP1_4_traits \
--merge-alleles ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/w_hm3.noMHC.snplist


~/projects/MVA_output/02_PD_gwas/ldsc/ldsc.py \
--h2 GIP1_4_traits.txt \
--ref-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--w-ld-chr ~/projects/MVA_output/02_PD_gwas/ldsc/eur_ldscores/ \
--out GIP1_4_traits_h2