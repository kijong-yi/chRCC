# Running IGV with batch file (screen shot saver)
# https://software.broadinstitute.org/software/igv/batch
new
genome /home/users/kjyi/igv/genomes/chRSRS.genome
expand
snapshotDirectory /home/users/kjyi/Projects/chromophobe/igv/RSRS/coding_mutations
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC01/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:8231-8914
snapshot chRCC01_100_MT-ATP6__8573A_G16D.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC02/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:14729-15412
snapshot chRCC02_098_MT-CYB___15071C_Y109H.png
goto chrRSRS:14160-14843
snapshot chRCC02_098_MT-ND6___14502C_I58V.png
goto chrRSRS:14876-15559
snapshot chRCC02_098_MT-CYB___15218G_T158A.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC03/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:12614-13297
snapshot chRCC03_099_MT-ND5___12956G_N207S.png
goto chrRSRS:3877-4560
snapshot chRCC03_100_MT-ND1___4219A_V305I.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC04/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:4839-5522
snapshot chRCC04_100_MT-ND2___5181G_T238A.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC05/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:8154-8837
snapshot chRCC05_100_MT-ATP8__8496C_M44T.png
goto chrRSRS:7014-7697
snapshot chRCC05_100_MT-CO1___7356A_V485M.png
goto chrRSRS:10809-11492
snapshot chRCC05_100_MT-ND4___11151T_A131V.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC07/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:15215-15898
snapshot chRCC07_079_MT-CYB___15557A_E271K.png
goto chrRSRS:7488-8171
snapshot chRCC07_100_MT-CO2___7830A_R82H.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC10/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:12721-13404
snapshot chRCC10_097_MT-ND5___13063A_V243I.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC12/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:6979-7662
snapshot chRCC12_025_MT-CO1___7321A_Stop-gain.png
goto chrRSRS:13907-14590
snapshot chRCC12_100_MT-ND6___14249A_A142V.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC17/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:8668-9351
snapshot chRCC17_023_MT-ATP6__9010A_A162T.png
new
load /home/users/kjyi/Projects/chromophobe/mtoolbox/OUT_chRCC18/OUT2-sorted.bam
collapse OUT2-sorted.bam
goto chrRSRS:8603-9286
snapshot chRCC18_100_MT-ATP6__8945C_M140T.png
