# под флагом а нужно поставить нужный файл с ридами. полученные данные перенести в папку 'dna_reads'
bedtools intersect -a K562_chr18_19.tab -b exon_bed.dms -wa > 'dna_read_exon'
bedtools intersect -a K562_chr18_19.tab -b intron_bed.dms -wa > 'dna_read_intron'
bedtools intersect -a K562_chr18_19.tab -b upstream_200_bed.dms -wa > 'dna_read_upstream'
bedtools intersect -a K562_chr18_19.tab -b downstream_200_bed.dms -wa > 'dna_read_downstream'
bedtools intersect -a K562_chr18_19.tab -b exon_bed.dms intron_bed.dms -v -wa > 'dna_read_nongenic'

