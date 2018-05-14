#Вставить под флаг -а название файла с рнк-днк контактами (можно получить в папке donut)
bedtools intersect -a K562_chr18_19.tab -b embl_edited_2 -wa -wb >rna_intersected
