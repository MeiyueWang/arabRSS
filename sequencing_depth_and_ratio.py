id=$1
strand=$2
file=${id}.uniq.cs.bam
bedfile=arab_${strand}.bed
test -e ${id}.${strand}.extract.bam || bedtools intersect -a ${bedfile} -b ${file} -wa -wb >tmp_${id}.${strand}.intersect.bed
test -e ${id}.${strand}.extract.bam || cat tmp_${id}.${strand}.intersect.bed |awk 'BEGIN{OFS="\t"}{if($4!=$10)print $8}' >tmp_${id}.${strand}.extract.list
test -e tmp_${id}.${strand}.extract.DMS.xls || samtools view ${file} |fgrep -w -f tmp_${id}.${strand}.extract.list |samtools view -Sh -b -t ~/genome/arabidopsis/TAIR10_chr_all.fas.fai -o ${id}.${strand}.extract.bam -
test -e tmp_${id}.${strand}.extract.DMS.xls || python DMS_ratio_20180709.py ${id}.${strand}.extract.bam tmp_${id}.${strand}.extract.DMS.xls
cat tmp_${id}.${strand}.extract.DMS.xls |sort -k1,1 |awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' |sort -k1,1 -k2,2n >tmp_${id}.${strand}.extract.DMS.sorted.bed
bedtools bamtobed -i ${id}.${strand}.extract.bam >tmp_${id}.${strand}.extract.bed
sort -k1,1 -k2,2n tmp_${id}.${strand}.extract.bed >${id}.${strand}.extract.sorted.bed
bedtools coverage -sorted -a tmp_${id}.${strand}.extract.DMS.sorted.bed -b ${id}.${strand}.extract.sorted.bed |cut -f1-6 |awk 'BEGIN{OFS="\t"}{print $0,$5/$6}' >${id}.${strand}.extract.DMSratio.bed
cat ${id}.${strand}.extract.DMSratio.bed |awk '$6>=5' |sort -k1,1 -k2,2n|tee ${id}.${strand}.extract.DMSratio.filterdp5.bed|cut -f4,5|awk 'BEGIN{OFS="\t"}{a[$1]+=$2}END{for(i in a)print i,a[i]}' >tmp_stat_${id}.${strand}.nuleocompos.txt
if [ ${strand} == "minus" ]
then
cat ${id}.${strand}.extract.DMSratio.filterdp5.bed |awk '$4=="G"||$4=="T"' >${id}.${strand}.extract.DMSratio.filterdp5.AC.bed
cat tmp_stat_${id}.${strand}.nuleocompos.txt |awk 'BEGIN{OFS="\t"}{if($1=="A")print "T",$2;else if($1=="T")print "A",$2;else if($1=="C")print "G",$2;else if($1=="G")print "C",$2}' |grep -v "N"|sort -k1,1 >stat_${id}.${strand}.nuleocompos.txt
else
cat ${id}.${strand}.extract.DMSratio.filterdp5.bed |awk '$4=="A"||$4=="C"' >${id}.${strand}.extract.DMSratio.filterdp5.AC.bed
cat tmp_stat_${id}.${strand}.nuleocompos.txt |grep -v "N"|sort -k1,1 >stat_${id}.${strand}.nuleocompos.txt
fi
#paste stat_${id}.plus.nuleocompos.txt stat_${id}.minus.nuleocompos.txt |awk 'BEGIN{OFS="\t"}{print $2+$4}'|sed '1i '${id}'' >tmp_${id}.nucleocompos.sta
