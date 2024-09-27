dir= #edit
species= #edit
cell= #edit
core= #edit
ref_dir= #edit
basalkit= #edit
fottrack_tool=  #edit

if [ "$" == "mouse" ]; then
    ref_base=$ref_dir/Mus_musculus/mm10
elif [ "$species" == "human" ]; then
    ref_base=$ref_dir/Homo_sapiens/GRCh38
else
    echo "Invalid species. Please specify 'mouse' or 'human'."
    exit 1
fi

bias=$fottrack_tool/bias/${species}/nakedDNA.pickle
score_bwinbed=$fottrack_tool/TFOS/score_bwinbed.py
motif=$fottrack_tool/motif/${species}_merge_nocomplex

bed_dir=${ref_base}/cell/${cell}/tf/cor
ATAC=${ref_base}/cell/${cell}/atac/atac.bed
CTCF=${ref_base}/cell/${cell}/tf/cor/ctcf/ctcf_raw.bed
CTCF_flank=${ref_base}/cell/${cell}/tf/cor/ctcf/ctcf_flank.bed
genome=${ref_base}/genome.fa
chrom=${ref_base}/chrom.sizes
TSS=${ref_base}/gene.bed
TSS50=${ref_base}/tss_50-50.bed
snp=${ref_base}/cell/${cell}/snp/snp.bed
lambdaDNA=$ref_dir/lambdaDNA

##############################

cat $dir/sample.txt | while read i; do

	# 00.rename
	mkdir -p $dir/00.rename/fastqc
	cd $dir/00.rename;
	fastqc -t ${core} -o $dir/00.rename/fastqc ${i}*.fq.gz;
	echo "****************** fastqc Done! *******************"

	# 01.filter
	mkdir -p $dir/01.filter
	cd $dir/01.filter;
	trim_galore --trim-n --clip_R1 3 --clip_R2 10 --three_prime_clip_R1 3 --three_prime_clip_R2 3 -j 7 --length 35 -q 20 --fastqc --paired $dir/00.rename/${i}.R1.fq.gz $dir/00.rename/${i}.R2.fq.gz -o $dir/01.filter >> $dir/01.filter/trim.log 2>&1;
	echo "****************** 01.filter Done! *******************"

	# 02.align
	mkdir -p $dir/02.align
	cd $dir/02.align;
	basal -a $dir/01.filter/${i}.R1_val_1.fq.gz -b $dir/01.filter/${i}.R2_val_2.fq.gz -d ${genome} -m 1 -x 1000 -p ${core} -r 1 -v 0.06 -s 16 -S 1 -n 0 -g 1 -M C:T -o $dir/02.align/${i}.bam 2> $dir/02.align/${i}_BSMAP_report.txt
	samtools sort -@ ${core} $dir/02.align/${i}.bam > $dir/02.align/${i}.sort.bam;
	rm $dir/02.align/${i}.bam
	samtools flagstat $dir/02.align/${i}.sort.bam > $dir/02.align/${i}.sort.flagstat.txt;
	samtools index $dir/02.align/${i}.sort.bam;
	# mark duplicates
	picard MarkDuplicates --INPUT $dir/02.align/${i}.sort.bam --OUTPUT $dir/02.align/${i}_mkdup.bam --METRICS_FILE $dir/02.align/${i}.markdup.qc --ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT;
	rm ${i}_mkdup.bam
	dup=$(awk 'NR==8 {print $10}' ${i}.markdup.qc)
	echo "****************** 02.align Done! *******************"

	# lambdaDNA
	mkdir -p $dir/lambdaDNA
	cd $dir/lambdaDNA
	bismark  --parallel 4 --output_dir $dir/lambdaDNA --gzip --nucleotide_coverage $lambdaDNA -1 $dir/01.filter/${i}.R1_val_1.fq.gz -2 $dir/01.filter/${i}.R2_val_2.fq.gz
	total_reads_number=$(awk 'NR==7 { match($0, /[0-9]+/); print substr($0, RSTART, RLENGTH); }' ${i}.R1_val_1_bismark_bt2_PE_report.txt)
	lambdaDNA_number=$(awk 'NR==8 { match($0, /[0-9]+/); print substr($0, RSTART, RLENGTH); }' ${i}.R1_val_1_bismark_bt2_PE_report.txt)
	total_sites=$(awk 'NR==24 { match($0, /[0-9]+/); print substr($0, RSTART, RLENGTH); }' ${i}.R1_val_1_bismark_bt2_PE_report.txt)
	methylation_sites=$(awk 'NR>=31 && NR<=34 { match($0, /[0-9]+/); sum += substr($0, RSTART, RLENGTH); } END { print sum; }' ${i}.R1_val_1_bismark_bt2_PE_report.txt)
	echo "$i lambdaDNA ratio:" >> lambdaDNA_stat.txt
	lambda_DNA_ratio1=${lambdaDNA_number}/${total_reads_number}
	echo $lambda_DNA_ratio1 >> lambdaDNA_stat.txt
	lambda_DNA_ratio2=$(echo "scale=8; ${lambdaDNA_number}/${total_reads_number}*100" | bc | awk '{printf "%.6f%%", $0}')
	echo "$lambda_DNA_ratio2" >> lambdaDNA_stat.txt
	echo "$i lambdaDNA conversion rate:" >> lambdaDNA_stat.txt
	lambda_DNA_conversion_rates=$(echo "scale=8; ${methylation_sites} / ${total_sites}*100" | bc | awk '{printf "%.2f%%", $0}')
	echo "$lambda_DNA_conversion_rates" >> lambdaDNA_stat.txt
	rm $dir/01.filter/${i}.R1_val_1.fq.gz
	rm $dir/01.filter/${i}.R2_val_2.fq.gz
	echo "****************** lambdaDNA stat Done! *******************"

	# 03.methratio
	mkdir -p $dir/03.methratio
	cd $dir/03.methratio;
	python3 ${basalkit} avgmod $dir/02.align/${i}.sort.bam ${genome} -r -m 1 -o $dir/03.methratio/${i}_basal_methratio -i correct
	awk 'BEGIN{OFS="\t"}{if(NR > 1 && $5 != "NA") print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"1-$5"\t"$8+$10}' $dir/03.methratio/${i}_basal_methratio_AvgMod.tsv > $dir/03.methratio/${i}_methratio.tsv
	rm $dir/03.methratio/${i}_basal_methratio_AvgMod.tsv
	bedtools intersect -v -a $dir/03.methratio/${i}_methratio.tsv -b $snp > $dir/03.methratio/${i}_methratio.noSNP.tsv
	rm $dir/03.methratio/${i}_methratio.tsv
	echo "****************** 03.methratio Done! *******************"

	# 04.depthProfile
	mkdir -p $dir/04.depthProfile/ordinary $dir/04.depthProfile/normalize
	cd $dir/04.depthProfile
	cut -f 1,2,3,7 $dir/03.methratio/${i}_methratio.noSNP.tsv | awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/' > $dir/04.depthProfile/${i}_depth.noSNP.bedgraph
	# average depth
	mean_depth=$(datamash --header-in mean 4 < $dir/04.depthProfile/${i}_depth.noSNP.bedgraph)
	echo "$i average depth: " >> $dir/stat1.txt
	echo $mean_depth >> $dir/stat1.txt
	bedGraphToBigWig $dir/04.depthProfile/${i}_depth.noSNP.bedgraph $chrom $dir/04.depthProfile/ordinary/${i}.bw
	awk -v var="$mean_depth" 'BEGIN{OFS="\t"}{$4 = $4 / var; print}' $dir/04.depthProfile/${i}_depth.noSNP.bedgraph > $dir/04.depthProfile/${i}.bedgraph
	bedGraphToBigWig $dir/04.depthProfile/${i}.bedgraph $chrom $dir/04.depthProfile/normalize/${i}.bw
	rm $dir/04.depthProfile/${i}_depth.noSNP.bedgraph
	rm $dir/04.depthProfile/${i}.bedgraph
	echo "****************** 04.depthProfile Done! *******************"

	# 05.conversionProfile
	mkdir -p $dir/05.conversionProfile/ordinary $dir/05.conversionProfile/normalize
	cd $dir/05.conversionProfile
	cut -f 1,2,3,6 $dir/03.methratio/${i}_methratio.noSNP.tsv | awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/' > $dir/05.conversionProfile/${i}_conversion.noSNP.bedgraph
	# average conversion rate
	mean_conversion_rate=$(datamash --header-in mean 4 < "$dir/05.conversionProfile/${i}_conversion.noSNP.bedgraph")
	mean_cr_o=$(printf "%.2f%%" $(echo "$mean_conversion_rate * 100" | bc -l))
	echo "$i average conversion rate: " >> $dir/stat1.txt
	echo $mean_conversion_rate >> $dir/stat1.txt
	bedGraphToBigWig $dir/05.conversionProfile/${i}_conversion.noSNP.bedgraph $chrom $dir/05.conversionProfile/ordinary/${i}.bw
	awk -v mean_conversion_rate="$mean_conversion_rate" 'BEGIN{OFS="\t"}{$4 = $4 - mean_conversion_rate; print}' "${i}_conversion.noSNP.bedgraph" > "${i}_avg_nor.bedgraph"
	bedGraphToBigWig "$dir/05.conversionProfile/${i}_avg_nor.bedgraph" "$chrom" "$dir/05.conversionProfile/normalize/${i}.bw"
	rm $dir/05.conversionProfile/${i}_conversion.noSNP.bedgraph
	rm "$dir/05.conversionProfile/${i}_avg_nor.bedgraph"
	rm $dir/03.methratio/*.tsv
	echo "****************** 05.conversionProfile Done! *******************"

	# 06.FootTrack
	# 01.BiasCorrect
	mkdir -p $dir/06.FootTrack/01.BiasCorrect; cd $dir/06.FootTrack/01.BiasCorrect
	FootTrack BiasCorrect -b $dir/05.conversionProfile/ordinary/${i}.bw -g ${genome} --cores ${core} --norm-off --k_flank 10 --score_mat PWM --window 100 --prefix ${i} --outdir $dir/06.FootTrack/01.BiasCorrect --extend 0 --bias-pkl ${bias} --track-off bias expected_gl expected_lc corrected_lc >> ${i}_BiasCorrect.log 2>&1;
	echo "****************** ${i} 01.BiasCorrect done *******************"
	# 02.ScoreBigwig
	mkdir -p $dir/06.FootTrack/02.ScoreBigwig; cd $dir/06.FootTrack/02.ScoreBigwig
	FootTrack ScoreBigwig --signal $dir/06.FootTrack/01.BiasCorrect/${i}_corrected_gl.bw --regions $dir/06.FootTrack/01.BiasCorrect/${i}_effective_ranges.bed --output ${i}_ScoreBigWig.bw --fp-min 10 --fp-max 10 --flank-min 20 --flank-max 20 --extend 0 --cores ${core} >> ${i}_ScoreBigwig.log 2>&1;
	echo "****************** ${i} 02.ScoreBigwig done *******************"


done

foottrack_dir=$dir/06.FootTrack

# 03.BINDetect
mkdir -p $foottrack_dir/03.BINDetect;cd $foottrack_dir/03.BINDetect
bw=()
for s in "${sample[@]}"; do
    bw+=("$foottrack_dir/02.ScoreBigwig/${s}_ScoreBigWig.bw")
done
joined_bw=$(IFS=" "; echo "${bw[*]}")
FootTrack BINDetect --motifs ${motif}/*.jaspar --signals $joined_bw --genome ${genome} --peaks ${ATAC} --outdir $foottrack_dir/03.BINDetect --cores ${core} --bound-pvalue 5e-2 --naming name >> BINDetect.log 2>&1;
echo "************************************ 03.BINDetect done ********************************************"

# TFOS
tfos_dir=$foottrack_dir/TFOS
rawdata_dir=$foottrack_dir/01.BiasCorrect
suffix=_corrected_gl.bw
for s in ${sample[@]}; do
    (
    mkdir -p $tfos_dir/${s}; cd $tfos_dir/${s}
    echo ${s} > $tfos_dir/${s}_TFOS.txt
    echo ${s} > $tfos_dir/${s}_flank.txt
    echo ${s} > $tfos_dir/${s}_postive_ratio.txt
    mkdir motif left_flank right_flank TFOS TFOS_noblack flank
    cat ${bed_dir}/tf.txt | while read i; do
        python ${score_bwinbed} --bed $bed_dir/${i}/${i}_raw.bed --bw $rawdata_dir/${s}${suffix} >  motif/${i}_motif.bed
        python ${score_bwinbed} --bed $bed_dir/${i}/${i}_left_flank.bed --bw $rawdata_dir/${s}${suffix} > left_flank/${i}_left_flank.bed
        python ${score_bwinbed} --bed $bed_dir/${i}/${i}_right_flank.bed --bw $rawdata_dir/${s}${suffix} > right_flank/${i}_right_flank.bed
        paste motif/${i}_motif.bed left_flank/${i}_left_flank.bed right_flank/${i}_right_flank.bed | sort | uniq > ${i}_tmp.txt
        # TFOS
        awk 'BEGIN{OFS="\t"} {if ($4=="NA" || $8=="NA" || $12=="NA") print $1, $2, $3, "NA"; else print $1, $2, $3, ($8+$12)/2-$4}' ${i}_tmp.txt > ${i}_tmp2.bed
        sort -k1,1 -k2,2n ${i}_tmp2.bed > TFOS/${i}_TFOS.bed
        rm ${i}_tmp2.bed
        TFOS=$(awk '$4 != "NA" {print $0}' TFOS/${i}_TFOS.bed | datamash mean 4)
        echo ${TFOS} >> $tfos_dir/${s}_TFOS.txt
        # flank
        awk 'BEGIN{OFS="\t"} {if ($8=="NA" && $12=="NA") print $1, $2, $3, "NA"; else print $1, $2, $3, ($8+$12)/2}' ${i}_tmp.txt > ${i}_tmp2.bed
        sort -k1,1 -k2,2n ${i}_tmp2.bed > flank/${i}_flank.bed
        rm ${i}_tmp2.bed
        flank=$(awk '$4 != "NA" {print $0}' flank/${i}_flank.bed | datamash mean 4)
        echo ${flank} >> $tfos_dir/${s}_flank.txt
        # postive_ratio
        all_line=$(wc -l TFOS/${i}_TFOS.bed | awk '{print $1}')
        awk '$4 > 0 && $4 != "NA"' TFOS/${i}_TFOS.bed > TFOS_noblack/${i}_noblack.bed
        black_line=$(wc -l TFOS_noblack/${i}_noblack.bed | awk '{print $1}')
        proportions=$(echo "scale=4; ${black_line} / ${all_line}" | bc)
        echo ${proportions} >> $tfos_dir/${s}_postive_ratio.txt
        rm ${i}_tmp.txt
    done
    ) &
done
wait
cd $tfos_dir
echo "tf" > chip_TFOS.txt
echo "tf" > chip_flank.txt
echo "tf" > postive_ratio.txt
cat ${bed_dir}/tf.txt | while read i; do
    echo ${i} >> chip_TFOS.txt
    echo ${i} >> chip_flank.txt
    echo ${i} >> postive_ratio.txt
done
for s in "${sample[@]}"; do
    paste chip_TFOS.txt "${s}_TFOS.txt" > tmp.txt
    mv tmp.txt chip_TFOS.txt
    paste chip_flank.txt "${s}_flank.txt" > tmp.txt
    mv tmp.txt chip_flank.txt
    paste postive_ratio.txt "${s}_postive_ratio.txt" > tmp.txt
    mv tmp.txt postive_ratio.txt
done
rm -f tmp.txt