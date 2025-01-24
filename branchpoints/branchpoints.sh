#!/bin/bash

#gencode_file=../references/gencode.v44.annotation.gtf.gz # hg38
gencode_file=../references/gencode.v44lift37.annotation.gtf.gz # hg19
gencode_file_basename=$(basename $gencode_file | cut -f1,2 -d'.')

# Make bed file from Supplementary Data 3
gunzip -c supp_gr.182899.114_Supplemental_DataS3.txt.gz | awk '{print $7,$7,$6,$5,$4}' OFS='\t' | sed -e 's/_/\t/1' -e 's/_/\t/1' -e 's/_.//1' | cut -f1 -d'.' > bp_mercer.hg19.bed
bedtools sort -i bp_mercer.hg19.bed > bp_mercer.hg19.sort.bed
mv bp_mercer.hg19.sort.bed bp_mercer.hg19.bed

# # Convert coordinates to hg38
# liftOver -bedPlus=3 bp_mercer.hg19.bed ~/Tools/hg19ToHg38.over.chain.gz bp_mercer.hg38.bed bp_mercer.hg38.unmapped.bed
# bedtools sort -i bp_mercer.hg38.bed > bp_mercer.hg38.sort.bed
# mv bp_mercer.hg38.sort.bed bp_mercer.hg38.bed

# Make Exon Map
gunzip -c $gencode_file | grep -v "^chrM" | awk '($3=="exon")' | cut -f9 | cut -f2 -d';' | cut -f2 -d'"' | cut -f1 -d'.' | sort | uniq -c | tr -d ' ' | sed 's/E/\tE/g' | awk '{print $2, $1}' OFS='\t' | sort -k1 > $gencode_file_basename.exonmap.tsv

## Make exon BED 
gunzip -c $gencode_file | grep -v "^chrM" | awk -F'\t' '($3=="exon") {print $1,$4,$5,$7,$9}' OFS='\t' | awk '{print $1,$2,$3,$4,$18,$8}' OFS='\t' | tr -d '";' | cut -f1 -d'.' | awk '{print $1,$2,$3,$6,$5,$4}' OFS='\t' | sort -k4 > $gencode_file_basename.exons.tsv
join -t $'\t' -1 4 -2 1 $gencode_file_basename.exons.tsv $gencode_file_basename.exonmap.tsv | awk '($7!="1") {print $2,$3,$4,$1,$5,$6}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n > $gencode_file_basename.exon_count.bed
awk -F'\t' '($5!="1")' OFS='\t' $gencode_file_basename.exon_count.bed | sort -k1,1 -k2,2n -k3,3n > $gencode_file_basename.exon_count.downstream.bed

## Get closest exon
bedtools closest -s -D a -io -iu -t "all" -a bp_mercer.hg19.bed -b $gencode_file_basename.exon_count.downstream.bed | awk -F'\t' '($7==$11) {if ($6=="+") print $1,$2,$3,$4,$5,$6,$7,$8":"$9,$14; else print $1,$2,$3,$4,$5,$6,$7,$8":"$10,$14}' OFS='\t' | sort -u | sort -k1,1 -k2,2n -k3,3n > bp_mercer.hg19.downstream_exon.bed
bedtools closest -s -D a -io -id -t "all" -a bp_mercer.hg19.downstream_exon.bed -b $gencode_file_basename.exon_count.bed | awk -F'\t' '($7==$13) {if ($6=="+") print $1,$2,$3,$4,$5,$6,"0",$10":"$12,$16,$8,$9; else print $1,$2,$3,$4,$5,$6,"0",$10":"$11,$16,$8,$9}' OFS='\t' | sort -u | sort -k1,1 -k2,2n -k3,3n > bp_mercer.hg19.downstream_upstream_exon.bed

cat bp_mercer.hg19.downstream_upstream_exon.bed | awk -F'\t' '{print $10,$4,$11,$6,$0}' OFS='\t' | sed 's/:/\t/1' | awk '{print $1,$2,$2,$0}' OFS='\t' | cut -f1-3,6- | sort -k1,1 -k2,2n -k3,3n > downstream_exon.bed
bedtools closest -s -D a -io -iu -t "first" -a downstream_exon.bed -b $gencode_file_basename.exon_count.downstream.bed | cut -f7-20,24 | awk -F'\t' '{if ($6=="+") print $12,$13,$13,$4,$15,$6,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12":"$13,$15; else print $12,$14,$14,$4,$15,$6,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12":"$14,$15}' OFS='\t' | sort -u | sort -k1,1 -k2,2n -k3,3n > downstream_exon2.bed
bedtools closest -s -D a -io -iu -t "first" -a downstream_exon2.bed -b $gencode_file_basename.exon_count.downstream.bed | cut -f7-22,26 | awk -F'\t' '{if ($6=="+") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14":"$15,$17; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14":"$16,$17}' OFS='\t' | sort -u | sort -k1,1 -k2,2n -k3,3n > bp_mercer.hg19.downstream3_upstream_exon.bed

bedtools intersect -c -a bp_mercer.hg19.downstream3_upstream_exon.bed -b ../references/U12.hg19.bed > bp_mercer.hg19.annotated.tsv
grep -vf fixes.tsv bp_mercer.hg19.annotated.tsv > bp_mercer.hg19.annotated2.tsv 
mv bp_mercer.hg19.annotated2.tsv bp_mercer.hg19.annotated.tsv 

runBPScript() {
    ref_genome=~/Projects/References/hg19.fa
    
    while read line; do
        chr=$(echo "$line" | cut -f1)
        bp_position=$(echo "$line" | cut -f3)
        bp_options=$(echo "$line" | cut -f5)
        strand=$(echo "$line" | cut -f6)
        
        # Surrounding exons
        upstream_exon_donor=$(echo "$line" | cut -f8 | cut -f2 -d ':')
        upstream_exon_distance=$(echo "$line" | cut -f9 | tr -d '-') # Reported as negative
        downstream_exon_acceptor=$(echo "$line" | cut -f10 | cut -f2 -d ':')
        downstream_exon_distance=$(echo "$line" | cut -f11)
        downstream_exon_acceptor2=$(echo "$line" | cut -f12 | cut -f2 -d ':')
        downstream_exon_distance2=$(echo "$line" | cut -f13)
        downstream_exon_acceptor3=$(echo "$line" | cut -f14 | cut -f2 -d ':')
        downstream_exon_distance3=$(echo "$line" | cut -f15)

        if [ $strand == "+" ]; then
            downstream_exon_distance=$(($downstream_exon_distance-1))
        fi

        ## U12 info
        U12_overlaps=$(echo "$line" | cut -f16)

        if [ "$U12_overlaps" == "0" ]; then
            intron_type="U2"
        else
            intron_type="U12"
        fi

        if [ "$strand" == "+" ]; then
            bp_seq=$(samtools faidx $ref_genome $chr:$(($bp_position-50))-$(($bp_position+50)) | grep -v "^>" | tr -d '\n')
            upstream_donor_seq=$(samtools faidx $ref_genome $chr:$(($upstream_exon_donor-50))"-"$(($upstream_exon_donor+50)) | grep -v "^>" | tr -d '\n')
            intron=$(echo $chr":"$(($upstream_exon_donor+1))"-"$(($downstream_exon_acceptor-1)))
            intron_seq=$(samtools faidx $ref_genome $intron | grep -v "^>" | tr -d '\n')
            acceptor_seq=$(samtools faidx $ref_genome $chr:$(($downstream_exon_acceptor-50))"-"$(($downstream_exon_acceptor+50)) | grep -v "^>" | tr -d '\n')
        else
            bp_seq=$(samtools faidx -i $ref_genome $chr:$(($bp_position-50))-$(($bp_position+50)) | grep -v "^>" | tr -d '\n')
            upstream_donor_seq=$(samtools faidx -i $ref_genome $chr:$(($upstream_exon_donor-50))"-"$(($upstream_exon_donor+50)) | grep -v "^>" | tr -d '\n')
            intron=$(echo $chr":"$(($downstream_exon_acceptor+1))"-"$(($upstream_exon_donor-1)))
            intron_seq=$(samtools faidx -i $ref_genome $intron | grep -v "^>" | tr -d '\n')
            acceptor_seq=$(samtools faidx -i $ref_genome $chr:$(($downstream_exon_acceptor-50))"-"$(($downstream_exon_acceptor+50)) | grep -v "^>" | tr -d '\n')
        fi

        if (( $((100 - $downstream_exon_distance - $downstream_exon_distance2)) > 0 )); then
            acceptor_distance2=$(($downstream_exon_distance + $downstream_exon_distance2))
            if [ $strand == "+" ]; then
                acceptor_seq2=$(samtools faidx $ref_genome $chr:$(($downstream_exon_acceptor2-50))"-"$(($downstream_exon_acceptor2+50)) | grep -v "^>" | tr -d '\n')
            else
                acceptor_seq2=$(samtools faidx -i $ref_genome $chr:$(($downstream_exon_acceptor2-50))"-"$(($downstream_exon_acceptor2+50)) | grep -v "^>" | tr -d '\n')
            fi
            
            if (( $((100 - $downstream_exon_distance - $downstream_exon_distance2 - $downstream_exon_distance3 )) > 0 )); then
                acceptor_distance3=$(($downstream_exon_distance + $downstream_exon_distance2 + $downstream_exon_distance3))
                if [ $strand == "+" ]; then
                    acceptor_seq3=$(samtools faidx $ref_genome $chr:$(($downstream_exon_acceptor3-50))"-"$(($downstream_exon_acceptor3+50)) | grep -v "^>" | tr -d '\n')
                else
                    acceptor_seq3=$(samtools faidx -i $ref_genome $chr:$(($downstream_exon_acceptor3-50))"-"$(($downstream_exon_acceptor3+50)) | grep -v "^>" | tr -d '\n')
                fi
            else
                acceptor_distance3="."
                acceptor_seq3="."
            fi
        else
            acceptor_distance2="."
            acceptor_seq2="."
            acceptor_distance3="."
            acceptor_seq3="."
        fi

        intron_length=$(echo $(( $(echo "$intron" | cut -f2 -d':')-1 )) | tr -d '-')
           
        echo "$chr:$bp_position"$'\t'"$bp_options"$'\t'"$strand"$'\t'"$bp_seq"$'\t'"$upstream_exon_donor"$'\t'"$upstream_exon_distance"$'\t'"$upstream_donor_seq"$'\t'"$intron"$'\t'"$intron_type"$'\t'"$intron_length"$'\t'"$downstream_exon_acceptor"$'\t'"$downstream_exon_distance"$'\t'"$acceptor_seq"$'\t'"$downstream_exon_acceptor2"$'\t'"$acceptor_distance2"$'\t'"$acceptor_seq2"$'\t'"$downstream_exon_acceptor3"$'\t'"$acceptor_distance3"$'\t'"$acceptor_seq3" > /var/tmp/branchpoints/$chr:$bp_position-$downstream_exon_acceptor.tsv
    done
}

export -f runBPScript

echo "coordinate"$'\t'"bp_options"$'\t'"strand"$'\t'"bp_seq"$'\t'"upstream_exon_donor"$'\t'"upstream_exon_distance"$'\t'"upstream_donor_seq"$'\t'"intron"$'\t'"intron_type"$'\t'"intron_length"$'\t'"downstream_exon_acceptor"$'\t'"downstream_exon_acceptor_distance"$'\t'"acceptor_seq"$'\t'"downstream_exon_acceptor2"$'\t'"downstream_exon_acceptor_distance2"$'\t'"acceptor_seq2"$'\t'"downstream_exon_acceptor3"$'\t'"downstream_exon_acceptor_distance3"$'\t'"acceptor_seq3" > branchpoints.header.tsv

rm -r /var/tmp/branchpoints/
mkdir /var/tmp/branchpoints
parallel -a bp_mercer.hg19.annotated.tsv --pipepart runBPScript

find /var/tmp/branchpoints/. -type f -exec cat {} + | sort -u > branchpoints.concat.tsv
cat branchpoints.header.tsv branchpoints.concat.tsv > branchpoints.hg19.fastas.tsv
rm branchpoints.concat.tsv