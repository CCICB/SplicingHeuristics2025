#!/bin/bash

## Script Definitions
gencode_file=../references/gencode.v44.annotation.gtf.gz
gencode_file_basename=$(basename $gencode_file | cut -f1,2 -d'.')

ref_genome=../references/hg38.fa

# Make Exon Map
gunzip -c $gencode_file | grep 'tag "Ensembl_canonical"' | grep -v "^chrM" | awk '($3=="exon")' | cut -f9 | cut -f2 -d';' | cut -f2 -d'"' | cut -f1 -d'.' | sort | uniq -c | tr -d ' ' | sed 's/E/\tE/g' | awk '{print $2, $1}' OFS='\t' | sort -k1 > $gencode_file_basename.exonmap.tsv

## Make exon bed
gunzip -c $gencode_file | grep 'tag "Ensembl_canonical"' | grep 'transcript_type "protein_coding' | grep -v "^chrM" | awk -F'\t' '($3=="exon") {print $1,$4,$5,$7,$9}' OFS='\t' | awk '{print $1,$2,$3,$4,$18,$8}' OFS='\t' | tr -d '";' | cut -f1 -d'.' | awk '{print $1,$2,$3,$6,$5,$4}' OFS='\t' | sort -k4 > $gencode_file_basename.exons.tsv
join -t $'\t' -1 4 -2 1 $gencode_file_basename.exons.tsv $gencode_file_basename.exonmap.tsv | awk '{print $2,$3,$4,$1,$5,$6,$7}' OFS='\'t | sort -k1,1 -k2,2n -k3,3n > $gencode_file_basename.exon_count.bed

## Get distance to nearby exons
# Ignore overlapping exons
bedtools closest -s -D a -io -id -t "all" -a $gencode_file_basename.exon_count.bed -b $gencode_file_basename.exon_count.bed | awk -F'\t' '($4==$11 || ($5=="1" && $7!="1")) {if ($6=="+") print $1,$2,$3,$4,$5,$6,$7,$8":"$10,$15; else print $1,$2,$3,$4,$5,$6,$7,$8":"$9,$15}' OFS='\t' | sort -u | sort -k1,1 -k2,2n -k3,3n > $gencode_file_basename.upstream_exon.bed
bedtools closest -s -D a -io -iu -t "all" -a $gencode_file_basename.upstream_exon.bed -b $gencode_file_basename.exon_count.bed | awk -F'\t' '($4==$13 || $5==$7) {if ($6=="+") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10":"$11,$17; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10":"$12,$17}' OFS='\t' | sort -u | sort -k1,1 -k2,2n -k3,3n > $gencode_file_basename.upstream_downstream_exon.bed

# ## Get distance to nearby U12 introns
# # Ignore overlapping U12 introns (removes exons that don't use the exact splice sites of the exon of interest)
# # Report only one entry if there is a tie
bedtools closest -s -D a -io -id -t "first" -a $gencode_file_basename.upstream_downstream_exon.bed -b ../references/U12.hg38.bed | awk '{if ($18=="-2") $18 = "U12"; else $18 = "U2"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$18,$10,$11}' OFS='\t' > $gencode_file_basename.upstream_U12.bed
bedtools closest -s -D a -io -iu -t "first" -a $gencode_file_basename.upstream_U12.bed -b U12.hg38.bed | awk '{if ($19=="2") $19 = "U12"; else $19 = "U2"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$19}' OFS='\t' > $gencode_file_basename.exons.canonical.proteincoding.annotated.tsv


runScript() {
    gencode_file=../references/gencode.v44.annotation.gtf.gz
    gencode_file_basename=$(basename $gencode_file | cut -f1,2 -d'.')

    ref_genome=../references/hg38.fa
    
    while read line; do
        ## Basic Info
        chr=$(echo "$line" | cut -f1)
        start=$(echo "$line" | cut -f2)
        end=$(echo "$line" | cut -f3)
        transcript=$(echo "$line" | cut -f4)
        strand=$(echo "$line" | cut -f6)

        ## Exon Info
        exon_number=$(echo "$line" | cut -f5)
        exon_count=$(echo "$line" | cut -f7)

        ## Nearby exon info
        upstream_exon_donor=$(echo "$line" | cut -f8 | cut -f2 -d':')
        upstream_intron_length=$(echo "$line" | cut -f9 | tr -d '-') # Reported as negative
        upstream_U12=$(echo "$line" | cut -f10)
        downstream_exon_acceptor=$(echo "$line" | cut -f11 | cut -f2 -d':')
        downstream_intron_length=$(echo "$line" | cut -f12)
        downstream_U12=$(echo "$line" | cut -f13)

        ## Set Exon Information based on exon location within gene
        if [ "$exon_number" == "1" ]; then
            # If more than one exon, get regions for donor only
            if [ "$exon_count" != "1" ]; then
                upstream_exon_donor="."
                upstream_intron_type="."
                upstream_intron_length="."
            else # Skip single exon genes
                continue
            fi
        elif [ "$exon_number" == "$exon_count" ]; then
            ## If last exon, get regions for acceptor only
            downstream_exon_acceptor="."
            downstream_intron_type="."
            downstream_intron_length="."
        fi

        exon_length=$(($end-$start+1))

        ## Get Intron and Exon Sequences
        if [ "$strand" == "+" ]; then
            acceptor_seq=$(samtools faidx $ref_genome $chr":"$(($start-50))"-"$(($start+50)) | grep -v "^>" | tr -d '\n')
            donor_seq=$(samtools faidx $ref_genome $chr":"$(($end-50))"-"$(($end+50)) | grep -v "^>" | tr -d '\n')
            
            if [ "$upstream_exon_donor" != "." ]; then
                upstream_exon_donor_seq=$(samtools faidx $ref_genome $chr":"$(($upstream_exon_donor-50))"-"$(($upstream_exon_donor+50)) | grep -v "^>" | tr -d '\n')
                upstream_intron=$(echo $chr":"$(($upstream_exon_donor+1))"-"$(($start-1)))
            else
                upstream_exon_donor_seq="."
                upstream_intron="."
                upstream_U12="."
            fi

            if [ "$downstream_exon_acceptor" != "." ]; then
                downstream_exon_acceptor_seq=$(samtools faidx $ref_genome $chr:$(($downstream_exon_acceptor-50))-$(($downstream_exon_acceptor+50)) | grep -v "^>" | tr -d '\n')
                downstream_intron=$(echo $chr":"$(($end+1))"-"$(($downstream_exon_acceptor-1)))
            else 
                downstream_exon_acceptor_seq="."
                downstream_intron="."
                downstream_U12="."
            fi
        else
            acceptor_seq=$(samtools faidx -i $ref_genome $chr":"$(($end-50))"-"$(($end+50)) | grep -v "^>" | tr -d '\n')
            donor_seq=$(samtools faidx -i $ref_genome $chr":"$(($start-50))"-"$(($start+50)) | grep -v "^>" | tr -d '\n')

            if [ "$upstream_exon_donor" != "." ]; then
                upstream_exon_donor_seq=$(samtools faidx -i $ref_genome $chr":"$(($upstream_exon_donor-50))"-"$(($upstream_exon_donor+50)) | grep -v "^>" | tr -d '\n')
                upstream_intron=$(echo $chr":"$(($end+1))"-"$(($upstream_exon_donor-1)))
            else 
                upstream_exon_donor_seq="."
                upstream_intron="."
                upstream_U12="."
            fi

            if [ "$downstream_exon_acceptor" != "." ]; then
                downstream_exon_acceptor_seq=$(samtools faidx -i $ref_genome $chr":"$(($downstream_exon_acceptor-50))"-"$(($downstream_exon_acceptor+50)) | grep -v "^>" | tr -d '\n')
                downstream_intron=$(echo $chr":"$(($downstream_exon_acceptor+1))"-"$(($start-1)))
            else 
                downstream_exon_acceptor_seq="."
                downstream_intron="."
                downstream_U12="."
            fi
        fi

        echo "$transcript"$'\t'"$strand"$'\t'"$exon_number"$'\t'"$exon_count"$'\t'"$upstream_exon_donor"$'\t'"$upstream_exon_donor_seq"$'\t'"$upstream_intron"$'\t'"$upstream_U12"$'\t'"$upstream_intron_length"$'\t'"$acceptor_seq"$'\t'"$chr:$start-$end"$'\t'"$exon_length"$'\t'"$donor_seq"$'\t'"$downstream_intron"$'\t'"$downstream_U12"$'\t'"$downstream_intron_length"$'\t'"$downstream_exon_acceptor"$'\t'"$downstream_exon_acceptor_seq" > /var/tmp/gencode/$chr:$start-$end-$upstream_exon-$downstream_exon.tsv
    done
}

export -f runScript

echo "transcript"$'\t'"strand"$'\t'"exon_number"$'\t'"exon_count"$'\t'"upstream_exon_donor_coords"$'\t'"upstream_exon_donor_seq"$'\t'"upstream_intron_coords"$'\t'"upstream_intron_type"$'\t'"upstream_intron_length"$'\t'"acceptor_seq"$'\t'"exon_coords"$'\t'"exon_length"$'\t'"donor_seq"$'\t'"downstream_intron_coords"$'\t'"downstream_intron_type"$'\t'"downstream_intron_length"$'\t'"downstream_exon_acceptor_coords"$'\t'"downstream_exon_acceptor_seq" > header.tsv

rm -r /var/tmp/gencode/
mkdir /var/tmp/gencode
parallel -a $gencode_file_basename.exons.canonical.proteincoding.annotated.tsv --pipepart runScript

find /var/tmp/gencode/. -type f -exec cat {} + > $gencode_file_basename.concat.tsv
cat header.tsv $gencode_file_basename.concat.tsv > $gencode_file_basename.exons.canonical.proteincoding.annotated.seq.tsv
rm $gencode_file_basename.concat.tsv