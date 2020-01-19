#!/usr/bin/env bash
# gtf_transcript_gene_map.sh
# Kamil Slowikowski
# November 26, 2015
#
# Create a map of transcripts to genes from a GTF file.

# if [[ ! -e "$1" ]]; then
#   echo "usage: ./$0 in.gtf > out.tsv"
#   exit 1
# fi

perl -lne '
BEGIN {
print "transcript_id\ttranscript_name\tgene_id\tgene_name"
}
($t) = /transcript_id "([^"]+)"/;
($tn) = /transcript_name "([^"]+)"/;
($g) = /gene_id "([^"]+)"/;
($gn) = /gene_name "([^"]+)"/;
if ($t && $tn && $g && $gn && !$H{$t}++) {
print "$t\t$tn\t$g\t$gn"
}' <(gunzip -c Homo_sapiens.GRCh38.89.gtf.gz) | gzip > Homo_sapiens.GRCh38.89.gene_names.tsv.gz

