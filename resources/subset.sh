bcftools view -s NA21110 -m2 -M2 \
shapeit5-phased-callset_final-vcf.phased.vcf.gz \
| bcftools view -i 'GT!="RR" && GT!="mis"' -Oz -o  NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz