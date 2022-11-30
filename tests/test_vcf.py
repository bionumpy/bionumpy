# from bionumpy.vcf import parse_header, parse_info_line
import pytest


lines = """\
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=27022019_15h52m43s
##source=IGSRpipeline
##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">
##contig=<ID=10>
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
##INFO=<ID=EX_TARGET,Number=0,Type=Flag,Description="indicates whether a variant is within the exon pull down target boundaries">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##bcftools_viewVersion=1.4.1+htslib-1.4.1
##bcftools_viewCommand=view -q 0.001 -O z ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz; Date=Fri Jan 22 15:41:31 2021
##contig=<ID=1>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=2>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=X>
##bcftools_concatVersion=1.4.1+htslib-1.4.1
##bcftools_concatCommand=concat -o all_variants_with_genotypes.vcf.gz -O z variants_with_genotypes_chr10_0.001.vcf.gz variants_with_genotypes_chr1_0.001.vcf.gz variants_with_genotypes_chr11_0.001.vcf.gz variants_with_genotypes_chr12_0.001.vcf.gz variants_with_genotypes_chr13_0.001.vcf.gz variants_with_genotypes_chr14_0.001.vcf.gz variants_with_genotypes_chr15_0.001.vcf.gz variants_with_genotypes_chr16_0.001.vcf.gz variants_with_genotypes_chr17_0.001.vcf.gz variants_with_genotypes_chr18_0.001.vcf.gz variants_with_genotypes_chr19_0.001.vcf.gz variants_with_genotypes_chr20_0.001.vcf.gz variants_with_genotypes_chr2_0.001.vcf.gz variants_with_genotypes_chr21_0.001.vcf.gz variants_with_genotypes_chr22_0.001.vcf.gz variants_with_genotypes_chr3_0.001.vcf.gz variants_with_genotypes_chr4_0.001.vcf.gz variants_with_genotypes_chr5_0.001.vcf.gz variants_with_genotypes_chr6_0.001.vcf.gz variants_with_genotypes_chr7_0.001.vcf.gz variants_with_genotypes_chr8_0.001.vcf.gz variants_with_genotypes_chr9_0.001.vcf.gz variants_with_genotypes_chrX_0.001.vcf.gz; Date=Fri Jan 22 16:54:29 2021
##bcftools_viewCommand=view -O z --regions 1:1-5000000 data/variants.vcf.gz; Date=Mon May 30 20:18:21 2022\
""".split("\n")

info_line = '''##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">'''


@pytest.mark.skip
def test_parse_info():
    assert parse_info_line(info_line) == {"ID": "AF", "Number": "A", "Type": "Float"}

@pytest.mark.skip("unimplemented")
def test_parse_header():
    assert parse_header(lines) == {"ID": "AF", "Number": "A", "Type": "Float"}
    
