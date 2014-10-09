import vcf
import sys


def run(filepath):
    vcf_reader = vcf.Reader(open(filepath, 'rb'))

    for rec in vcf_reader:
        print 'chr'+rec.CHROM+'\t'+str(rec.POS)+'\t'+str(rec.POS+len(rec.REF))

if __name__ == '__main__':
    run(sys.argv[1])
