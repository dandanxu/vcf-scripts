import vcf
import sys
import time

MOTHER_SAMPLE = 9
FATHER_SAMPLE = 10
TOTAL_AD = 100


def de_novo_one_parent(progeny, parent):
    for n in progeny:
        if n in parent:
            return False
    return True


def de_novo_both_parents(progeny, parents):
    possible_progeny = [[parents[0][0], parents[1][0]],
                        [parents[0][0], parents[1][1]],
                        [parents[0][1], parents[1][0]],
                        [parents[0][1], parents[1][1]]]
    if (progeny in possible_progeny) or (progeny[::-1] in possible_progeny):
        return False
    return True


def run(filepath):
    # speed testing
    t0 = time.time()
    n_records = 0

    # opens file with PyVCF
    vcf_reader = vcf.Reader(open(filepath, 'rb'))
    output = 'output_de_novo_100ad_filter.vcf'
    vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)

    SAMPLES = ['701-1-1', '701-1-3', '701-1-6',
               '701-1-7', '701-2-1', '701-2-2',
               '701-2-3', '701-2-4', '701-3-2',
               '701-Female', '701-Male']
    de_novo_counts = {sample: 0 for sample in SAMPLES}

    for record in vcf_reader:
        n_records += 1
        match_type = 0
        de_novo = False

        mother = record.samples[MOTHER_SAMPLE]
        father = record.samples[FATHER_SAMPLE]

        if mother['GT'] is None:
            if father['GT'] is None:
                match_type = 0
            else:
                one_parent = [father['GT'][0], father['GT'][2]]
                # print 'father:', one_parent
                match_type = 1
        else:
            if father['GT'] is None:
                one_parent = [mother['GT'][0], mother['GT'][2]]
                # print 'mother:', one_parent
                match_type = 1
            else:
                both_parents = [[mother['GT'][0], mother['GT'][2]],
                                [father['GT'][0], father['GT'][2]]]
                # print 'both:', both_parents
                match_type = 2

        try:
            if sum(mother['AD']) < TOTAL_AD:
                match_type = 0
        except:
            pass

        try:
            if sum(father['AD']) < TOTAL_AD:
                match_type = 0
        except:
            pass

        if match_type != 0:
            for sample in range(0, 9):
                progeny = record.samples[sample]
                if progeny['GT'] is not None:
                    progeny_gt = [progeny['GT'][0], progeny['GT'][2]]
                try:
                    if sum(progeny['AD']) > TOTAL_AD:
                        if match_type == 1:
                            if de_novo_one_parent(progeny_gt, one_parent):
                                de_novo_counts[progeny.sample] += 1
                                de_novo = True
                                print progeny, one_parent, mother['AD'], father['AD']
                        elif match_type == 2:
                            if de_novo_both_parents(progeny_gt, both_parents):
                                de_novo = True
                                de_novo_counts[progeny.sample] += 1
                                print progeny, both_parents, mother['AD'], father['AD']
                        else:
                            pass
                except:
                    pass
        if de_novo:
            vcf_writer.write_record(record)

        if n_records % 10000 == 0:
            print 'Processed', n_records, 'records...'
            # break
    print de_novo_counts, 'total records:', n_records
    t1 = time.time()
    print 'Buona notte!', t1-t0, 'seconds have passed.'

if __name__ == '__main__':
    run(sys.argv[1])
