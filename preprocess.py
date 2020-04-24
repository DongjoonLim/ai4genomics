import pysam
import numpy
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

# chromosome = '3'
start = 10000
end = 30000
record_dict = SeqIO.to_dict(SeqIO.parse("data/GRCh37.p13.genome.fa", "fasta"))
# print(record_dict["chr1"].seq[30090])
        
# These are all children
samfile1 = pysam.AlignmentFile('data/MG17-5220.sorted.dedup.bam', "rb")
samfile2 = pysam.AlignmentFile('data/MG18-3089.sorted.dedup.bam', "rb")
samfile3 = pysam.AlignmentFile('data/MG18-3272.sorted.dedup.bam', "rb")
samfile4 = pysam.AlignmentFile('data/MG18-2009.sorted.dedup.bam', "rb")

# These are all parents
samfile5 = pysam.AlignmentFile('data/MG18-2010.sorted.dedup.bam', "rb")
samfile6 = pysam.AlignmentFile('data/MG18-2632.sorted.dedup.bam', "rb")
samfile7 = pysam.AlignmentFile('data/MG18-3013.sorted.dedup.bam', "rb")
samfile8 = pysam.AlignmentFile('data/MG18-3014.sorted.dedup.bam', "rb")
samfile9 = pysam.AlignmentFile('data/MG18-3088.sorted.dedup.bam', "rb")
samfile10 = pysam.AlignmentFile('data/MG18-3094.sorted.dedup.bam', "rb")
samfile11 = pysam.AlignmentFile('data/MG18-3270.sorted.dedup.bam', "rb")
samfile12 = pysam.AlignmentFile('data/MG18-3271.sorted.dedup.bam', "rb")

samfiles =[samfile1,samfile2,samfile3,samfile4,samfile5,samfile6,samfile7,samfile8,samfile9,samfile10,samfile11,samfile12]

def posDict(samfile, start, end, chromosome):
    result = {}
    for pileupcolumn in tqdm(samfile.pileup(chromosome)):
        
        # print ("/ncoverage at base %s = %s" %
        #     (pileupcolumn.pos, pileupcolumn.n))
        nucList = []
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                # print ('/tbase in read %s = %s' %
                #     (pileupread.alignment.query_name,
                #     pileupread.alignment.query_sequence[pileupread.query_position]))
                nucList.append(pileupread.alignment.query_sequence[pileupread.query_position])
        # print(nucList)
        try:
            result[pileupcolumn.pos] = max(set(nucList), key = nucList.count)
        except ValueError:
            continue
    return result

def genDictList(start, end, chromosome):
    dictList = []
    for i in range(12):
        dictList.append(posDict(samfiles[i], start, end, chromosome))
    return dictList



def genData(raw, start, end, chromosome):
    posList = set(list(raw[0].keys()))
    for i in range(1, 12):
        posList =posList.intersection(list(raw[i].keys()))
    posList = sorted(posList)
    data = {}
    print('This is posList : ', posList)
    majority = []
    for item in tqdm(posList):
        data[item] = []
        for sample in raw:
            if sample[item] == str(record_dict['chr{}'.format(chromosome)].seq[item]):
                data[item].append(0)
            else:
                data[item].append(1)
        if len(set(data[item])) == 1:
            del data[item]
    df = pd.DataFrame.from_dict(data)
    df.to_csv("{}_{}_{}.csv".format('chr{}'.format(chromosome), start, end))

for chro in range(1,23):
    chromosome = str(chro)
    result = genDictList(start,end, chromosome)
    print(result)
    genData(result, start, end, chromosome)
    
samfile1.close()
samfile2.close()
samfile3.close()
samfile4.close()
samfile5.close()
samfile6.close()
samfile7.close()
samfile8.close()
samfile9.close()
samfile10.close()
samfile11.close()
samfile12.close()