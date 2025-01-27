import pyBigWig
import pybedtools
import argparse

parser = argparse.ArgumentParser(description="Calculate the mean value of BigWig files over regions specified in a BED file.")
parser.add_argument("--bw", type=str, required=True, help="Path to the BigWig file.")
parser.add_argument("--bed", type=str, required=True, help="Path to the BED file.")
args = parser.parse_args()

bw = pyBigWig.open(args.bw)
bed = pybedtools.BedTool(args.bed)

results = []
for region in bed:
    start = max(0, region.start)
    end = max(start + 1, region.end)

    try:
        mean_val = bw.stats(region.chrom, start, end, type='mean')[0]
        if mean_val is not None:
            results.append(f"{region.chrom}\t{start}\t{end}\t{mean_val:.8f}")
        else:
            results.append(f"{region.chrom}\t{start}\t{end}\tNA")
    except RuntimeError:
        results.append(f"{region.chrom}\t{start}\t{end}\tNA")
bw.close()

for result in results:
    print(result)