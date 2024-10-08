import glob
import os
import argparse
import sys
import pandas as pd
import subprocess
from pybedtools import BedTool
from foottrack.parsers import add_submerge_arguments
from foottrack.utils.logger import foottrackLogger


def run_submerge(args):

    logger = foottrackLogger("SubMerge", args.verbosity)
    logger.begin()
    parser = add_submerge_arguments(argparse.ArgumentParser())
    logger.arguments_overview(parser, args)

    # Get all TFBS files
    tfbs_files = glob.glob(os.path.join(args.tfbs, "*/*_overview.txt"))

    # check if all required columns are in regions file
    regions_obj = BedTool(args.regions)
    regions_columns = regions_obj.field_count()

    if regions_columns < 6:

        # convert to dataframe
        regions_df = regions_obj.to_dataframe()

        # add missing columns
        while regions_columns < 6:

            match regions_columns:
                case 3:
                    regions_df.insert(3, 'name', 'region')  # add identifier
                case 4:
                    regions_df.insert(4, 'score', '.')      # add score placeholder
                case 5:
                    regions_df.insert(5, 'strand', '.')     # add strand placeholder

            regions_columns = len(regions_df.columns)

    # convert df back to bedtool
    regions_obj = BedTool.from_dataframe(regions_df)

    # Read TFs
    if args.tf is not None:
        logger.debug("Reading TFs from: " + args.tf)
        tfs = pd.read_csv(args.tf, header=None, names=["TF"])
        tfs = tfs["TF"].tolist()

        # get only those TFs that were provided
        tfbs_files = [tfbs for tfbs in tfbs_files if os.path.basename(os.path.dirname(tfbs)) in tfs]

    # intersect each TFBS file with the regions
    logger.debug('Intersecting query with all TBFS files')
    all_intersections = ''
    for file in tfbs_files:

        # remove header of TF file and intersect with regions
        command = f'sed "1d" {file} | bedtools intersect -a {regions_obj.fn} -b stdin -wa -wb'
        intersection = subprocess.check_output(command, shell=True).decode("utf-8")
        all_intersections += intersection

    # make list out of lines
    all_intersections = all_intersections.strip().split("\n")

    # get tfbs header
    with open(tfbs_files[0]) as f:
        tfbs_header = f.readline()

    # merge all intersection files
    logger.debug('Merging all intersections')
    with open(args.output, "w") as f:

        query_header = "query chr\tquery start\tquery end\tquery id\tquery score\tquery strand\t"
        header = query_header + tfbs_header
        f.write(header)

        for line in all_intersections:
            f.write(line + "\n")

    # filter output file
    df = pd.read_csv(args.output, sep="\t")

    logger.debug('Sorting')

    # Sort the dataframe
    df.sort_values(by=["query chr", "query start", "TFBS_name", "TFBS_chr", "TFBS_start"],
                   key=lambda x: x.str.replace("chr", "").astype(int) if x.name in ["query chr", "TFBS_chr"] else x,
                   inplace=True)

    if args.output.endswith(".xlsx"):
        df.to_excel(args.output, index=False)
    else:
        df.to_csv(args.output,
                  sep="\t", 
                  index=False,
                  header=not args.output.endswith(".bed"))

    logger.info("Output written to: " + args.output)
    logger.end()


def main():

    parser = argparse.ArgumentParser()
    parser = add_submerge_arguments(parser)
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    run_submerge(args)


if __name__ == '__main__':
    main()
