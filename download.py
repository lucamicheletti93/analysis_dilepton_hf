"""
Simple script to download and merge train outputs
"""

import argparse
import os
import numpy as np


def merge(nfiles, copy_ao2ds, outputdir, outsuffix, merge_ao2ds):
    """
    main function to merge files downloaded from alien
    """
    if nfiles > 1 and merge_ao2ds:
        os.system(f"hadd -f -n 100 {os.path.join(outputdir, f'AnalysisResults{outsuffix}.root')} "
                  "AnalysisResults_*.root")
        os.system("rm AnalysisResults_*.root")
        if copy_ao2ds:
            if nfiles < 200:
                os.system("o2-aod-merger --input ao2ds_to_merge.txt --max-size 10000000000"
                          f" --output {os.path.join(outputdir, f'AO2D{outsuffix}.root')}")
                os.system("rm ao2ds_to_merge.txt")
                os.system("rm AO2D_*.root")
            else: # do partial merging
                nchunks = int(np.ceil(nfiles / 200))
                os.system(f"mkdir {os.path.join(outputdir, f'AO2D{outsuffix}')}")
                for ichunk in range(nchunks):
                    row_min = ichunk * 200
                    row_max = (ichunk + 1) * 200
                    os.system(f"awk 'NR>={row_min} && NR<{row_max}' ao2ds_to_merge.txt"
                              f" > ao2ds_to_merge_chunk{ichunk}.txt")
                    subdir = os.path.join(outputdir, f"AO2D{outsuffix}", f"{ichunk:03d}")
                    os.system(f"mkdir {subdir}")
                    os.system(f"o2-aod-merger --input ao2ds_to_merge_chunk{ichunk}.txt --max-size 10000000000"
                              f" --output {os.path.join(subdir, f'AO2D{outsuffix}.root')}")
    else:
        if nfiles == 1:
            os.system(f"mv AnalysisResults_0.root {os.path.join(outputdir, f'AnalysisResults{outsuffix}.root')}")
            os.system(f"mv AO2D_0.root {os.path.join(outputdir, f'AO2D{outsuffix}.root')}")
        else:
            # we anyway merge AnalysisResults.root files
            os.system(f"hadd -f -n 100 {os.path.join(outputdir, f'AnalysisResults{outsuffix}.root')} "
                      "AnalysisResults_*.root")
            os.system("rm AnalysisResults_*.root")
            os.system(f"mkdir {os.path.join(outputdir, f'AO2D{outsuffix}')}")
            for ifile in range(nfiles):
                subdir = os.path.join(outputdir, f"AO2D{outsuffix}", f"{ifile:03d}")
                os.system(f"mkdir {subdir}")
                os.system(f"mv AO2D_{ifile}.root {subdir}")


def download(infile, copy_ao2ds, is_unmerged, outputdir, outsuffix, merge_ao2ds):
    """
    main function to download files from alien and merge them
    """
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    os.system(f"rm {os.path.join(outputdir, f'AnalysisResults{outsuffix}.root')}")
    os.system(f"rm {os.path.join(outputdir, f'AO2D{outsuffix}.root')}")

    nfiles = 0
    with open(infile, "r") as file:
        if not is_unmerged:
            for iline, line in enumerate(file):
                nfiles += 1
                line = line.replace('\n', '')
                os.system(f"alien_cp alien://{os.path.join(line, 'AnalysisResults.root')} "
                          f"file:AnalysisResults_{iline}.root")
                if copy_ao2ds:
                    os.system(f"alien_cp alien://{os.path.join(line, 'AO2D.root')} "
                              f"file:AO2D_{iline}.root")
            if copy_ao2ds:
                os.system("ls AO2D_*.root > ao2ds_to_merge.txt")
        else:
            ifile_analysis, ifile_ao2d = 0, 0
            for iline, line in enumerate(file):
                line = line.replace('\n', '')
                os.system(f"alien_find alien://{line} AnalysisResults.root > infile_analysisresults.txt")
                with open("infile_analysisresults.txt", "r") as file2:
                    for line2 in file2:
                        nfiles += 1
                        line2 = line2.replace('\n', '')
                        os.system(f"alien_cp alien://{line2} file:AnalysisResults_{ifile_analysis}.root")
                        ifile_analysis += 1
                os.system("rm infile_analysisresults.txt")
                os.system("ls AnalysisResults_*.root > analysisresults_to_merge.txt")
                if copy_ao2ds:
                    os.system(f"alien_find alien://{line} AO2D.root > infile_ao2ds.txt")
                    with open("infile_ao2ds.txt", "r") as file2:
                        for line2 in file2:
                            line2 = line2.replace('\n', '')
                            os.system(f"alien_cp alien://{line2} file:AO2D_{ifile_ao2d}.root")
                            ifile_ao2d += 1
                    os.system("rm infile_ao2ds.txt")
                    os.system("ls AO2D_*.root > ao2ds_to_merge.txt")

    merge(nfiles, copy_ao2ds, outputdir, outsuffix, merge_ao2ds)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("--infile", "-i", metavar="text",
                        default="files.txt", help="input file with list of paths (w/o file name)",
                        required=True)
    PARSER.add_argument("--aod2s", "-a", action="store_true", default=False,
                        help="copy also AO2D.root files", required=False)
    PARSER.add_argument("--unmerged", "-u", action="store_true", default=False,
                        help="flag for unmerged train outputs", required=False)
    PARSER.add_argument("--outputdir", "-o", metavar="text", default=".",
                        help="output directory", required=False)
    PARSER.add_argument("--outputsuffix", "-s", metavar="text", default="",
                        help="output suffix", required=False)
    PARSER.add_argument("--merge_ao2ds", "-m", action="store_true", default=False,
                        help="merge downloaded ao2ds", required=False)
    ARGS = PARSER.parse_args()

    download(ARGS.infile, ARGS.aod2s, ARGS.unmerged,
             ARGS.outputdir, ARGS.outputsuffix, ARGS.merge_ao2ds)
