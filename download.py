"""
Simple script to download and merge train outputs
"""

import argparse
import os

def download(infile, copy_ao2ds=False, is_unmerged=False, outputdir=".", outsuffix=""):
    """
    main function to download files from alien and merge them
    """
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    os.system(f"rm {os.path.join(outputdir, f'AnalysisResults_{outsuffix}.root')}")
    os.system(f"rm {os.path.join(outputdir, f'AO2D_{outsuffix}.root')}")

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
        else:
            for iline, line in enumerate(file):
                os.system(f"alien_find alien://{line} AnalysisResults.root > infile_analysisresults.txt")
                with open("infile_analysisresults.txt", "r") as file2:
                    for iline2, line2 in enumerate(file2):
                        nfiles += 1
                        line2 = line2.replace('\n', '')
                        os.system(f"alien_cp alien://{line2} file:AnalysisResults_{iline2}.root")
                os.system("rm infile_analysisresults.txt")
                if copy_ao2ds:
                    os.system(f"alien_find alien://{line} AO2D.root > infile_ao2ds.txt")
                    with open("infile_ao2ds.txt", "r") as file2:
                        for iline2, line2 in enumerate(file2):
                            line2 = line2.replace('\n', '')
                            os.system(f"alien_cp alien://{line2} file:AO2D_{iline2}.root")
                    os.system("rm infile_ao2ds.txt")
                    os.system("ls AO2D_*.root > ao2ds_to_merge.txt")

    if nfiles > 1:
        os.system(f"hadd -f {os.path.join(outputdir, f'AnalysisResults{outsuffix}.root')} AnalysisResults_*.root")
        os.system("rm AnalysisResults_*.root")
        if copy_ao2ds:
            os.system("o2-aod-merger --input ao2ds_to_merge.txt --max-size 10000000000"
                    f" --output {os.path.join(outputdir, f'AO2D{outsuffix}.root')}")
            os.system("rm ao2ds_to_merge.txt")
            os.system("rm AO2D_*.root")
    else:
        os.system(f"mv AnalysisResults_0.root {os.path.join(outputdir, f'AnalysisResults{outsuffix}.root')}")
        os.system(f"mv AO2D_0.root {os.path.join(outputdir, f'AO2D{outsuffix}.root')}")


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("--infile", "-i", metavar="text",
                        default="files.txt", help="input file with list of paths (w/o file name)", required=True)
    PARSER.add_argument("--aod2s", "-a", action="store_true", default=False,
                        help="copy also AO2D.root files", required=False)
    PARSER.add_argument("--unmerged", "-u", action="store_true", default=False,
                        help="flag for unmerged train outputs", required=False)
    PARSER.add_argument("--outputdir", "-o", metavar="text", default=".",
                        help="output directory", required=False)
    PARSER.add_argument("--outputsuffix", "-s", metavar="text", default=".",
                        help="output suffix", required=False)
    ARGS = PARSER.parse_args()

    download(ARGS.infile, ARGS.aod2s, ARGS.unmerged, ARGS.outputdir, ARGS.outputsuffix)
