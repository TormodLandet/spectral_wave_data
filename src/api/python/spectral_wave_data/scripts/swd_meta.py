import argparse
import sys

from spectral_wave_data import (
    SpectralWaveData,
    SwdFileCantOpenError,
    SwdFileBinaryError,
    SwdFileDataError,
    SwdInputValueError,
    version_full,
)


def main():
    """
    Show the metadata from an SWD file
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(prog="swd_meta", description="Display metadata for SWD files")
    parser.add_argument("file_swd", nargs="*", help="SWD file to analyze")
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    args = parser.parse_args()

    if args.version:
        print(version_full)
        return

    if not args.file_swd:
        print("ERROR: No SWD file provided. Please specify a file.")
        parser.print_usage()
        sys.exit(1)

    for file_swd in args.file_swd:
        swd_meta_main(file_swd)


def swd_meta_main(file_swd):
    """
    Process a single SWD file and print its metadata
    """
    try:
        swd = SpectralWaveData(file_swd)
    except SwdFileCantOpenError as e:
        print(f"Not able to open SWD file ({e}): {file_swd}")
        sys.exit(1)
    except SwdFileBinaryError as e:
        print(f"This file does not have the correct binary convention ({e}): {file_swd}")
        sys.exit(2)
    except SwdFileDataError as e:
        print(f"This file does not look like a SWD-file ({e}): {file_swd}")
        sys.exit(3)
    except SwdInputValueError as e:
        print(f"This file does not accept current input values ({e}): {file_swd}")
        sys.exit(3)
    except Exception as e:
        print(f"This file has unexpected error ({e}): {file_swd}")
        sys.exit(3)

    def write_swd_tag(tag):
        print(f"{tag + ':':<8} {swd[tag]}")

    print(f"Metadata for SWD file: {file_swd}\n{'=' * 40}\n")
    write_swd_tag("version")
    write_swd_tag("prog")
    write_swd_tag("date")
    write_swd_tag("fmt")
    write_swd_tag("shp")
    write_swd_tag("amp")
    write_swd_tag("tmax")
    write_swd_tag("dt")
    write_swd_tag("nsteps")
    write_swd_tag("nstrip")
    write_swd_tag("order")
    write_swd_tag("d")

    shp = swd["shp"]
    if shp in [1, 2, 3]:
        # Long-crested seas
        write_swd_tag("n")
        if shp == 3:
            write_swd_tag("nh")
        write_swd_tag("sizex")
        write_swd_tag("lmax")
        write_swd_tag("lmin")
        write_swd_tag("dk")

    if shp in [4, 5]:
        # Short-crested seas
        write_swd_tag("nx")
        write_swd_tag("ny")
        write_swd_tag("sizex")
        write_swd_tag("sizey")
        write_swd_tag("lmax")
        write_swd_tag("lmin")
        write_swd_tag("dkx")
        write_swd_tag("dky")

    if shp in [6]:
        # Airy waves
        write_swd_tag("n")

    write_swd_tag("cid")


if __name__ == "__main__":
    main()
