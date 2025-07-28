import argparse
import sys

from spectral_wave_data import (
    SpectralWaveData,
    SwdFileCantOpenError,
    SwdFileBinaryError,
)


def main():
    """
    Show the metadata from an SWD file
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        prog="swd_meta", description="Display metadata for SWD files"
    )
    parser.add_argument("file_swd", help="SWD file to analyze")
    args = parser.parse_args()
    file_swd = args.file_swd

    try:
        swd = SpectralWaveData(file_swd, x0=0.0, y0=0.0, t0=0.0, beta=0.0)
    except SwdFileCantOpenError:
        print(f"Not able to open: {file_swd}")
        sys.exit(1)
    except SwdFileBinaryError:
        print(f"This SWD file don't have the correct binary convention: {file_swd}")
        sys.exit(2)
    except SwdFileBinaryError:
        print(f"This file don't look like a SWD-file: {file_swd}")
        sys.exit(3)

    def write_swd_tag(tag):
        print(f"{tag + ':':<8} {swd[tag]}")

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
