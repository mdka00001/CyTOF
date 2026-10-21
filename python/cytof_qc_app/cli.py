"""Command-line options for the Python QC pipeline."""
import argparse

from .pipeline import Pipeline


def main():
    parser = argparse.ArgumentParser(description="Scanpy/AnnData CyTOF QC")
    parser.add_argument("--input", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--intensities", default="intensities")
    parser.add_argument("--regionprops", default="regionprops")
    parser.add_argument("--images", default="img")
    parser.add_argument("--masks", default="masks")
    parser.add_argument("--transforms", default="arcsinh,arcsinh_zscore")
    parser.add_argument("--cofactor", type=float, default=1)
    parser.add_argument("--min-area", type=float, default=5)
    parser.add_argument("--pixel-size", type=float, default=1)
    parser.add_argument("--seed", type=int, default=220225)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--perplexity", type=float, default=30)
    parser.add_argument("--snr-cells", type=int, default=30000)
    parser.add_argument("--low-snr", type=float, default=2)
    parser.add_argument("--high-snr", type=float, default=10)
    parser.add_argument("--min-signal", type=float, default=2)
    parser.add_argument("--exclude", default="")
    parser.add_argument("--image-channels", default="")
    parser.add_argument("--skip-images", action="store_true")
    parser.add_argument("--skip-embeddings", action="store_true")
    args = parser.parse_args()
    try:
        Pipeline(args).run()
    except Exception as error:
        parser.exit(1, f"ERROR: {error}\n")

