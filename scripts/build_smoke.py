import argparse
import os
import sys
from dataclasses import replace
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src_py"
sys.path.insert(0, str(SRC_DIR))

import database as database_module
from model_config import load_model_config
from utils import makeModel


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess one model configuration for a compile smoke test.")
    parser.add_argument("--model", required=True)
    parser.add_argument("--nlayers", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--h2-spin", type=int, choices=(0, 1), required=True)
    parser.add_argument("--layer-thickness", type=float)
    parser.add_argument("--thin-ice-approximation", type=int, choices=(0, 1))
    return parser.parse_args()


def main():
    args = parse_args()
    os.chdir(ROOT_DIR)
    makeModel(args.model)

    overrides = {
        "nlayers": args.nlayers,
        "h2_spin": bool(args.h2_spin),
        "multiprocessing": False,
    }
    if args.layer_thickness is not None:
        overrides["layer_thickness"] = args.layer_thickness
    if args.thin_ice_approximation is not None:
        overrides["thin_ice_approximation"] = bool(
            args.thin_ice_approximation)

    config = replace(load_model_config(), **overrides)
    database_module.load_model_config = lambda: config

    model = database_module.database()
    model.run(run=False)


if __name__ == "__main__":
    main()
