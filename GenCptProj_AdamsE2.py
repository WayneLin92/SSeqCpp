"""Generate a compact CMake project for AdamsE2"""
import argparse
import os
import subprocess

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('-o', help='output folder')
    parser.add_argument('-o', help='output folder')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions