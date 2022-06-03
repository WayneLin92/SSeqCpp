"""Generate a compact CMake project for AdamsE2"""
import argparse
import os
import subprocess

filenames = [
"./include/algebras/benchmark.h",
"./include/algebras/database.h",
"./include/algebras/myexception.h",
"./include/algebras/myio.h",
"./include/algebras/steenrod.h",
"./include/algebras/utility.h",
"./src/database.cpp",
"./src/myio.cpp",
"./src/steenrod.cpp",
"./src/utility.cpp",
"./AdamsE2/groebner_steenrod.h",
"./AdamsE2/groebner_steenrod.cpp",
"./AdamsE2/AdamsE2.cpp",
]

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('-o', help='output folder')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    for fn in filenames:
        with open(fn) as file:
            content = file.read()
            content = content.replace("#include \"algebras/", "#include \"")
            content = content.replace("#include <sqlite3.h>", "#include \"sqlite3.h\"")
            with open(os.path.join(args.o, os.path.basename(fn)), "w") as file_out:
                file_out.write(content)