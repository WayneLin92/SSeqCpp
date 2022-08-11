"""Generate a compact CMake project for AdamsE2"""
import argparse
import os
import subprocess

filenames_Adams = [
"./Adams/AdamsRes.cpp",
"./Adams/AdamsResExport.cpp",
"./Adams/AdamsResProd.cpp",
"./Adams/groebner_steenrod_const.h",
"./Adams/groebner_steenrod.cpp",
"./Adams/groebner_steenrod.h",
"./Adams/main.h",
"./Adams/main.cpp",
"./include/algebras/algebras.h",
"./include/algebras/benchmark.h",
"./include/algebras/database.h",
"./include/algebras/dbAdamsSS.h",
"./include/algebras/linalg.h",
"./include/algebras/myexception.h",
"./include/algebras/myio.h",
"./include/algebras/steenrod.h",
"./include/algebras/utility.h",
"./src/algebras.cpp",
"./src/database.cpp",
"./src/linalg.cpp",
"./src/myio.cpp",
"./src/steenrod.cpp",
"./src/utility.cpp",
"./src/dbAdamsSS.cpp",
]

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('-o', help='output folder')
    parser.add_argument('-p', help='name of the project')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    if args.o is None or args.p is None:
        print("Need -o -p options")
        exit()

    if args.p == "Adams":
        filenames = filenames_Adams
    else:
        exit()
    for fn in filenames:
        with open(fn) as file:
            content = file.read()
            content = content.replace("#include \"algebras/", "#include \"")
            content = content.replace("#include <sqlite3.h>", "#include \"sqlite3.h\"")
            content = content.replace("//#define MYDEPLOY", "#define MYDEPLOY")
            with open(os.path.join(args.o, os.path.basename(fn)), "w") as file_out:
                file_out.write(content)