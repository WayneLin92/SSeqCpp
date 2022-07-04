"""Generate a compact CMake project for AdamsE2"""
import argparse
import os
import subprocess

filenames_AdamsE2 = [
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

filenames_AdamsE2Prod = [
"./include/algebras/benchmark.h",
"./include/algebras/database.h",
"./include/algebras/myexception.h",
"./include/algebras/myio.h",
"./include/algebras/linalg.h",
"./include/algebras/steenrod.h",
"./include/algebras/utility.h",
"./src/linalg.cpp",
"./src/database.cpp",
"./src/myio.cpp",
"./src/steenrod.cpp",
"./src/utility.cpp",
"./AdamsE2/groebner_steenrod_const.h",
"./AdamsE2/groebner_steenrod_const.cpp",
"./AdamsE2/AdamsE2Prod.cpp",
]

filenames_AdamsE2Export = [
"./include/algebras/benchmark.h",
"./include/algebras/database.h",
"./include/algebras/dbAdamsSS.h",
"./include/algebras/myexception.h",
"./include/algebras/myio.h",
"./include/algebras/linalg.h",
"./include/algebras/steenrod.h",
"./include/algebras/utility.h",
"./include/algebras/algebras.h",
"./src/algebras.cpp",
"./src/linalg.cpp",
"./src/database.cpp",
"./src/dbAdamsSS.cpp",
"./src/myio.cpp",
"./src/steenrod.cpp",
"./src/utility.cpp",
"./AdamsE2/AdamsE2Export.cpp",
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

    if args.p == "AdamsE2":
        filenames = filenames_AdamsE2
    elif args.p == "AdamsE2Prod":
        filenames = filenames_AdamsE2Prod
    elif args.p == "AdamsE2Export":
        filenames = filenames_AdamsE2Export
    for fn in filenames:
        with open(fn) as file:
            content = file.read()
            content = content.replace("#include \"algebras/", "#include \"")
            content = content.replace("#include <sqlite3.h>", "#include \"sqlite3.h\"")
            content = content.replace("//#define MYDEPLOY", "#define MYDEPLOY")
            with open(os.path.join(args.o, os.path.basename(fn)), "w") as file_out:
                file_out.write(content)