import argparse
from miscpy.utils.gitutils import makenewgit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create new git repo"
    )
    parser.add_argument("path", nargs=1, type=str, help="Repo path (str).")
    args = parser.parse_args()

    makenewgit(args.path[0])
