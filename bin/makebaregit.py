import argparse
from miscpy.utils.gitutils import makebaregit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create new bare git repo"
    )
    parser.add_argument("name", nargs=1, type=str, help="Repo name (str).")
    args = parser.parse_args()

    makebaregit(args.name[0])


