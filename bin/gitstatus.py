import os
import subprocess
import glob

basedir = os.path.join(os.environ['HOME'],'Documents')
gitdirs = ['Proposals', 'Talks', 'Teaching', 'TeXStuff', 'MATLAB', 'Notes',\
           'BragSheets', 'Reviews']
gitdirs = [os.path.join(basedir,d) for d in gitdirs]
gitdirs += glob.glob(os.path.join(basedir,'gitrepos/*'))

for d in gitdirs:
    if "gitrepos" in d:
        tmp = subprocess.run(["git", "remote", "-v"], cwd=d,
                             capture_output=True).stdout.decode()
        if "dsavransky" not in tmp:
            continue

    res = subprocess.run(["git", "status"], cwd=d, capture_output=True).stdout.decode()
    if ("Your branch is up to date" in res) and ("nothing to commit" in res) and\
            ("working tree clean" in res):
        continue
    print("{}\n".format(d))
    print(res)
    print("\n")

