import os
import os.path
import subprocess
import shutil
import inspect
import miscpy


def makebaregit(name):
    """
    Create new bare git repo
    """
    basedir = os.path.join(os.sep, "data", "git")
    gitdir = name + ".git"
    dd = os.path.join(basedir, gitdir)
    os.mkdir(dd)
    res = subprocess.run(["git", "init", "--bare"], cwd=dd, capture_output=True)
    assert res.returncode == 0

    res2 = subprocess.run(
        ["chown", "-R", "git.git", gitdir], cwd=basedir, capture_output=True
    )
    assert res2.returncode == 0


def makenewgit(path, gitignore=None):
    """
    Create new git repo with either default or provided gitignore
    path - full path to repo
    gitignore - full path to gitignore file (or None for default)
    """

    path = os.path.abspath(
        os.path.normpath(os.path.expandvars(os.path.expanduser(path)))
    )
    os.mkdir(path)

    if gitignore is None:
        gitignore = os.path.normpath(
            os.path.join(os.path.dirname(inspect.getfile(miscpy)), "..", ".gitignore")
        )

    assert os.path.exists(gitignore)

    _ = shutil.copy(gitignore, path)

    res1 = subprocess.run(["git", "init"], cwd=path, capture_output=True)
    assert res1.returncode == 0

    res2 = subprocess.run(["git", "add", "."], cwd=path, capture_output=True)
    assert res2.returncode == 0
    res3 = subprocess.run(
        ["git", "commit", "-m 'initial commit'"], cwd=path, capture_output=True
    )
    assert res3.returncode == 0
