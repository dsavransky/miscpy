import setuptools
import os, re

with open("README.md", "r") as fh:
    long_description = fh.read()


with open(os.path.join("miscpy", "__init__.py"), "r") as f:
    version_file = f.read()

version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)

if version_match:
    version_string = version_match.group(1)
else:
    raise RuntimeError("Unable to find version string.")


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="miscpy",
    version=version_string,
    author="Dmitry Savransky",
    author_email="ds264@cornell.edu",
    license="MIT",
    description="Miscellaneous Python utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dsavransky/miscpy",
    packages=setuptools.find_packages(),
    install_requires=["sympy", "matplotlib", "numpy", "scipy"],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
