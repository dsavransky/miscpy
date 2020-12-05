import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="admissions",
    version="1.0.0",
    author="Dmitry Savransky",
    author_email="ds264@cornell.edu",
    license='MIT',
    description="Miscellaneous Python utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dsavransky/miscpy",
    packages=["miscpy"],
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
