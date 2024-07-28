import setuptools
with open("README.md", "r") as readme:
    long_description = readme.read()

setuptools.setup(
    name="Leftover Reads Inspector",
    version="1.0.0",
    author="Thomas Collins",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Linux"
    ],
    install_requires=[
        "subprocess",
        "csv",
        "sys",
        "itertools",
        "pandas",
        "argparse",
        "traceback",
        "matplotlib.pyplot",
        "matplotlib",
        "difflib",
        "math",
        "collections",
        "numpy",
        "seaborn",
        "random"
        ],
    packages=setuptools.find_packages(),
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "LeftoverReadsInspector = src.driver:main"
        ]
    }
)
