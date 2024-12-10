from setuptools import setup, find_packages
import re

try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "See project details on GitHub: https://github.com/Bennibraun/bifidotyper"

setup(
    name="bifidotyper",
    version=re.search(r'__version__ = "(.*?)"', open("src/bifidotyper/__init__.py", "r").read()).group(1),
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "bifidotyper": ["data/reference/*", "data/bin/*"],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'bifidotyper=bifidotyper.cli:main',
        ],
    },
	install_requires=[
        'numpy',
        'pandas',
        'seaborn',
        'matplotlib',
		'tqdm',
		'scikit-learn',
		'natsort',
		'palettable',
    ],
	url='https://github.com/Bennibraun/bifidotyper',
	description='A bioinformatics tool for analyzing Bifidobacteria in sequencing data.',
	long_description=long_description,
	long_description_content_type="text/markdown",
    author="Ben Braun",
    author_email="braun.ben-1@colorado.edu",
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
	python_requires=">=3.7",
)