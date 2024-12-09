from setuptools import setup, find_packages

setup(
    name="bifidotyper",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "bifidotyper": ["data/reference/*"],
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
)
