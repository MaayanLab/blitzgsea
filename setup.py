import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blitzgsea",
    version="0.1.",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="Package for fast calculation of GSEA similar to prerank.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maayanlab/blitzgsea",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_data={
        "prismx": ["data/*"]
    },
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy',
        'scikit-learn',
        'progress',
        'loess',
        'tqdm',
        'statsmodels'
    ],
    python_requires='>=3.6',
)