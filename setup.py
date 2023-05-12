import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blitzgsea",
    version="1.3.31",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="Package for fast and accurate calculation of Gene Set Enrichment Analysis (GSEA) similar to prerank using gamma distribution approximation.",
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
        "blitzgsea": ["data/*"]
    },
    include_package_data=True,
    install_requires=[
        'pandas>=1.1.5',
        'numpy',
        'scikit-learn',
        'progress',
        'tqdm',
        'statsmodels',
        'mpmath'
    ],
    python_requires='>=3.6',
)
