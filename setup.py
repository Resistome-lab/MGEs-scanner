import io
from setuptools import setup, find_packages

setup(
    name="mges_scanner",
    version="0.1.0",
    description="Metagenomic Mobile Genetic Element (MGE) profiler and Mobility Risk Scorer",
    long_description=io.open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Mahmoud Faheem",
    packages=find_packages(),
    package_data={
        "mges_scanner": ["db/*"],
    },
    include_package_data=True,
    install_requires=[
        "pyhmmer>=0.10.0",
        "pyrodigal>=3.0.0",
    ],
    entry_points={
        "console_scripts": [
            "mges-scanner=mges_scanner.cli:main",
        ],
    },
    python_requires=">=3.8",
)
