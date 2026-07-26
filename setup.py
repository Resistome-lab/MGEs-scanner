from setuptools import setup, find_packages

setup(
    name="mges-scanner",
    version="0.1.0",
    author="Mahmoud Faheem",
    description="Mobile Genetic Element Profiler & Mobility Risk Scoring Engine",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Resistome-lab/MGEs-scanner",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pandas>=1.3.0",
        "biopython>=1.79",
    ],
    entry_points={
        "console_scripts": [
            "mges-scanner=mges_scanner.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
