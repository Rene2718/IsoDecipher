from setuptools import setup, find_packages

setup(
    name="IsoDecipher",
    version="0.1",
    packages=find_packages(include=["IsoDecipher", "IsoDecipher.*"]),
    python_requires=">=3.10",
)