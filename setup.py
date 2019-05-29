import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ColourWorkflow",
    version="0.0.1",
    author="Gianluca Prandini, Gian-Marco Rignanese, Nicola Marzari",
    author_email="gianluca.prandini@epfl.ch",
    description="ColourWorkflow repository",
    long_description="This repository contains the workflow and the tools necessary in order to calculate from first principles the reflectivity and colour of metals within the independent-particle approximation (IPA)",
    long_description_content_type="text/markdown",
    url="",
    packages=['tools'], #setuptools.find_packages(),
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT and GNU licenses",
        "Operating System :: OS Independent",
    ],
)
