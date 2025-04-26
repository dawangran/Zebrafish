from setuptools import setup, find_packages

setup(
    name='Zebrafish',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'anndata',
        'scanpy',
        'omicverse',

    ],
)