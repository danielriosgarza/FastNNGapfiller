from setuptools import setup, find_packages

setup(
    name='FastNNGapfiller',
    version='0.1.0',
    packages=find_packages(include=['FastNNGapfiller', 'FastNNGapfiller.*']),
    install_requires=[
        'cobra',
        'numpy>=1.14.5',
        'gurobipy'
    ]
)
