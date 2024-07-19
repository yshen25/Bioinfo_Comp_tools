from setuptools import setup, find_packages

setup(
    name='Bioinfo_Comp_tools',  # Replace with your package name
    version='1.0',
    packages=find_packages(),
    install_requires=[],  # List your package dependencies here
    author='Yue (Shawn) Shen',
    author_email='yshen25@alum.utk.edu',
    description='Simple but useful python scripts for Bioinformatics and Computational Biology tasks',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yshen25/Bioinfo_Comp_tools',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GPL-3.0 license',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.9',
)