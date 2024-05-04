from setuptools import setup, find_packages

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='pyadic',
    version='0.2.1',
    license='GNU General Public License v3.0',
    description='p-Adic numbers and finite fields in Python',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/pyadic',
    download_url='https://github.com/GDeLaurentis/pyadic/archive/v0.2.0.tar.gz',
    keywords=['p-adic numbers', 'finite-fields', 'rational reconstruction'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'sympy',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.9',
    ],
)
