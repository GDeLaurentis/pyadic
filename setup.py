from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with (this_directory / "pentagon_functions" / "version.py").open() as f:
    version = f.read().split(" = '")[1].split("'\n")[0]


setup(
    name='pyadic',
    version=version,
    license='GNU General Public License v3.0',
    description='p-Adic numbers and finite fields in Python',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    url='https://github.com/GDeLaurentis/pyadic',
    download_url=f'https://github.com/GDeLaurentis/pyadic/archive/v{version}.tar.gz',
    project_urls={
        'Documentation': 'https://gdelaurentis.github.io/pyadic/',
        'Issues': 'https://github.com/GDeLaurentis/pyadic/issues',
    },
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
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
)
