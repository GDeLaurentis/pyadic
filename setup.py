from setuptools import setup, find_packages


setup(
    name='pyadic',
    version='0.1.0',
    author='Giuseppe De Laurentis',
    author_email='g.dl@hotmail.it',
    description='p-Adic numbers and finite fields in Python',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['numpy', 'sympy'],
)
