import setuptools
from distutils.core import setup

setup(
    name='Localization',
    version='0.1.6',
    author='Kamal Shadi',
    author_email='kamal.shadi85@gmail.com',
    packages=['localization', 'localization.test'],
    scripts=['bin/sample.py'],
    url='https://github.com/kamalshadi/Localization',
    license='LICENSE.txt',
    description='Multilateration and Triangulation.',
    long_description=open('README.md').read(),
)
