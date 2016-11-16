from setuptools import setup, find_packages

setup(
    name='PySWaN',
    version='0.0',
    author='Gerben de Boer',
    author_email='gerben.deboer@vanoord.com',
    packages=find_packages(),
    description='Generic toolbox for spectral oceanwaves plus SWAN IO toolbox',
    install_requires=[
        'numpy',
        'scipy',
    ],
    #setup_requires=[
    #    'sphinx',
    #    'sphinx_rtd_theme'
    #],
)
