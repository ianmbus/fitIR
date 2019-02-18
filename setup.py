from setuptools import setup, find_packages


# with open('README.rst') as f:
#     readme = f.read()
# 
# with open('LICENSE') as f:
#     license = f.read()

setup(
    name='fitIR',
    version='0.1.0',
    description='IR SED fitting using MCMC',
    # long_description=readme,
    author='Stephen Wilkins',
    author_email='s.wilkins@sussex.ac.uk',
    url='https://github.com/stephenmwilkins/fitIR',
    # license=license,
    packages=find_packages(exclude=('examples', 'docs'))
)