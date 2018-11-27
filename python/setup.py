from distutils.core import setup

with open('README.rst') as f:
    readme = f.read()

setup(
    name='iid-python',
    version='0.1',
    description='Pure python mathemetical and statistical functions.',
    url='https://github.com/drtconway/iid',
    long_description=readme,
    author='Thomas Conway',
    author_email='drtomc@gmail.com',
    packages=['iid', 'iid.tests'],
    classifiers=[
        'License :: OSI Approved :: Apache 2.0',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
    requires=['pyyaml']
)
