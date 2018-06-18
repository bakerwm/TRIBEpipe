from setuptools import setup, find_packages

def readme():
    with open('README.md', 'r') as f:
        return f.read()

def license():
    with open('LICENSE', 'r') as f:
        return f.read()

setup(
    name='TRIBEpipe',
    version='0.1',
    description='TRIBE analysis pipeline',
    long_description=readme(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Text Processing :: Linguistic',
        ],
    keywords='TRIBE',
    url='http://github.com/',
    author='Ming Wang',
    author_email='wangm08@hotmail.com',
    license=license(),
    packages=find_packages(), # or ['goldclip'],
    install_requires=[
        'pysam >= 0.11.1',
        'numpy >= 1.13.1',
        'pandas >= 0.18.1',
        'pybedtools >= 0.7.10',
        ],
    # data_files=[('bin', ['bin/*.sh'])],
    #scripts=['bin/funniest-joke'],
    entry_points={
        'console_scripts': ['TRIBEpipe=TRIBEpipe.TRIBEpipe:main'],
    },
    test_suite='nose.collector',
    tests_require=['nose'],
    include_package_data=True,
    zip_safe=False
    )
