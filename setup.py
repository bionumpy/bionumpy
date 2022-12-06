#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy>=1.20',
                'npstructures>=0.2.3']
# 'npstructures @ git+https://github.com/knutdrand/npstructures.git']

test_requirements = ['pytest>=3', ]

setup(
    author="Knut Rand",
    author_email='knutdrand@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Library for working with biological sequence data as numpy arrays.",
    entry_points={
        'console_scripts': [
            'bionumpy=bionumpy.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='bionumpy',
    name='bionumpy',
    packages=find_packages(include=['bionumpy', 'bionumpy.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/knutdrand/bionumpy',
    version='0.2.8',
    zip_safe=False,
)

# python -m build
# twine upload dist/*
