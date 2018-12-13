import setuptools

with open('README.md','r') as fh:
    long_description = fh.read()

setuptools.setup(
    name         = 'xover',
    version      = '0.0.0',
    author       = 'Matthew P. Humphreys',
    author_email = 'm.p.humphreys@cantab.net',
    description  = 'Cross-over analysis of hydrographic variables',
    url          = 'https://github.com/mvdh7/xover',
    packages     = setuptools.find_packages(),
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    classifiers = (
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',),)
