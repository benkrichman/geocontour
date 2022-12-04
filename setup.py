import setuptools as st

with open('README.md','r') as ldf:
    longdesc=ldf.read()

st.setup(
    name='geocontour',
    version='1.0.0',
    license='MIT',
    url='https://github.com/benkrichman/geocontour',
    description='Functions for applying contour tracing to gridded data',
    long_description=longdesc,
    long_description_content_type="text/markdown",
    author='Benjamin Krichman',
    author_email='benkrichman@gmail.com',
    packages=st.find_packages(),
    install_package_data=True,
    package_data={'geocontour':['geocontour/data/*']},
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'shapely',
        'datascale'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: Hydrology"
    ]
)
