import setuptools as st

with open('README.md','r') as ldf:
    longdesc=ldf.read()

st.setup(
    name='geocontour',
    version='1.3.1',
    license='MIT',
    url='https://github.com/benkrichman/geocontour',
    description='Utilities for masking, contour tracing, and geocontour construction with gridded geographic data',
    long_description=longdesc,
    long_description_content_type="text/markdown",
    author='Benjamin Krichman',
    author_email='benkrichman@gmail.com',
    packages=st.find_packages(),
    include_package_data=True,
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
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Scientific/Engineering :: Visualization"
    ]
)


