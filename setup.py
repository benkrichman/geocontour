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
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'shapely',
        'datascale'
    ],
)
