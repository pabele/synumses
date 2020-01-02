import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(name='synumses-pkg-pabele',
                 version='0.55',
                 description='Numerical simulation package for semiconductor devices',
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 url='',
                 author='Peter Abele',
                 author_email='ppabele@web.de',
                 license='MIT',
                 packages=setuptools.find_packages(),
                 classifiers=[
                     "Development Status :: 3 - Alpha",
                     "Programming Language :: Python :: 3.7",
                     "License :: OSI Approved :: MIT License",
                     "Topic :: Scientific/Engineering",
                     "Topic :: Scientific/Engineering :: Physics"
                 ],
                 python_requires='>=3.7',
)
