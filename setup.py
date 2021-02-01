import setuptools

with open("DESCRIPTION.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(name='synumses-pkg-pabele',
                 version='0.8.3',
                 description='Numerical solver for semiconductor devices in Python',
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 url='',
                 author='Peter Abele',
                 author_email='ppabele@web.de',
                 license='GNU GPLv3',
                 packages=setuptools.find_packages(),
                 classifiers=[
                     "Development Status :: 4 - Beta",
                     "Programming Language :: Python :: 3.7",
                     "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
                     "Topic :: Scientific/Engineering",
                     "Topic :: Scientific/Engineering :: Mathematics",
                     "Topic :: Scientific/Engineering :: Physics",
                     "Topic :: Scientific/Engineering :: Visualization",
                     "Topic :: Software Development :: Libraries :: Python Modules",
                 ],
                 python_requires='>=3.7',
                 install_requires=[
                     'sympy',
                     'scipy',
                     'numpy',
                 ]

)
