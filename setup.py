import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


from setuptools import setup

setuptools.setup(name='synumses-pkg-pabele',
      version='0.1',
      description='Symbolic and numerical simulation of semiconductor devices',
      url='',
      author='Peter Abele',
      author_email='ppabele@web.de',
      license='MIT',
      packages=setuptools.find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
