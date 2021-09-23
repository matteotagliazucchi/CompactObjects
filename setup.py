from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

install_requires = [
    "numpy",
    "matplotlib",
    "scipy",
    "astropy"
]

setup(
    name='compactobjects',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    version='0.1.0',
    description='A package for studying astrophysical compact objects.',
    author='Matteo Tagliazucchi',
    author_email="matteo.tagliazucchi@protonmail.com",
    url="https://github.com/matteotagliazucchi/compactobjects",
    license='None',
    install_requires=install_requires,
    classifiers=['Development Status :: 1 - Beta',
                 'Intended Audience :: Developers',
                 'Topic :: Software Development :: Libraries',
                 'Programming Language :: Python :: 3.8',        
                 "Operating System :: OS Independent"],
)


