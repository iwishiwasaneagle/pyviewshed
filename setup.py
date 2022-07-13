import sys
from setuptools import find_packages
from skbuild import setup

if __name__ == "__main__":
    try:
        setup(
            cmake_install_dir="src/pyviewshed",
            cmake_args=["-DBUILD_TESTS=OFF"],
            packages=["pyviewshed"],
            package_dir={"": "src"},
        )
    except:  # noqa
        print(
            "\n\nAn error occurred while building the project, "
            "please ensure you have the most updated version of setuptools, "
            "setuptools_scm, skbuild, cmake, and pybind11 with:\n"
            "   pip install -U setuptools setuptools_scm scikit-build cmake pybind11"
            "\n\n",
            file=sys.stderr,
        )
        raise