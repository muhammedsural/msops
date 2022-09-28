#from setuptools import find_packages, setup
from distutils.core import setup
import os
import sys


def setup_package():
    src_path = os.path.dirname(os.path.abspath(__file__))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    
    from setuptools import setup

    try:
      setup(
    name='msops',
    packages=["msops"],
    version='0.4.1',
    description='Python library for main shock after shock analysis with openseespy ',
    author='Muhammed Sural',
    license='MIT',
    author_email = "muhammedsural@gmail.com",
    url = 'https://github.com/muhammedsural/msops',
    download_url = 'https://github.com/muhammedsural/msops/archive/refs/tags/0.4.1.tar.gz',
    install_requires=[            
          "EzGM==1.6.5.5",
          "numpy==1.22.4",
          "openseespy==3.4.0.1",
          "openseespywin==3.4.0.1.1",
          "opsvis==1.0.12",
          "pandas==1.2.5",
          "scipy==1.6.3"
      ],
    classifiers=['Development Status :: 3 - Alpha',    
                 'Intended Audience :: Developers',      
                 'Topic :: Software Development :: Build Tools',
                 'License :: OSI Approved :: MIT License',   
                 'Programming Language :: Python :: 3.9'
                ]
  )
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return


if __name__ == '__main__':
    setup_package()

    



