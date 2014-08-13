from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='PySeidon',
      version='v0.1',
      description='Description of the package goes here.',
      long_description=readme(),
      url='https://github.com/GrumpyNounours/PySeidon',
      author='Thomas Roc, Wesley Bowman, Jon Smith',
      author_email='thomas.roc@acadiau.ca,wesley.bowman23@gmail.com,lavieenroux20@gmail.com',
      maintainer='Thomas Roc',
      license='MIT',
      packages=['pyseidon'],
      #data_files=[('utide', ['utide/ut_constants.mat'])],
      zip_safe=False)
