from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(name='sherali_adams',
      version='0.2',
      description='run k rounds of SA hierarchy',
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords='integer programming mip sherali-adams optimization',
      url='http://github.com/knavely/sherali_adams',
      author='Matthew Drescher',
      author_email='knavely@gmail.com',
      license='Apache 2.0',
      test_suite='nose.collector',
      tests_require=['nose'],
      packages=['sherali_adams'],
      install_requires=['numpy'],
      zip_safe=False)
