from setuptools import setup

setup(name='sherali_adams',
      version='0.1',
      description='run k rounds of SA hierarchy',
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
