from setuptools import setup, find_packages


# def readme():
#     with open('README.rst') as f:
#         return f.read()


setup(name='ngsCAT2',
      version='0.1',
      description='Next-Generation Sequencing Capture Assessment Tool 2',
      long_description='',
      classifiers=[
        'Development Status :: 1 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Enrichment analysis',
      ],
      keywords='targeted sequencing  analysis ',
      url='http://github.com/storborg/funniest',
      author='Alejandro García Sánchez',
      author_email='alegarsan2@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
                        'certifi==2018.10.15',
                        'chardet==3.0.4',
                        'decorator==4.3.0',
                        'idna==2.7',
                        'ipython-genutils==0.2.0',
                        'jsonschema==2.6.0',
                        'jupyter-core==4.4.0',
                        'nbformat==4.4.0',
                        'numexpr==2.6.9',
                        'numpy==1.22.0',
                        'pandas==0.23.4',
                        'plotly==3.3.0',
                        'pysam==0.15.1',
                        'python-dateutil==2.7.5',
                        'pytz==2018.6',
                        'requests==2.20.0',
                        'retrying==1.3.3',
                        'scipy==1.1.0',
                        'six==1.11.0',
                        'tables==3.4.4',
                        'toolz==0.9.0',
                        'traitlets==4.3.2',
                        'urllib3==1.24',
                        'xlwt==1.3.0',
      ],
      #test_suite='',
      #tests_require=['nose', 'nose-cover3'],
      entry_points={
          'console_scripts': ['ngscat2=ngscat2.main:main'],
      },
      include_package_data=True,
      zip_safe=False)
