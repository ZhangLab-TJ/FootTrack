import os
import sys
import re
from setuptools import setup, Extension, dist, find_packages

#Test if numpy is installed
try:
	import numpy as np
except:
	#Else, fetch numpy if needed
	dist.Distribution().fetch_build_eggs(['numpy'])
	import numpy as np

#Add cython modules depending on the availability of cython
cmdclass = {}
try:
	from Cython.Distutils import build_ext
	cmdclass = {'build_ext': build_ext}
except ImportError:
	use_cython = False
else:
	use_cython = True

#To compile or not to compile
if use_cython:
	ext_modules = [Extension("foottrack.utils.ngs", ["foottrack/utils/ngs.pyx"], include_dirs=[np.get_include()]),
				Extension("foottrack.utils.sequences", ["foottrack/utils/sequences.pyx"], include_dirs=[np.get_include()]),
				Extension("foottrack.utils.signals", ["foottrack/utils/signals.pyx"], include_dirs=[np.get_include()])]

else:
	if os.path.exists("foottrack/utils/ngs.c"):
		ext_modules = [Extension("foottrack.utils.ngs", ["foottrack/utils/ngs.c"], include_dirs=[np.get_include()]),
						Extension("foottrack.utils.sequences", ["foottrack/utils/sequences.c"], include_dirs=[np.get_include()]),
						Extension("foottrack.utils.signals", ["foottrack/utils/signals.c"], include_dirs=[np.get_include()])] 
	else: #Only happens when installing directly from source
		sys.exit("Cython is needed to compile FootTrack.")

#Path of setup file to establish version
setupdir = os.path.abspath(os.path.dirname(__file__))

def find_version(init_file):
	version_file = open(init_file).read()
	version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
	if version_match:
		return version_match.group(1)
	else:
		raise RuntimeError("Unable to find version string.")

#Readme from git
def readme():
	with open('README.md') as f:
		return f.read()

setup(name='foottrack',
		version=find_version(os.path.join(setupdir, "foottrack", "__init__.py")),	#get version from __init__.py
		description='Transcription factor Occupancy prediction By Investigation of cFOOT-seq Signal',
		long_description=readme(),
		long_description_content_type='text/markdown',
		license='MIT',
		packages=find_packages(),
		entry_points={
			'console_scripts': ['FootTrack=foottrack.FootTrack:main']
		},
		ext_modules=ext_modules,
		cmdclass=cmdclass,
		setup_requires=["numpy"],
		install_requires=[
			'numpy',
			'scipy',
			'pysam',
			'pybedtools',
			'matplotlib>=2',
			'scikit-learn',
			'pandas',
			'pypdf2>=1.28.0',       #PdfFileMerger -> PdfMerger
			'xlsxwriter',
			'adjustText>=0.8',
			'pyBigWig>=0.3',		#to use values(numpy=True)
			'MOODS-python',
			'svist4get>=1.2.24',	#pin due to changes in API
			#'gimmemotifs',			#only used for MotifClust; now warns for this tool if it is not there
			'logomaker',
			'seaborn>=0.9.1',
			'boto3',
			'pyyaml>5.1',			#requirement for norns full_load
			'kneed'
		],
		scripts=["foottrack/scripts/filter_important_factors.py",
				 "foottrack/scripts/cluster_sites_by_overlap.py"],
		classifiers=[
			'License :: OSI Approved :: MIT License',
			'Intended Audience :: Science/Research',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			'Programming Language :: Python :: 3'
		],
		zip_safe=False	#gives cython import error if True
		)
