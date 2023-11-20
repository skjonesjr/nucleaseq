from setuptools import setup
from nucleaseq.constants import VERSION
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np


if __name__ == '__main__':
    setup(
        name='nucleaseq',
        packages=['nucleaseq'],
        version=VERSION,
        include_dirs = [np.get_include()],
        entry_points={
          'console_scripts': [
              'nucleaseq = nucleaseq.main:main'
          ]
        },
        ext_modules=cythonize('nucleaseq/align_cy.pyx'),
        include_package_data=True,
        zip_safe=False,
        data_files=[
            ('resources', ['resources/base_logp.pkl']),
            ('params', [
                'params/cas12a_all_but_insweight_params.txt',
                'params/cas12a_pampwm_subpen_params.txt',
                'params/cas12a_pampwm_subpen_subtrans_params.txt',
                'params/cas12a_pampwm_subtrans_params.txt',
                'params/cas12a_single_effects_params.txt',
                'params/enh_all_but_insweight_params.txt',
                'params/enh_pampwm_subpen_params.txt',
                'params/enh_pampwm_subpen_subtrans_params.txt',
                'params/enh_pampwm_subtrans_params.txt',
                'params/enh_single_effects_params.txt',
                'params/hf1_all_but_insweight_params.txt',
                'params/hf1_pampwm_subpen_params.txt',
                'params/hf1_pampwm_subpen_subtrans_params.txt',
                'params/hf1_pampwm_subtrans_params.txt',
                'params/hf1_single_effects_params.txt',
                'params/hypa_all_but_insweight_params.txt',
                'params/hypa_pampwm_subpen_params.txt',
                'params/hypa_pampwm_subpen_subtrans_params.txt',
                'params/hypa_pampwm_subtrans_params.txt',
                'params/hypa_single_effects_params.txt',
                'params/wt_all_but_insweight_params.txt',
                'params/wt_pampwm_subpen_params.txt',
                'params/wt_pampwm_subpen_subtrans_params.txt',
                'params/wt_pampwm_subtrans_params.txt',
                'params/wt_single_effects_params.txt',
            ])],
        description='Process NucleaSeq data',
        url='http://www.finkelsteinlab.org',
        keywords=['DNA', 'NGS', 'bioinformatics', 'barcodes', 'Nuclease'],
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'License :: Freely Distributable',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
