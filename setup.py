from setuptools import setup

setup(
    name='CLAUDIO',
    version='v0.9.0',
    packages=['claudio', 'claudio.utils', 'claudio.module01', 'claudio.module01.src', 'claudio.module01.src.io',
              'claudio.module01.src.algorithm', 'claudio.module02',
              'claudio.module02.src_structure_search', 'claudio.module02.src_structure_search.io',
              'claudio.module02.src_structure_search.algorithm', 'claudio.module02.src_distance_reevaluation',
              'claudio.module02.src_distance_reevaluation.io', 'claudio.module02.src_distance_reevaluation.algorithm',
              'claudio.module03', 'claudio.module03.src', 'claudio.module03.src.io', 'claudio.module03.src.algorithm',
              'claudio.module04', 'claudio.module04.src', 'claudio.module04.src.io', 'claudio.module04.src.algorithm'],
    url='https://github.com/KohlbacherLab/CLAUDIO',
    license='MIT',
    author='Alexander Röhl',
    author_email='alexander.roehl@uni-tuebingen.de',
    description='Structural analysis, mapping, validation, visualization, and modeling of protein cross-links on '
                'protein and protein-protein interaction.',
    install_requires=['biopython==1.79', 'click==8.1.3', 'matplotlib==3.6.3', 'pandas==1.5.3', 'requests==2.28.2']
)
