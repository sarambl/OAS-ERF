import setuptools

setuptools.setup(
    name='OAS-ERF',
    version='1',
    packages=['oas_erf', 'oas_erf.util', 'oas_erf.util.Nd', 'oas_erf.util.Nd.sizedist_class_v2',
              'oas_erf.util.plot', 'oas_erf.util.imports', 'oas_erf.util.collocate', 'oas_erf.util.EBAS_data',
              'oas_erf.util.eusaar_data', 'oas_erf.util.slice_average', 'oas_erf.util.slice_average.avg_pkg',
              'oas_erf.data_info'],
    url='https://github.com/sarambl/OAS-ERF',
    license='MIT',
    author='Sara Blichner (sarambl)',
    author_email='sara.blichner@aces.su.se',
    description='Analysis code for acp publication'
)
