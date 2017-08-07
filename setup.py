from setuptools import setup, find_packages

# get the version information from the relevant file
version = {}
with open('./cubric_mrs/_version.py') as f:
    exec(f.read(), version)


setup(
        name='cubric-mrs',
        version=version['__version__'],
        packages=find_packages(),
        url='https://github.com/bennyrowland/cubric_mrs.git',
        license='MIT',
        author='bennyrowland',
        author_email='bennyrowland@mac.com',
        description='a collection of processing scripts for handling MRS data at CUBRIC',
        entry_points={
            "console_scripts": [
                "mrs_mega = cubric_mrs.megapress:megapress_script",
                "mrs_press = cubric_mrs.press:press_script"
            ]
        },
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'Topic :: Scientific/Engineering :: Physics',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
        ],
        install_requires=['suspect', 'numpy', 'Pillow', 'pyx', 'nibabel'],
        test_requires=['pytest']
)
