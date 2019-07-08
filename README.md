This a trivial version of SPheRIO with a graphical interface. The code is written in python 2.7. 

To run the code, one need to install python 2.7, together with some libraries which can be installed by the following commands

pip install scipy
pip uninstall matplotlib
pip install matplotlib==2.0.2
pip install cython

In Debian or Ubuntu one should uses apt-get
apt-get install python-scipy
apt-get install python-matplotlib
apt-get install cython cython3

python setup.py build_ext --inplace

after this you should be able to run
python sph.py
and get a window to play with

version at
https://github.com/PhMota/sph/releases/tag/1.0.2
repository with the working version at
https://github.com/PhMota/sph

Ph. Mota
