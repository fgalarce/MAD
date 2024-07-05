# =============================================================================
#  This file is part of the code MAD
#  Multy-physics for biomedicAl engineering and Data assimilation.
#  Copyright (C) 2017-2024,
#    
#     Felipe Galarce
#
#  MAD is free software; you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation; either version 2.1 of the License, or (at your option)
#  any later version.
#
#  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
#  more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with MAD. If not, see http://www.gnu.org/licenses/.
#  ============================================================================= 

# Download and install standalone Petsc/Mpich/Slepc last version
cd ../
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  echo "export MAD_ROOT=$(pwd)/source" >> ~/.bashrc
  echo "export PETSC_ROOT=$(pwd)/petsc" >> ~/.bashrc
  echo "alias mpirun=$(pwd)/petsc/build/bin/mpirun" >> ~/.bashrc
  echo "alias mad=cd $(pwd)/source" >> ~/.bashrc
  source ~/.bashrc
  sudo apt install build-essential
elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo "export MAD_ROOT=$(pwd)/source" >> ~/.zshrc
  echo "export PETSC_ROOT=$(pwd)/petsc" >> ~/.zshrc
  echo "alias mpirun=$(pwd)/petsc/build/bin/mpirun" >> ~/.zshrc
  echo "alias mad=cd $(pwd)/source" >> ~/.zshrc
  source ~/.zshrc
  xcode-select --install
fi
git clone -b release https://gitlab.com/petsc/petsc.git petsc
mkdir -p ./petsc/build
cd petsc
./configure --prefix=$PWD/build/ --download-openmpi=yes --download-cmake=yes --download-slepc=yes --download-superlu=yes --download-superlu_dist=yes --download-make=yes --with-debugging=no -download-f2cblaslapack=1
make all
make install

# Compile main software trunk
cd ../source
make

# MAD gui
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  sudo apt install glfw
elif [[ "$OSTYPE" == "darwin"* ]]; then
  brew reinstall glfw
fi

# Install required python modules
#python -m pip install numpy matplotlib vtk torch nibabel pykeops
