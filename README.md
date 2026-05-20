# NGS-myassembling

Using this addon: https://jschoeberl.github.io/iFEM/FEM/myassembling.html

follow the installation instruction from project https://github.com/TUWien-ASC/NGS-myfe.

    git clone https://github.com/TUWien-ASC/NGS-myassembling.git
    cd NGS-myassembling
    mkdir build
    cd build
    cmake ..
    make -j4 install


* Install to `conda` base environment:
  ```
  cd NGS-myassembling
  mkdir build
  cd build

  cmake .. \
  -DPython3_EXECUTABLE=/Users/myh/anaconda3/bin/python3.11 \
  -DCMAKE_INSTALL_PREFIX=/Users/myh/anaconda3/lib/python3.11/site-packages

  make -j4 install
  ```
* Check the source of myassembling
  ```
  import myassembling

  print(myassembling)
  print(myassembling.__file__)
  print(dir(myassembling))
  ```
