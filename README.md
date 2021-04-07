# pytheas
Software to identify and map RNA modifications via LC-MS/MS


Make sure to have at least 30GB of dedicated HDD space (especially important if running Linux on a Virtual Machine), otherwise the wxpython installation will fail
$ sudo apt-get update
$ sudo apt-get install python3
#$ sudo apt install -y make gcc libgtk-3-dev freeglut3 freeglut3-dev python3-#gst-1.0 libglib2.0-dev ubuntu-restricted-extras libgstreamer-plugins-#base1.0-dev
$ sudo apt install -y python3-pip
$ wget  -P ~ "https://extras.wxpython.org/wxPython4/extras/linux/gtk3/ubuntu-20.04/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl
$ pip3 install ~/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl
$ pip3 install -r requirements.txt
