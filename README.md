The majority of the documentation is locate at:
https://vdbpm.org/miseqpipeline/Wiki

# Installation

1. Clone

  ```
  git clone https://github.com/VDBWRAIR/miseqpipeline.git
  cd miseqpipeline
  ```

2. Ensure pre-requisites

  Read the REQUIREMENTS.md file

3. Setup virtualenv

  ```
  wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz#md5=f61cdd983d2c4e6aeabb70b1060d6f49 -O- | tar xzf -
  $HOME/bin/python virtualenv-1.11.6/virtualenv.py env 
  . env/bin/activate
  ```

4. Install into virtualenv

  ```
  python setup.py install
  ```
