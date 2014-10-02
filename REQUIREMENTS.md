# Requirements outside of python

- Python >=2.7.3
- samtools ==0.1.19[commit 96b5f2294ac0054230e88913c4983d548069ea4e]
- bwa ==0.7.6a
- Trimmomatic ==0.32

# Installation

  You need to ensure that you have a python virtualenv setup and activated prior to installing as these
  instructions assume that you do. The [README.md](README.md) has instructions on how to do that.

## Trimmomatic

  ```
  mkdir -p ~/src && pushd ~/src
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
  unzip -o Trimmo*.zip
  cp -R Trimmo*/ ${VIRTUAL_ENV}/lib/
  popd
  ```

  - Verify Trimmomatic

  ```
  [[ $(java -jar ${VIRTUAL_ENV}/lib/Trimmomatic-0.32/trimmomatic-0.32.jar 2>&1 | egrep 'SE|PE' | wc -l) == 2 ]] && echo "Trimmomatic installed" || echo "Trimmomatic not installed"
  ```

Back to [README.md](README.md)
