# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

# Do all the install steps
def provision_pipeline( config )
    # Clone the repo
    config.vm.provision "shell", privileged: false,
        inline: "echo 'Cloning'; [ -e ~/ngs_mapper ] && (cd ~/ngs_mapper && git pull && cd ..) || git clone /vagrant ~/ngs_mapper"

    # Setup config.yaml
    config.vm.provision "shell", privileged: false,
        inline: "echo 'Creating config.yaml'; mkdir -p ~/NGSDATA; sed 's|/path/to/NGSDATA|/home/vagrant/NGSDATA|' ~/ngs_mapper/ngs_mapper/config.yaml.default > ~/ngs_mapper/ngs_mapper/config.yaml"

    # Ensure that setup.py and vagrant-provision.sh are up to date in case they are not committed
    config.vm.provision "shell", privileged: false,
        inline: "cp -f /vagrant/setup.py ~/ngs_mapper; cp -f /vagrant/vagrant-provision.py ~/ngs_mapper;"

    # Ensure pipeline.log is owned by vagrant and 644
    config.vm.provision "shell", privileged: false,
        inline: "chown vagrant:vagrant /home/vagrant/ngs_mapper/pipeline.log"

    # Run the provisioning system-packages
    config.vm.provision "shell", privileged: false,
        inline: "echo 'Installing system packages'; cd ~/ngs_mapper; sudo python vagrant-provision.py --install-system-packages"

    # Run the provisioning pipeline install
    config.vm.provision "shell", privileged: false,
        inline: "echo 'Installing pipeline'; cd ~/ngs_mapper; python vagrant-provision.py --install-pipeline"
end

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
    config.vm.define "ubuntu" do |ubuntu|
        ubuntu.vm.box = "hashicorp/precise64"
        # Use a different mirror in case the ubuntu one is blocked...
        ubuntu.vm.provision "shell", privileged: true,
            inline: "sed -i -e 's/us.archive.ubuntu.com/ubuntu.osuosl.org/' -e 's|security.ubuntu.com|ubuntu.osuosl.org|' /etc/apt/sources.list"

        # Update
        ubuntu.vm.provision "shell", privileged: true,
            inline: "apt-get update"

        # Ensure git installed
        ubuntu.vm.provision "shell", privileged: true,
            inline: "apt-get install -y git"

        # Install and test pipeline
        provision_pipeline(ubuntu)
    end

    config.vm.define "centos" do |centos|
        centos.vm.box = "chef/centos-6.5"
        # Ensure git installed
        config.vm.provision "shell", privileged: true,
            inline: "yum install -y git"

        # Install setuptools so we can get argparse
        config.vm.provision "shell", privileged: true,
            inline: "yum install -y python-setuptools; easy_install argparse"

        # Install and test pipeline
        provision_pipeline(centos)
    end


    # Set the memory to 4GB
    config.vm.provider :virtualbox do |vb|
        vb.customize ["modifyvm", :id, "--memory", "4024"]
    end
end
