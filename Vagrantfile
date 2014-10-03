# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
    #config.vm.box = "ubuntu/trusty64"
    config.vm.box = "chef/centos-6.5"

    # The url from where the 'config.vm.box' box will be fetched if it
    # doesn't already exist on the user's system.
    # config.vm.box_url = "http://domain.com/path/to/above.box"

    # Set the memory to 4GB
    config.vm.provider :virtualbox do |vb|
        vb.customize ["modifyvm", :id, "--memory", "4024"]
    end
end
