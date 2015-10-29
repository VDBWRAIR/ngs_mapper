from toolz.dicttoolz import valmap

env = jinja2.Environment(loader=jinja2.loaders.FileSystemLoader('.')
plate = env.get_template('Makefile')
cfg = yaml.load(open('config.yaml'))
string = plate.render(cfg.yaml)

# write it to the sample outdir 
ngs_data = cfg.yaml.pop('NGSDATA')
# Extract all the "default" values
conf = valmap(partial(valmap, lambda x: x['default']), cfg.yaml)


conf['args'] = {}
# have to change the argparse destination to match
# overwrite with args from command-line 
for argname, argval in args.__dict__:
    for cmd, option in conf:
    #this will work because the args are already over-written by the defaults in config.yaml if they were not suppplied. otherwise this would over-write with null values.
        if option == argname:
            conf[cmd][option] = argval
        else:
            #don't accidentally overwrite a command
            assert not (argname in conf.values())
            conf['args'][argname] = argval


# Jinja won't error out if a key doesn't exist if it doesn't have a dot call . . . 


#    def update_cmd(cmd, field): cfg.yaml[cmd][field] =  getattr(args, field)
