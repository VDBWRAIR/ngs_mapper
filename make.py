from toolz.dicttoolz import valmap
import jinja2
import os, sys
import yaml
from functools import partial
# write it to the sample outdir
default_args= {
    'prefix' : 947,
    'reference' : '947.ref.fasta'
}


#@typecheck(int, str, str)
def make_template(template, dest, args=default_args, config='config.yaml'):
    env = jinja2.Environment(loader=jinja2.loaders.FileSystemLoader(os.path.dirname(template)))
    plate = env.get_template(template)
    settings = get_settings(config, args)
    print settings
    string = plate.render(settings)
    with open(dest, 'w') as out:
        out.write(string)
    return dest

# Extract all the "default" values
def get_settings(confpath, args):
    cfg = yaml.load(open(confpath))
    _ = cfg.pop('NGSDATA')
    conf = valmap(partial(valmap, lambda x: x['default']), cfg)
    conf['args'] = {}
    # have to change the argparse destination to match
    # overwrite with args from command-line
    for argname, argval in args.__dict__.items():
        for cmd, option in conf.items():
        #this will work because the args are already over-written by the defaults in config.yaml if they were not suppplied. otherwise this would over-write with null values.
            if option == argname:
                conf[cmd][option] = argval
            else:
                #don't accidentally overwrite a command
                assert not (argname in conf.values())
                conf['args'][argname] = argval
    return conf

usage = "<template> <dest>"

if __name__ == '__main__':
    try:
        make_template(*sys.arv[1:])
    except Exception, e:
        print "Usage: " + usage


# Jinja won't error out if a key doesn't exist if it doesn't have a dot call . . .


#    def update_cmd(cmd, field): cfg.yaml[cmd][field] =  getattr(args, field)
