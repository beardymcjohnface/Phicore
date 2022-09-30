"""
Entrypoint for Strand Swaps
"""

import sys
import os
import subprocess
import yaml
from shutil import copyfile
from time import localtime, strftime

import click


"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
- Set a different default system config file in copy_config() if you need to"""


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def print_version():
    with open(snake_base('strand_swaps.VERSION'), 'r') as f:
        version = f.readline()
    click.echo('\n' + 'Strand Swaps version ' + version + '\n')


def msg(err_message):
    tstamp = strftime('[%Y:%m:%d %H:%M:%S] ', localtime())
    click.echo(tstamp + err_message)


def msg_box(splash, errmsg=None):
    msg('-' * (len(splash) + 4))
    msg(f'| {splash} |')
    msg(('-' * (len(splash) + 4)))
    if errmsg:
        click.echo('\n' + errmsg)


def copy_config(local_config, system_config=snake_base(os.path.join('config', 'config.yaml'))):
    if not os.path.isfile(local_config):
        msg(f'Copying system default config to {local_config}')
        copyfile(system_config, local_config)
    else:
        msg(f'Config file {local_config} already exists. Using existing config file.')


def read_config(file):
    with open(file, 'r') as stream:
        _config = yaml.safe_load(stream)
    return _config


def write_config(_config, file):
    msg(f'Writing runtime config file to {file}')
    with open(file, 'w') as stream:
        yaml.dump(_config, stream)


class OrderedCommands(click.Group):
    """Preserve the order of subcommands when printing --help"""
    def list_commands(self, ctx: click.Context):
        return list(self.commands)


"""RUN A SNAKEFILE
Hopefully you shouldn't need to tweak this function at all.
- You must provide a Snakefile, all else is optional
- Highly recommend supplying a configfile and the default snakemake args"""


def run_snakemake(configfile=None, snakefile_path=None, merge_config=None, profile=None, threads=1, use_conda=False,
                  conda_frontend=None, conda_prefix=None, outdir=None, snake_default_args=None, snake_extra=None):
    """Run a Snakefile"""
    snake_command = ['snakemake', '-s', snakefile_path]

    # if using a configfile
    if configfile:
        # copy sys default config if needed
        copy_config(configfile)

        # read the config
        snake_config = read_config(configfile)

        # merge in command line config if provided
        if merge_config:
            snake_config.update(merge_config)

        # create runtime config file for Snakemake execution
        if outdir:
            runtime_config = os.path.join(outdir, 'strand_swaps.config.yaml')
            if not os.path.exists(os.path.normpath(outdir)):
                os.makedirs(os.path.normpath(outdir))
        else:
            runtime_config = 'strand_swaps.config.yaml'
        write_config(snake_config, runtime_config)
        snake_command += ['--configfile', runtime_config]

        # display the runtime configuration
        msg_box('Runtime config', errmsg=yaml.dump(snake_config, Dumper=yaml.Dumper))

    # add --profile [profile]
    if profile:
        snake_command += ['--profile', profile]

    # add threads
    snake_command += ['--cores', threads]

    # add conda args if using conda
    if use_conda:
        snake_command += ['--use-conda']
        if conda_frontend:
            snake_command += ['--conda-frontend', conda_frontend]
        if conda_prefix:
            snake_command += ['--conda-prefix', conda_prefix]

    # add snakemake default args
    if snake_default_args:
        snake_command += snake_default_args

    # add any additional snakemake commands
    if snake_extra:
        snake_command += list(snake_extra)

    # Run Snakemake!!!
    snake_command = ' '.join(str(s) for s in snake_command)
    msg_box('Snakemake command', errmsg=snake_command)
    if not subprocess.run(snake_command, shell=True).returncode == 0:
        msg('Error: Snakemake failed')
        sys.exit(1)
    else:
        msg('Snakemake finished successfully')
    return 0


"""COMMON OPTIONS
Any click options that will be used in more than one subcommand can be defined here"""


def common_options(func):
    """Common options decorator for use with click commands."""
    options = [
        click.option('--output', help='Output directory', type=click.Path(),
                     default='strand_swaps.out', show_default=True),
        click.option('--configfile', default='config.yaml', help='Custom config file', show_default=True),
        click.option('--threads', help='Number of threads to use', default=8, show_default=True),
        click.option('--profile', help='Run on cluster with a Snakemake profile', type=str),
        click.option('--use-conda/--no-use-conda', default=True, help='Use conda for Snakemake rules',
                     show_default=True),
        click.option('--conda-frontend',
                     type=click.Choice(['mamba', 'conda'], case_sensitive=True),
                     default='mamba', help='Specify Conda frontend', show_default=True),
        click.option('--conda-prefix', default=snake_base(os.path.join('workflow', 'conda')),
                     help='Custom conda env directory', type=click.Path(), show_default=False),
        click.option('--snake-default', multiple=True,
                     default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs', '--notemp'],
                     help="Customise Snakemake runtime args", show_default=True),
    ]
    for option in reversed(options):
        func = option(func)
    return func


"""COMMAND LINE ARGS
Customise your launcher!
"""


@click.group(cls=OrderedCommands)
def cli():
    """For more options, run:
    strand_swaps command --help"""
    pass


"""MAIN SCRIPT
Add or edit any command line args related to running the main script here,
customise the epilog etc.
"""


EPILOG = """
\b
RUN EXAMPLES:
  Required:             strand_swaps run --input [file]
  Specify threads:      strand_swaps run ... --threads [threads]
  Run on cluster:       strand_swaps run ... --profile [profile]
  Disable conda:        strand_swaps run ... --no-use-conda 
  Change defaults:      strand_swaps run ... --snake-default="-k --nolock"
  Add Snakemake args:   strand_swaps run ... --dry-run --keep-going --touch
  Specify targets:      strand_swaps run ... all print_targets
  Available targets:
        all             Run everything (default)
        print_targets   List available targets
"""


@click.command(epilog=EPILOG, context_settings={"ignore_unknown_options": True})
@click.option('--input', '_input', help='Input file/directory', type=str, required=True)
@click.argument('snake_args', nargs=-1)
@common_options
def run(_input, configfile, output, threads, profile, use_conda, conda_frontend, conda_prefix, snake_default,
        snake_args, **kwargs):
    """Run Strand Swaps"""
    # Config to add or update in configfile
    merge_config = {
        'input': _input,
        'output': output,}

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'Snakefile')),   # Full path to Snakefile
        configfile=configfile,
        outdir=output,
        merge_config=merge_config,
        threads=threads,
        profile=profile,
        use_conda=use_conda,
        conda_frontend=conda_frontend,
        conda_prefix=conda_prefix,
        snake_default_args=snake_default,
        snake_extra=snake_args,
    )


"""SUBCOMMAND EXAMPLE
Copy the system default config file to working directory.
"""


@click.command()
@click.option('--configfile', default='config.yaml', help='Custom config file', show_default=True)
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


cli.add_command(run)
cli.add_command(config)


def main():
    print_version()
    cli()


if __name__ == '__main__':
    main()
