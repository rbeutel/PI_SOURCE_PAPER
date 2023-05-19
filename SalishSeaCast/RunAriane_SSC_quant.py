"""
Loop over dates, setting up a series of forcing links, running Ariane starting particles over multiple days
and saving the statistic results"""

import argparse
import arrow
import os
import subprocess

TARGET_TMPL = 'SalishSea_1h_{:06d}_grid_{:s}.nc' # :03d - d part converts a number into decimal format, 03 part states that is has to be 3 digits (if not lead with 0s); :s - assuming for strings
FILENAME_TMPL = 'SalishSea_1h_{:%Y%m%d}_{:%Y%m%d}_grid_{:s}.nc'
SUBDIR_TMPL = '{:%d%b%y}'
NEW_SUBDIR_TMPL = 'for_{:%d%b%y}'


def make_links(rundate, runlength):
    dir = '/results2/SalishSea/nowcast-green.201905/'
    tardir = 'Links'
    for grid in ['T', 'U', 'V', 'W']:
        for fileno in range(runlength):
            target = TARGET_TMPL.format(fileno+1, grid)
            date = rundate.shift(days=+fileno).datetime #rundate is an 'arrow' object (date thingy), .shift moves it forward by the amount of days in brakets, .datetime makes datetime object
            print (date, dir)
            link = FILENAME_TMPL.format(date, date, grid)
            subdir = SUBDIR_TMPL.format(date).lower()
            try:
                os.unlink(os.path.join(tardir, target))
            except FileNotFoundError:
                pass
            os.symlink(os.path.join(dir, subdir, link),
                       os.path.join(tardir, target))


def run_ariane():
    with open('babypoo', 'wt') as stdout:
        with open('errpoo', 'wt') as stderr:
            subprocess.run(
                "/ocean/rbeutel/MOAD/ariane-2.3.0_03/bin/ariane",
                stdout=stdout, stderr=stderr,
                universal_newlines=True)


def rename_results(rundate=arrow.utcnow(), subdir='',
                   nday=1, labeltype='date'):
    finaldir = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/'
    filelist = ['ariane_memory.log', 'ariane_statistics_quantitative.nc',
                'final_pos.txt', 'init_pos.txt', 'output',
                'ariane_positions_quantitative.nc', 'final.sav',
                'init.sav', 'stats.txt']
    for filename in filelist:
        if labeltype == 'date':
            newdir = NEW_SUBDIR_TMPL.format(rundate.datetime).lower()
        else:
            newdir = ''
        print(os.path.join(finaldir, subdir, newdir, filename))
        if not os.path.exists(os.path.join(finaldir, subdir, newdir)):
                os.makedirs(os.path.join(finaldir, subdir, newdir))
        os.rename(filename, os.path.join(finaldir, subdir, newdir, filename))


def main(args):
    initialrundate = arrow.get(args.initialrundate, 'YYYY-MM-DD')
    for nday in range(args.numberofdays):
        rundate = initialrundate.shift(days=+nday)
        print ('Start', rundate)
        if args.forback == 'forward':
            startfile = rundate
        elif args.forback == 'backward':
            startfile = rundate.shift(days=-(args.runlength+args.releaseparticledays-1))
        make_links(startfile, args.runlength+args.releaseparticledays)
        print ('Startfile', startfile)
        run_ariane()
        rename_results(rundate=rundate, subdir=args.runtype)
        print ('End', rundate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'initialrundate', help='Date to start from as YYYY-MM-DD', type=str)
    parser.add_argument(
        'numberofdays', help='Number of different dates to do', type=int)
    parser.add_argument(
        'runlength', help='Number of days to track particles', type=int)
    parser.add_argument(
        'releaseparticledays', help='Number of days over which to release particles for one run', type=int)
    parser.add_argument('forback', help='Run forward or backward', type=str)
    parser.add_argument('runtype',
                        help='FluxesSouth FluxesNorth BackSouth BackNorth',
                        type=str)
    args = parser.parse_args()
    main(args)
